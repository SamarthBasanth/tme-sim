/**
 * Procedural, GBM-referenced tumor point clouds + vessel curves.
 *
 * ⚠ FALLBACK ONLY, AND ITS PROVENANCE IS FROM THE WRONG PRODUCT.
 *
 * This does NOT render while `public/models/tumor.glb` is present:
 * `TumorModel` resolves to `GLTFTumor` and the real exported Xenium section
 * ships instead. It survives as the no-asset fallback path.
 *
 * The references below are to MRI segmentation and 3-D-printed SURGICAL
 * PLANNING literature — the product this page used to describe, before it
 * was corrected to a drug-transport simulator. They describe how a
 * glioblastoma reads on contrast-enhanced MRI, which is a genuine and
 * correctly-cited basis for a *procedural stand-in*, but it is not what
 * Vivum reconstructs: Vivum builds tissue from spatial transcriptomics,
 * cell by cell, and never touches an MRI.
 *
 * So: the citations are accurate about what they claim, and irrelevant to
 * the shipped product. If this placeholder is ever restored to the render
 * path it needs the same retarget the rest of the scene got (see
 * ConcentrationField / DrugFront / FieldContours), and the BraTS
 * NCR/ET/ED label convention below should go with it.
 *
 * Grounded in how glioblastoma actually reads on contrast-enhanced MRI and
 * in gross/3D-printed surgical-planning specimens (BraTS challenge
 * documentation; see PMC11945730 "Advancing Precision: A Comprehensive
 * Review of MRI Segmentation Datasets from BraTS Challenges", and the
 * multi-material 3D-print surgical planning literature, e.g. PubMed
 * 34656248): a GBM mass is irregular and lobulated (not a smooth sphere),
 * and is famous for *ring enhancement* — a necrotic, non-enhancing core
 * (hypointense on T1Gd) surrounded by a brighter, irregular enhancing rim
 * (BraTS label convention: NCR=necrotic core, ET=enhancing tumor,
 * ED=peritumoral edema).
 *
 * Rendered as a point cloud instead of a displaced mesh — matching how the
 * real product (app/web/js/gl/viewer.js + shaders.js) actually renders a
 * tumour specimen: individual per-cell points, not a smooth surface. Each
 * point samples a category (0 = core tissue, 1 = necrotic, 2 =
 * enhancing/margin rim) so a single draw can be coloured like the real
 * viewer's category-attribute channel. Points fill a *volume*, not just a
 * shell: density is biased toward the interior and thinned toward an
 * irregular, noise-jittered outer boundary, so the point cloud itself
 * carries the tumour's lobulated silhouette rather than relying on an
 * invisible mesh boundary underneath it.
 *
 * All noise is evaluated once per point inside `useMemo` (never per frame)
 * using the deterministic CPU noise in `lib/noise.ts`.
 */
import * as THREE from "three";
import { fbm3, warpedFbm3 } from "@/lib/noise";
import { mulberry32 } from "@/lib/prng";

function smoothstep(edge0: number, edge1: number, x: number): number {
  const t = THREE.MathUtils.clamp((x - edge0) / (edge1 - edge0), 0, 1);
  return t * t * (3 - 2 * t);
}

export interface NecroticPit {
  center: THREE.Vector3;
  angularRadius: number;
  seed: number;
}

/** Two off-center, differently sized cavities so the result reads as organic, not symmetric. */
export const DEFAULT_PITS: readonly NecroticPit[] = [
  { center: new THREE.Vector3(0.55, 0.4, 0.62).normalize(), angularRadius: 0.46, seed: 11.3 },
  { center: new THREE.Vector3(-0.48, -0.32, 0.36).normalize(), angularRadius: 0.32, seed: 47.1 },
  { center: new THREE.Vector3(-0.2, 0.62, -0.5).normalize(), angularRadius: 0.26, seed: 83.7 },
];

function pitInfluence(dir: THREE.Vector3, pits: readonly NecroticPit[]): number {
  let best = 0;
  for (const pit of pits) {
    const d = THREE.MathUtils.clamp(dir.dot(pit.center), -1, 1);
    const angle = Math.acos(d);
    const jitter = fbm3(dir.x * 2.4 + pit.seed, dir.y * 2.4, dir.z * 2.4 + pit.seed, 3) * 0.16;
    const edge = pit.angularRadius + jitter;
    const infl = 1 - smoothstep(edge * 0.35, edge, angle);
    if (infl > best) best = infl;
  }
  return best;
}

/** Uniform random direction on the unit sphere, seeded PRNG in, unit Vector3 out. */
function randomDirection(rand: () => number, out: THREE.Vector3): THREE.Vector3 {
  const z = rand() * 2 - 1;
  const theta = rand() * Math.PI * 2;
  const s = Math.sqrt(Math.max(0, 1 - z * z));
  return out.set(s * Math.cos(theta), s * Math.sin(theta), z);
}

/** Category codes shared by the point-cloud shader (see pointCloudShader.ts). */
export const CATEGORY_CORE = 0;
export const CATEGORY_NECROTIC = 1;
export const CATEGORY_ENHANCING = 2;

export interface PointCloudData {
  positions: Float32Array;
  /** Per-point category, 0..2 — see CATEGORY_* above. */
  categories: Float32Array;
  /** Per-point 0..1 seeded random value, used for size/opacity variance. */
  jitter: Float32Array;
}

export interface MassCloudOptions {
  count?: number;
  seed?: number;
  lobeAmp?: number;
  fineAmp?: number;
  pitDepth?: number;
  pits?: readonly NecroticPit[];
}

/**
 * The tumor mass as a point cloud: category 0 (core tissue) fills most of
 * the volume, category 1 (necrotic) clusters near the pit centres close to
 * the surface — the point-cloud analogue of the old mesh's necrotic
 * dimples — and category 2 (enhancing/margin rim) is denser right at the
 * pit boundary and at the mass's own outer margin, standing in for BraTS
 * ring enhancement. The outer boundary is rejection-thinned with a
 * noise-jittered threshold so the silhouette itself is ragged, not a hard
 * sphere clipped by an invisible surface.
 */
export function sampleTumorMassPointCloud(radius: number, opts: MassCloudOptions = {}): PointCloudData {
  const {
    count = 9000,
    seed = 5081,
    lobeAmp = 0.26,
    fineAmp = 0.06,
    pitDepth = 0.17,
    pits = DEFAULT_PITS,
  } = opts;

  const rand = mulberry32(seed);
  const positions = new Float32Array(count * 3);
  const categories = new Float32Array(count);
  const jitter = new Float32Array(count);
  const dir = new THREE.Vector3();

  const maxAttempts = count * 8;
  let written = 0;
  let attempts = 0;

  while (written < count && attempts < maxAttempts) {
    attempts++;
    randomDirection(rand, dir);

    const coarse = warpedFbm3(dir.x * 1.1, dir.y * 1.1, dir.z * 1.1, 0.7, 4);
    const fine = fbm3(dir.x * 5.5, dir.y * 5.5, dir.z * 5.5, 3);
    const pit = pitInfluence(dir, pits);
    const localRadius = radius * (1 + coarse * lobeAmp + fine * fineAmp - pit * pitDepth);
    if (localRadius <= 0) continue;

    // Roughly-uniform volumetric fill (cbrt keeps number density constant
    // with r), then a noise-jittered rejection band near the boundary so
    // the margin itself reads as irregular rather than a clipped sphere.
    const t = Math.cbrt(rand());
    const marginNoise = fbm3(dir.x * 4 + 12.7, dir.y * 4, dir.z * 4 + 12.7, 3);
    const marginThreshold = 0.86 + marginNoise * 0.08;
    if (t > marginThreshold) {
      const acceptProb = 1 - (t - marginThreshold) / (1 - marginThreshold);
      if (rand() > acceptProb * acceptProb) continue;
    }

    const r = localRadius * t;
    positions[written * 3] = dir.x * r;
    positions[written * 3 + 1] = dir.y * r;
    positions[written * 3 + 2] = dir.z * r;

    const necroticAmount = smoothstep(0.4, 0.8, pit) * smoothstep(1, 0.7, t);
    const rimBand = smoothstep(0.08, 0.42, pit) * smoothstep(0.95, 0.42, pit);
    const marginBand = smoothstep(marginThreshold - 0.18, marginThreshold, t);

    let category = CATEGORY_CORE;
    const roll = rand();
    if (necroticAmount > 0.3 && roll < necroticAmount * 0.95) {
      category = CATEGORY_NECROTIC;
    } else if ((rimBand > 0.35 || marginBand > 0.35) && roll < Math.max(rimBand, marginBand) * 0.9) {
      category = CATEGORY_ENHANCING;
    }

    categories[written] = category;
    jitter[written] = rand();
    written++;
  }

  if (written === count) return { positions, categories, jitter };
  return {
    positions: positions.subarray(0, written * 3),
    categories: categories.subarray(0, written),
    jitter: jitter.subarray(0, written),
  };
}

export interface HaloCloudOptions {
  count?: number;
  seed?: number;
  /** Domain-warp frequency + offset for the bounding-volume noise field. */
  freq?: number;
  offset?: number;
  ampNoise?: number;
  /** [inner, outer] fraction of the local radius points are sampled between. */
  shell?: [number, number];
  shellGamma?: number;
}

/**
 * A sparse, shell-biased point cloud for a diffuse halo volume (edema, or
 * the faint outer margin shell) — deliberately lower density and a
 * different radial fill than the mass cloud so it reads as a diffuse
 * swelling wrapping the mass, not a second solid cell population.
 */
function sampleHaloPointCloud(radius: number, opts: HaloCloudOptions = {}): PointCloudData {
  const {
    count = 2200,
    seed = 6203,
    freq = 0.75,
    offset = 40,
    ampNoise = 0.16,
    shell = [0.55, 1],
    shellGamma = 0.6,
  } = opts;

  const rand = mulberry32(seed);
  const positions = new Float32Array(count * 3);
  const categories = new Float32Array(count);
  const jitter = new Float32Array(count);
  const dir = new THREE.Vector3();
  const [shellLo, shellHi] = shell;

  for (let i = 0; i < count; i++) {
    randomDirection(rand, dir);
    const n = warpedFbm3(dir.x * freq + offset, dir.y * freq, dir.z * freq + offset, 0.8, 3);
    const localRadius = radius * (1 + n * ampNoise);
    const t = shellLo + Math.pow(rand(), shellGamma) * (shellHi - shellLo);
    const r = localRadius * t;

    positions[i * 3] = dir.x * r;
    positions[i * 3 + 1] = dir.y * r;
    positions[i * 3 + 2] = dir.z * r;
    categories[i] = CATEGORY_CORE;
    jitter[i] = rand();
  }

  return { positions, categories, jitter };
}

/** Sparser, differently-coloured point cloud standing in for peritumoral edema. */
export function sampleEdemaPointCloud(radius: number, count = 2200): PointCloudData {
  return sampleHaloPointCloud(radius, {
    count,
    seed: 6203,
    freq: 0.75,
    offset: 40,
    ampNoise: 0.16,
    shell: [0.55, 1],
    shellGamma: 0.6,
  });
}

/** Very sparse, faint point cloud standing in for surrounding-tissue margin. */
export function sampleMarginPointCloud(radius: number, count = 1200): PointCloudData {
  return sampleHaloPointCloud(radius, {
    count,
    seed: 9012,
    freq: 0.6,
    offset: 80,
    ampNoise: 0.08,
    shell: [0.82, 1],
    shellGamma: 0.75,
  });
}

export interface VesselCurveOptions {
  count?: number;
  seed?: number;
  radius?: number;
}

/**
 * A handful of thin, branching vessel-like curves hugging just outside the
 * mass surface — a cheap stand-in for surface vascularity (real GBM is
 * famously hypervascular). Built from a seeded PRNG so it is deterministic
 * across renders.
 */
export function createVesselCurves(opts: VesselCurveOptions = {}): THREE.CatmullRomCurve3[] {
  const { count = 5, seed = 2024, radius = 0.66 } = opts;
  const rand = mulberry32(seed);
  const curves: THREE.CatmullRomCurve3[] = [];

  for (let i = 0; i < count; i++) {
    const theta0 = rand() * Math.PI * 2;
    const phi0 = Math.acos(2 * rand() - 1);
    const start = new THREE.Vector3(
      Math.sin(phi0) * Math.cos(theta0),
      Math.sin(phi0) * Math.sin(theta0),
      Math.cos(phi0),
    );

    const points: THREE.Vector3[] = [];
    const segments = 4 + Math.floor(rand() * 3);
    let dir = start.clone();
    for (let s = 0; s <= segments; s++) {
      const wobble = new THREE.Vector3(
        (rand() - 0.5) * 0.22,
        (rand() - 0.5) * 0.22,
        (rand() - 0.5) * 0.22,
      );
      dir = dir.clone().add(wobble).normalize();
      // Deliberately well inside the nominal mass radius: the mass surface
      // itself is heavily noise-displaced (see sampleTumorMassPointCloud)
      // and can dip well below `radius` locally, so a curve pinned just
      // outside `radius` would poke through empty space wherever the
      // surface dips — read as a stray floating line rather than a vessel
      // under the tissue. Sitting further inside keeps it plausibly "just
      // under the surface" everywhere instead.
      const r = radius * (0.62 + rand() * 0.12);
      points.push(dir.clone().multiplyScalar(r));
    }
    curves.push(new THREE.CatmullRomCurve3(points));
  }

  return curves;
}

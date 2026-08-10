/**
 * The SOLVED interstitial fluid pressure field, as iso-contours over the same
 * measured section the hero point cloud shows.
 *
 * `public/models/ifp_contours.json` is written by
 * `experiments/export_landing_asset.py`, which builds it from the SAME
 * reconstruction and the SAME coordinate transform as `tumor.glb` — one
 * `reconstruct_real()` call, one scale-and-centre, both outputs. That shared
 * transform is the reason these contours register with the cells rather than
 * merely sitting near them. The vertices arrive already in the page's
 * standing-up scene frame (the exporter lays the section in the glTF XZ plane
 * at y == 0; GLTFTumor turns it -90 degrees about X, which sends the asset's z
 * to world y), so nothing here transforms anything — it only reshapes the JSON
 * into the tuple arrays drei's `<Line>` wants.
 *
 * WHAT IS REAL, AND THE ONE DISPLAY CHOICE
 * ----------------------------------------
 * The field is a genuine steady-state solve —
 * `tmesim.physics.pressure.solve_ifp`, a screened-Poisson Darcy model with
 * leaky-vessel Starling sources driven by the section's MEASURED endothelial
 * cell positions and lymphatic drainage suppressed inside the tumour. Nothing
 * about it is authored.
 *
 * One display choice is applied between the solver and these vertices, and it
 * is stated on the page: the contoured quantity is the solved field's LOCAL
 * MEAN over a 250 um window (`coarseGrainUm`). That is forced by resolution,
 * not by taste. With a discrete vessel source the field's screening length is
 * of order 60 um, so the raw solution is thousands of per-vessel haloes, each
 * well under one screen pixel when 8.1 mm of tissue is fitted into ~2.6 scene
 * units — drawn literally it is speckle, not a field. Because it is a mean,
 * the level values here (~0.2-0.65 mmHg) are far below the raw field's peak
 * (`rawMaxMmHg`, ~12.9 mmHg); the page must never quote one as the other.
 *
 * The structure that survives the averaging is real: one dominant lobe in the
 * lower-left of the section, sitting on the section's measured
 * endothelial-density maximum — which `shaders/viridis.ts` derived
 * independently as `ENTRY_POINTS[0]`, from the GLB, with no knowledge of this
 * solve. High pressure where the vessels actually are is the model's central
 * mechanism, and here it falls out of the measured geometry.
 */
import raw from "@/public/models/ifp_contours.json";

/** One iso-level of the solved field, with every closed/open path at it. */
export interface IfpLevel {
  /** Iso-level in mmHg, of the coarse-grained field. */
  valueMmHg: number;
  /** `valueMmHg / displayedMaxMmHg` — also this level's position on the ramp. */
  fractionOfMax: number;
  /**
   * `viridis(fractionOfMax)` under the product's own ramp table
   * (`app/web/js/gl/ramps.js`, `STOPS.viridis`) — the table
   * `shaders/viridis.ts` subsamples for the concentration field, so a value
   * is coloured here exactly as the Studio would colour it. Baked by the
   * exporter rather than sampled in TypeScript so there is one ramp, not two.
   */
  colorHex: string;
  /** Polylines in scene units, ready for drei's `<Line points>`, z left to the caller. */
  paths: Array<Array<[number, number, number]>>;
}

interface RawLevel {
  fraction_of_max: number;
  value_mmhg: number;
  color_srgb_hex: string;
  paths: number[][][];
}

const doc = raw as unknown as {
  quantity: string;
  units: string;
  solver: string;
  displayed_field_max_mmhg: number;
  raw_field_max_mmhg: number;
  display_coarse_grain_um: number;
  um_per_scene_unit: number;
  levels: RawLevel[];
};

/**
 * Outermost (lowest pressure) level first, matching the JSON's stated
 * `level_order`. Consumers that thin the set for a reduced detail budget
 * therefore slice from the low end.
 */
export const IFP_LEVELS: readonly IfpLevel[] = doc.levels.map((level) => ({
  valueMmHg: level.value_mmhg,
  fractionOfMax: level.fraction_of_max,
  colorHex: level.color_srgb_hex,
  paths: level.paths.map((path) => path.map(([x, y]) => [x, y, 0] as [number, number, number])),
}));

/** Provenance the page is allowed to state, kept next to the numbers it describes. */
export const IFP_FIELD = {
  quantity: doc.quantity,
  units: doc.units,
  solver: doc.solver,
  /** Peak of the CONTOURED (coarse-grained) field, i.e. the ramp's top. */
  displayedMaxMmHg: doc.displayed_field_max_mmhg,
  /** Peak of the RAW solved field, which is ~20x higher. Not what is drawn. */
  rawMaxMmHg: doc.raw_field_max_mmhg,
  coarseGrainUm: doc.display_coarse_grain_um,
  umPerSceneUnit: doc.um_per_scene_unit,
} as const;

/** Lowest and highest drawn iso-level, for the on-canvas legend. */
export const IFP_LEVEL_RANGE_MMHG: readonly [number, number] = [
  IFP_LEVELS.length ? IFP_LEVELS[0].valueMmHg : 0,
  IFP_LEVELS.length ? IFP_LEVELS[IFP_LEVELS.length - 1].valueMmHg : 0,
];

/**
 * `"0.22–0.65 mmHg"` — the range actually drawn, in units, for the annotation.
 * Two decimals because the levels are ~0.2-0.7 mmHg; one would collapse two of
 * them onto the same printed value.
 */
export const IFP_LEVEL_RANGE_LABEL = `${IFP_LEVEL_RANGE_MMHG[0].toFixed(2)}–${IFP_LEVEL_RANGE_MMHG[1].toFixed(2)} ${IFP_FIELD.units}`;

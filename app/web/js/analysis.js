// analysis.js — measurement over the payloads that are already on the client.
//
// Everything in here is a deterministic function of the TMSG cell array and
// the TMSF field grid the backend served. Nothing is fetched, nothing is
// modelled, nothing is invented: a number that cannot be derived from those
// two buffers is returned as null with a reason, exactly as the API does.
//
// Estimator notes, because these numbers get read as if they were exact:
//   * Quantiles come from a 1024-bin histogram CDF with linear interpolation
//     inside the containing bin, NOT from a sort. Worst-case error is one bin
//     width, i.e. (max-min)/1024 of the range. This is stated in the UI.
//   * Field statistics over a region are computed on the served grid nodes
//     inside the region, so they inherit the run's dx. A region narrower than
//     one dx contains no node and reports null, never an interpolated value.
//   * Per-cell exposure is the field bilinearly interpolated AT the measured
//     cell centroid. That is an interpolation of a computed field, so it is
//     computed — never measured.

const HIST_BINS = 1024;

/* ====================================================== histogram ======= */

/**
 * @param {Float32Array} values
 * @returns {{bins:Uint32Array, n:number, min:number, max:number, total:number,
 *            nonFinite:number, peak:number}}
 */
export function buildHistogram(values, min, max, n = HIST_BINS) {
  const bins = new Uint32Array(n);
  let lo = min, hi = max;
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) {
    // Degenerate range (constant field, or an all-NaN grid). Recover the real
    // extremes rather than fabricating a 0..1 axis.
    lo = Infinity; hi = -Infinity;
    for (let i = 0; i < values.length; i++) {
      const v = values[i];
      if (!Number.isFinite(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    if (!Number.isFinite(lo)) return { bins, n, min: 0, max: 0, total: 0, nonFinite: values.length, peak: 0 };
    if (hi === lo) hi = lo + 1e-12;
  }
  const s = n / (hi - lo);
  let total = 0, nonFinite = 0;
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (!Number.isFinite(v)) { nonFinite++; continue; }
    let b = ((v - lo) * s) | 0;
    if (b < 0) b = 0; else if (b >= n) b = n - 1;
    bins[b]++;
    total++;
  }
  let peak = 0;
  for (let i = 0; i < n; i++) if (bins[i] > peak) peak = bins[i];
  return { bins, n, min: lo, max: hi, total, nonFinite, peak };
}

/** Quantile in [0,1] from a histogram CDF. Returns null for an empty histogram. */
export function quantile(hist, q) {
  if (!hist || !hist.total) return null;
  const want = q * hist.total;
  let acc = 0;
  const w = (hist.max - hist.min) / hist.n;
  for (let i = 0; i < hist.n; i++) {
    const c = hist.bins[i];
    if (acc + c >= want) {
      const f = c > 0 ? (want - acc) / c : 0;
      return hist.min + (i + f) * w;
    }
    acc += c;
  }
  return hist.max;
}

/** Re-bin a histogram down to `m` display columns (max, not mean — a single
 *  tall spike must not be averaged out of existence at display width). */
export function rebin(hist, m) {
  const out = new Float64Array(m);
  if (!hist) return out;
  const per = hist.n / m;
  for (let j = 0; j < m; j++) {
    const a = Math.floor(j * per), b = Math.min(hist.n, Math.floor((j + 1) * per));
    let s = 0;
    for (let i = a; i < Math.max(a + 1, b); i++) s += hist.bins[i] || 0;
    out[j] = s;
  }
  return out;
}

/* ==================================================== line profile ====== */

/**
 * Sample a field along a line in um. Returns distances from p0 and values;
 * a sample outside the grid is NaN and counted, never silently dropped.
 */
export function sampleLine(viewer, x0, y0, x1, y1, n) {
  const dx = x1 - x0, dy = y1 - y0;
  const len = Math.hypot(dx, dy);
  const N = Math.max(2, Math.min(1024, n | 0));
  const d = new Float64Array(N);
  const v = new Float64Array(N);
  let outside = 0;
  let min = Infinity, max = -Infinity, sum = 0, cnt = 0;
  for (let i = 0; i < N; i++) {
    const t = i / (N - 1);
    d[i] = len * t;
    const s = viewer.sampleField(x0 + dx * t, y0 + dy * t);
    if (s === null || !Number.isFinite(s)) { v[i] = NaN; outside++; continue; }
    v[i] = s;
    if (s < min) min = s;
    if (s > max) max = s;
    sum += s; cnt++;
  }
  return {
    d, v, n: N, lengthUm: len, outside,
    min: cnt ? min : null, max: cnt ? max : null,
    mean: cnt ? sum / cnt : null,
    stepUm: N > 1 ? len / (N - 1) : 0,
    valid: cnt,
  };
}

/* ============================================================ ROI ======= */

/** Grid-node index window covering a um rectangle, using the same node
 *  mapping as Viewer.sampleField (node i sits at i*extent/(cols-1)). */
function nodeWindow(rect, field, ex, ey) {
  const cx = (field.cols - 1) / ex, cy = (field.rows - 1) / ey;
  return {
    i0: Math.max(0, Math.ceil(rect.x0 * cx)),
    i1: Math.min(field.cols - 1, Math.floor(rect.x1 * cx)),
    j0: Math.max(0, Math.ceil(rect.y0 * cy)),
    j1: Math.min(field.rows - 1, Math.floor(rect.y1 * cy)),
  };
}

/**
 * Field statistics over the served grid nodes inside a um rectangle.
 * @returns {null|{n, mean, sd, min, max, p10, p50, p90, nonFinite}}
 */
export function rectFieldStats(field, ex, ey, rect) {
  if (!field) return null;
  const { i0, i1, j0, j1 } = nodeWindow(rect, field, ex, ey);
  if (i1 < i0 || j1 < j0) return { n: 0, reason: 'region is narrower than one grid step (dx)' };
  let n = 0, nonFinite = 0, sum = 0, sum2 = 0, min = Infinity, max = -Infinity;
  for (let j = j0; j <= j1; j++) {
    const row = j * field.cols;
    for (let i = i0; i <= i1; i++) {
      const v = field.values[row + i];
      if (!Number.isFinite(v)) { nonFinite++; continue; }
      n++; sum += v; sum2 += v * v;
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  if (!n) return { n: 0, nonFinite, reason: 'no finite grid node inside the region' };
  const mean = sum / n;
  const varr = Math.max(0, sum2 / n - mean * mean);
  // second pass only for quantiles, and only over the region
  const hist = buildHistogram(subGrid(field, i0, i1, j0, j1), min, max, 512);
  return {
    n, nonFinite, mean, sd: Math.sqrt(varr), min, max,
    p10: quantile(hist, 0.10), p50: quantile(hist, 0.50), p90: quantile(hist, 0.90),
  };
}

function subGrid(field, i0, i1, j0, j1) {
  const w = i1 - i0 + 1, h = j1 - j0 + 1;
  const out = new Float32Array(w * h);
  for (let j = 0; j < h; j++) {
    const src = (j0 + j) * field.cols + i0;
    out.set(field.values.subarray(src, src + w), j * w);
  }
  return out;
}

/**
 * Cell composition and per-cell field exposure inside a um rectangle.
 * `types` counts are MEASURED (segmented phenotype calls); exposure is the
 * field interpolated at those positions, so it is COMPUTED.
 */
export function rectCellStats(viewer, rect, nTypes = 5) {
  const g = viewer.geom;
  if (!g || !viewer.index) return null;
  const counts = new Int32Array(nTypes);
  let total = 0;
  const exposure = [];
  const hasField = !!viewer.field;
  let sum = 0, sum2 = 0, min = Infinity, max = -Infinity, nExp = 0;
  const tumorSum = { s: 0, n: 0 };
  viewer.index.forEachInRect(rect.x0, rect.y0, rect.x1, rect.y1, (i) => {
    const t = g.types[i];
    if (t < nTypes) counts[t]++;
    total++;
    if (!hasField) return;
    const v = viewer.sampleField(g.xs[i], g.ys[i]);
    if (v === null || !Number.isFinite(v)) return;
    nExp++; sum += v; sum2 += v * v;
    if (v < min) min = v;
    if (v > max) max = v;
    if (t === 0) { tumorSum.s += v; tumorSum.n++; }
    if (exposure.length < 400000) exposure.push(v);
  });
  let exp = null;
  const values = Float32Array.from(exposure);
  if (nExp) {
    const mean = sum / nExp;
    const hist = buildHistogram(values, min, max, 512);
    exp = {
      n: nExp, mean, sd: Math.sqrt(Math.max(0, sum2 / nExp - mean * mean)),
      min, max,
      p10: quantile(hist, 0.10), p50: quantile(hist, 0.50), p90: quantile(hist, 0.90),
      tumorMean: tumorSum.n ? tumorSum.s / tumorSum.n : null,
      tumorN: tumorSum.n,
    };
  }
  // `values` is handed back so the caller can histogram the same samples
  // against a different range (the field-wide axis, for the overlay) without
  // walking the region a second time.
  return { counts, total, exposure: exp, values };
}

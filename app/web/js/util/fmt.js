// fmt.js — number formatting.
//
// Rule of the house: a number is never rendered without a unit, and a missing
// number is rendered as an explicit em-dash placeholder, never as 0, never as
// "N/A" hidden in prose. `null` from the backend means "we do not know", and
// the UI says so.

export const NULL_GLYPH = '—'; // em dash

/** Fixed-decimal with thousands separators; null-safe. */
export function num(v, dp = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  if (!Number.isFinite(v)) return v > 0 ? '∞' : '-∞';
  return v.toLocaleString('en-US', { minimumFractionDigits: dp, maximumFractionDigits: dp });
}

/** Integer with thousands separators. */
export function int(v) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  return Math.round(v).toLocaleString('en-US');
}

/** Scientific notation, e.g. 1.20e7 — for AUC-scale quantities. */
export function sci(v, dp = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  if (v === 0) return (0).toFixed(dp);
  return v.toExponential(dp).replace('e+', 'e');
}

/**
 * Pick a readable representation automatically: scientific outside
 * [1e-3, 1e6), fixed inside. Keeps a metrics column from mixing widths badly.
 */
export function auto(v, dp = 3) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return sci(v, 2);
  if (a >= 1000) return num(v, 0);
  if (a >= 100) return num(v, 1);
  if (a >= 1) return num(v, dp - 1);
  return num(v, dp);
}

/** Fraction in [0,1] -> percent string. */
export function pct(v, dp = 1) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  return (v * 100).toFixed(dp);
}

/** Seconds -> compact h/m/s. */
export function dur(s) {
  if (s === null || s === undefined || Number.isNaN(s)) return NULL_GLYPH;
  s = Math.max(0, s);
  if (s < 60) return s.toFixed(s < 10 ? 1 : 0) + 's';
  const m = Math.floor(s / 60), r = Math.round(s % 60);
  if (m < 60) return `${m}m ${String(r).padStart(2, '0')}s`;
  const hh = Math.floor(m / 60);
  return `${hh}h ${String(m % 60).padStart(2, '0')}m`;
}

/** Bytes -> KiB/MiB. */
export function bytes(b) {
  if (b === null || b === undefined) return NULL_GLYPH;
  if (b < 1024) return b + ' B';
  if (b < 1048576) return (b / 1024).toFixed(1) + ' KiB';
  return (b / 1048576).toFixed(2) + ' MiB';
}

/** ISO-ish local timestamp for the run table. */
export function ts(x) {
  if (!x) return NULL_GLYPH;
  const d = typeof x === 'number' ? new Date(x * (x > 1e11 ? 1 : 1000)) : new Date(x);
  if (Number.isNaN(d.getTime())) return String(x);
  const p = (n) => String(n).padStart(2, '0');
  return `${d.getFullYear()}-${p(d.getMonth() + 1)}-${p(d.getDate())} ${p(d.getHours())}:${p(d.getMinutes())}:${p(d.getSeconds())}`;
}

/** Signed delta, used by the comparison view. */
export function delta(a, b) {
  if (a === null || a === undefined || b === null || b === undefined) return null;
  if (Number.isNaN(a) || Number.isNaN(b)) return null;
  return b - a;
}

export function signed(v, dp = 3) {
  if (v === null || v === undefined || Number.isNaN(v)) return NULL_GLYPH;
  const s = auto(Math.abs(v), dp);
  return (v > 0 ? '+' : v < 0 ? '-' : '') + s;
}

/** Relative change as a percentage of the reference value. */
export function relPct(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b) || a === 0) return null;
  return (b - a) / Math.abs(a);
}

/** 1-2-5 nice number at or below x — used by the scale bar and axis ticks. */
export function nice125(x) {
  if (!(x > 0)) return 1;
  const e = Math.pow(10, Math.floor(Math.log10(x)));
  const f = x / e;
  return (f >= 5 ? 5 : f >= 2 ? 2 : 1) * e;
}

/** Micrometre length: switches to mm above 1000 um, keeps the unit explicit. */
export function um(v, dp = 1) {
  if (v === null || v === undefined || Number.isNaN(v)) return { v: NULL_GLYPH, u: 'um' };
  if (Math.abs(v) >= 1000) return { v: num(v / 1000, 2), u: 'mm' };
  return { v: num(v, dp), u: 'um' };
}

export function clamp(v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; }

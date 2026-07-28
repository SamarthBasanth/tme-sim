// ramps.js — perceptual colour ramps for scalar field overlays.
//
// These are the real matplotlib control points (17 stops each), interpolated
// to a 256-entry LUT and uploaded as a 256x1 RGBA8 texture. They are
// perceptually uniform (viridis / inferno / magma / cividis) so that equal
// steps in the scalar map to equal perceived steps — the whole point of using
// them instead of a hue rainbow. `cividis` is the CVD-safe default offer;
// `grey` is there because sometimes the honest thing is no colour at all.
//
// A ramp is never shown without a legend carrying its units and its min/max.

const STOPS = {
  viridis: [[68,1,84],[72,24,106],[71,45,123],[66,64,134],[59,82,139],[51,99,141],[44,114,142],[38,130,142],[33,145,140],[31,160,136],[40,174,128],[63,188,115],[94,201,98],[132,212,75],[173,220,48],[216,226,25],[253,231,37]],
  inferno: [[0,0,4],[11,7,36],[33,12,74],[61,9,101],[87,16,110],[113,25,110],[138,34,106],[163,44,97],[188,55,84],[210,70,68],[228,90,49],[241,115,29],[249,142,9],[252,172,17],[249,203,53],[242,234,105],[252,255,164]],
  magma: [[0,0,4],[10,8,34],[29,17,71],[54,16,107],[81,18,124],[106,28,129],[131,38,129],[156,46,127],[183,55,121],[208,65,111],[231,82,99],[245,107,92],[252,137,97],[254,167,114],[254,196,136],[253,226,163],[252,253,191]],
  cividis: [[0,34,78],[0,46,106],[26,56,111],[50,67,109],[67,78,108],[83,90,109],[97,101,111],[111,112,115],[125,124,120],[140,136,120],[155,148,118],[171,160,114],[188,174,108],[205,187,99],[222,201,88],[240,216,70],[254,232,56]],
  coolwarm: [[59,76,192],[78,104,216],[98,130,234],[119,154,247],[141,176,254],[163,194,254],[185,208,249],[204,217,237],[221,220,220],[236,211,197],[245,196,172],[247,176,147],[244,152,122],[235,125,98],[221,95,75],[202,59,55],[180,4,38]],
  grey: [[8,8,8],[24,24,24],[40,40,40],[56,56,56],[72,72,72],[88,88,88],[104,104,104],[120,120,120],[136,136,136],[152,152,152],[168,168,168],[184,184,184],[200,200,200],[216,216,216],[232,232,232],[244,244,244],[255,255,255]],
};

export const RAMP_NAMES = Object.keys(STOPS);

/** Human labels + a one-line "why you would pick this". */
export const RAMP_META = {
  viridis:  { label: 'viridis',  note: 'perceptually uniform, CVD-safe' },
  inferno:  { label: 'inferno',  note: 'uniform, high dynamic range' },
  magma:    { label: 'magma',    note: 'uniform, dark-biased' },
  cividis:  { label: 'cividis',  note: 'optimised for deuteranopia' },
  coolwarm: { label: 'coolwarm', note: 'diverging — use with a centred scale' },
  grey:     { label: 'grey',     note: 'no hue; luminance only' },
};

const cache = new Map();

/** 256*4 Uint8Array RGBA LUT for a named ramp. */
export function lut(name) {
  if (cache.has(name)) return cache.get(name);
  const stops = STOPS[name] || STOPS.viridis;
  const out = new Uint8Array(256 * 4);
  const n = stops.length - 1;
  for (let i = 0; i < 256; i++) {
    const t = (i / 255) * n;
    const i0 = Math.min(n, Math.floor(t));
    const i1 = Math.min(n, i0 + 1);
    const f = t - i0;
    const a = stops[i0], b = stops[i1];
    out[i * 4 + 0] = Math.round(a[0] + (b[0] - a[0]) * f);
    out[i * 4 + 1] = Math.round(a[1] + (b[1] - a[1]) * f);
    out[i * 4 + 2] = Math.round(a[2] + (b[2] - a[2]) * f);
    out[i * 4 + 3] = 255;
  }
  cache.set(name, out);
  return out;
}

/** css rgb() for a single sample in [0,1] — used by legends and read-outs. */
export function sample(name, t) {
  const l = lut(name);
  const i = Math.max(0, Math.min(255, Math.round(t * 255))) * 4;
  return `rgb(${l[i]},${l[i + 1]},${l[i + 2]})`;
}

/** A `linear-gradient(...)` string for painting a legend swatch in CSS. */
export function gradientCss(name, steps = 24) {
  const parts = [];
  for (let i = 0; i < steps; i++) {
    const t = i / (steps - 1);
    parts.push(`${sample(name, t)} ${(t * 100).toFixed(1)}%`);
  }
  return `linear-gradient(90deg, ${parts.join(',')})`;
}

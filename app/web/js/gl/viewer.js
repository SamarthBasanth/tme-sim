// viewer.js — WebGL2 cell viewer.
//
// Budget: 850,000 cells at >= 30 fps with continuous zoom/pan and per-cell
// hover picking. The design decisions that buy that:
//
//  1. ONE draw call for the whole population. gl.POINTS, a single static VBO,
//     no instancing (instancing a quad would quadruple the vertex work for no
//     gain — a point sprite already gives us a screen-space-sized quad free).
//  2. Positions uploaded as normalised uint16 over the section extent rather
//     than f32: 3.4 MB instead of 6.8 MB of vertex fetch, quantisation 0.12 um
//     at an 8.1 mm extent (~1/80 of a cell), i.e. sub-visible.
//  3. Channel culling happens in the vertex shader by emitting a zero-area
//     point off-frustum, so hidden cell types cost no rasterisation.
//  4. Field colouring samples an R32F texture in the VERTEX shader, so there
//     is no per-cell colour buffer to re-upload when the overlay changes.
//  5. Picking is a CPU bin grid (picker.js), never readPixels.
//  6. Rendering is dirty-flag driven; the loop idles at 0% GPU when nothing
//     moves. An adaptive LOD (draw a prefix of a spatially shuffled buffer)
//     engages only while the camera is actually in motion, and only if the
//     measured frame time exceeds the budget.

import { POINTS_VS, POINTS_FS, FIELD_VS, FIELD_FS, LINE_VS, LINE_FS } from './shaders.js';
import { GridIndex } from './picker.js';
import { lut } from './ramps.js';

import { MECH_VARS } from '../mechanism.js';

export const CELL_TYPES = ['TUMOR', 'STROMAL', 'IMMUNE', 'ENDOTHELIAL', 'NECROTIC'];
export const CELL_TYPE_VARS = ['--ch-tumor', '--ch-stromal', '--ch-immune', '--ch-endothelial', '--ch-necrotic'];

/** Colour modes accepted by setColorMode(). */
export const COLOR_MODE = { TYPE: 0, FIELD: 1, LABEL: 2 };

/** Size of the u_labelColor / u_labelVis uniform arrays in shaders.js. */
const MECH_SLOTS = 8;

// Constant "everything visible" uniform payloads for the hover/selection
// marker pass (see drawMark() in render()): a pinned or hovered cell is drawn
// with every channel and every mechanism category forced on, and that mask
// never varies, so it is allocated once here rather than with `new
// Float32Array(...)` on every draw call of every frame a marker is showing.
const ALL_TYPES_VISIBLE = new Float32Array([1, 1, 1, 1, 1]);
const ALL_LABELS_VISIBLE = new Float32Array(MECH_SLOTS).fill(1);

const TARGET_MS = 1000 / 55;   // engage LOD above this while interacting
const RELAX_MS = 1000 / 90;    // relax LOD below this

/** Drawing-buffer area budget in device pixels. See Viewer.resize(). */
export const MAX_DRAW_PX = 4.2e6;

/**
 * Floor on the drawn cell diameter, in CSS pixels. See Viewer.pointSizePx.
 *
 * A cell is ~9 um. Fitting an 8.1 x 5.1 mm Xenium section to an 840 px panel
 * is ~10 um per CSS pixel, so a to-scale cell is 0.9 CSS px — under one pixel,
 * and at dpr 2 that is a 1.8 device-px speck with no round-sprite (the circle
 * discard only engages above 2.5 px). 40,887 of those cover 4.7% of the
 * drawing buffer: individually visible, collectively not a tissue.
 *
 * Below this floor the dots are drawn LARGER THAN LIFE. That is unavoidable —
 * you cannot draw a sub-pixel object — but it is a real distortion of apparent
 * density, so it is disclosed rather than hidden: `pointOversize` reports the
 * factor and the scale bar prints it whenever it exceeds 1.1x.
 *
 * 1.75 CSS px is the largest floor that still leaves the cloud sparse at fit:
 * pi * (1.75 * dpr / 2)^2 * 40,887 is ~18% of an 840x657 dpr-2 buffer, so
 * overlap stays the exception and the density gradient, the tears and the
 * section outline all survive as differences in local coverage rather than
 * flattening into a solid field.
 */
export const MIN_POINT_CSS_PX = 1.75;

/* ------------------------------------------------------------------ util */

let probe = null;
/** Resolve a CSS custom property to [r,g,b] in 0..1, theme-aware. */
export function cssRGB(varName) {
  if (!probe) {
    probe = document.createElement('span');
    probe.style.cssText = 'position:absolute;left:-9999px;top:0;width:0;height:0';
    document.body.appendChild(probe);
  }
  probe.style.color = `var(${varName})`;
  const c = getComputedStyle(probe).color;
  const m = /rgba?\(([^)]+)\)/.exec(c);
  if (!m) return [1, 0, 1];
  const p = m[1].split(',').map(Number);
  return [p[0] / 255, p[1] / 255, p[2] / 255];
}

function compile(gl, type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
    const log = gl.getShaderInfoLog(s);
    gl.deleteShader(s);
    throw new Error('shader compile failed: ' + log);
  }
  return s;
}

function program(gl, vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vs));
  gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
    throw new Error('program link failed: ' + gl.getProgramInfoLog(p));
  }
  const u = {};
  const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < n; i++) {
    const info = gl.getActiveUniform(p, i);
    const name = info.name.replace(/\[0\]$/, '');
    u[name] = gl.getUniformLocation(p, info.name);
  }
  return { p, u };
}

/* ==================================================================== */

export class Viewer {
  constructor(canvas, opts = {}) {
    this.canvas = canvas;
    this.gl = canvas.getContext('webgl2', {
      alpha: false,
      antialias: false,          // we draw points; MSAA buys nothing and costs fill
      depth: false,
      stencil: false,
      desynchronized: true,
      powerPreference: 'high-performance',
      preserveDrawingBuffer: false,
    });
    if (!this.gl) throw new Error('WebGL2 unavailable');

    const gl = this.gl;
    this.prog = {
      points: program(gl, POINTS_VS, POINTS_FS),
      field: program(gl, FIELD_VS, FIELD_FS),
      line: program(gl, LINE_VS, LINE_FS),
    };

    this.floatLinear = !!gl.getExtension('OES_texture_float_linear');

    // static unit quad, shared by the field pass and the outline pass
    this.quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quad);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
    this.quadLoop = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quadLoop);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 1, 1, 0, 1]), gl.STATIC_DRAW);

    this.quadVao = gl.createVertexArray();
    gl.bindVertexArray(this.quadVao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quad);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    this.loopVao = gl.createVertexArray();
    gl.bindVertexArray(this.loopVao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quadLoop);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    gl.bindVertexArray(null);

    // ramp LUT texture (256x1 RGBA8)
    this.rampTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    this.rampName = null;
    this.setRamp(opts.ramp || 'viridis');

    this.fieldTex = null;
    this.field = null;             // {name, rows, cols, dx, min, max, values, unit}

    // geometry
    this.geom = null;              // {n, ex, ey, xs, ys, types}
    this.vao = null;
    this.posBuf = null;
    this.typeBuf = null;
    this.index = null;             // GridIndex

    // view state
    this.dpr = Math.min(window.devicePixelRatio || 1, 2);
    this.cam = { cx: 0, cy: 0, zoom: 1 };   // zoom = css px per um
    this.cellDiamUm = 9;
    this.sizeBoost = 1;
    this.roundPoints = true;
    this.pointAlpha = 1;
    this.colorMode = 0;            // 0 type, 1 field
    this.fieldVisible = true;
    this.fieldOpacity = 0.55;
    this.fieldGamma = 1;
    // Alpha of the field underlay tracks the field VALUE rather than being
    // flat, so a zero concentration paints nothing instead of painting the
    // ramp's darkest colour across the whole extent. See FIELD_FS in
    // shaders.js for the full argument and the one case where it can mislead.
    this.fieldAlphaRamp = true;
    // Window / level, in FIELD UNITS. null = "not set yet"; setField installs
    // the full data range. The ramp is mapped across [c - w/2, c + w/2] and
    // clamped outside, which is what lets an operator pull structure out of a
    // compressed tail instead of looking at a flat wash.
    this.win = null;
    // Float32Array, not a plain array: this is uploaded to u_typeVis on every
    // dirty frame, and while the camera is moving that is every presented
    // frame. A plain array would need `new Float32Array(this.visible)` at
    // upload time — a fresh 5-element allocation per frame for the whole
    // length of a pan/zoom gesture. Mutating it in place in setVisible() costs
    // nothing extra and removes the allocation entirely.
    this.visible = new Float32Array([1, 1, 1, 1, 1]);
    // Per-cell colour-by-attribute channel (colorMode 2). Currently the
    // mechanism-attribution label. `labels` holds the payload in the ORIGINAL
    // (geometry) cell order for read-out; the GPU buffer holds it in the
    // viewer's internal shuffled order.
    this.labels = null;
    this.labelBuf = null;
    this.labelVis = new Float32Array(MECH_SLOTS).fill(1);   // same reasoning as `visible`
    this.labelCounts = null;
    this.hovered = -1;
    this.selected = -1;
    this.lod = 1;
    this.lodAuto = true;
    this.interacting = false;
    this._interactTimer = 0;

    this.stats = {
      fps: 0, frameMs: 0, drawn: 0, total: 0, lod: 1,
      gpuBytes: 0, indexBytes: 0, mode: 'idle',
    };
    this._frames = [];
    this._lastT = 0;
    this._ema = 0;
    this.dirty = true;
    this.continuous = false;
    this._running = false;
    this.onFrame = null;
    this._syncPx = new Uint8Array(4);

    // Cached CSS box of the canvas, kept current by a ResizeObserver instead
    // of a live getBoundingClientRect() read at every use. cssSize is read on
    // the hot path — once per dirty frame via resize(), and again per raw
    // pointermove via screenToWorld() for the cursor read-out — and every
    // stage overlay (scale bar text, hover read-out) rewrites its DOM text on
    // that same frame cadence while the operator is panning. A direct
    // getBoundingClientRect() there would force a synchronous layout on
    // almost every call: read, then write, then read again next frame, with
    // no frame in the gesture where the browser gets to batch the two. A
    // ResizeObserver on the canvas itself keeps `_cssW`/`_cssH` correct on
    // every real size change — window resize, rail/dock splitter drag,
    // initial mount — without a read anywhere in the render or pointer path.
    this._cssW = 1; this._cssH = 1;
    this._measureCss();
    this._ro = null;
    if (typeof ResizeObserver === 'function') {
      this._ro = new ResizeObserver(() => this._measureCss());
      this._ro.observe(this.canvas);
    }

    this.refreshTheme();
    this._loop = this._loop.bind(this);
  }

  _measureCss() {
    const r = this.canvas.getBoundingClientRect();
    this._cssW = Math.max(1, r.width);
    this._cssH = Math.max(1, r.height);
    this.dirty = true;
  }

  /* ------------------------------------------------------------- theme -- */
  refreshTheme() {
    this.typeColors = CELL_TYPE_VARS.map(cssRGB);
    this.bg = cssRGB('--bg-void');
    this.outline = cssRGB('--line-2');
    this.markColor = cssRGB('--focus');
    this.flatTypeColors = new Float32Array(15);
    this.typeColors.forEach((c, i) => this.flatTypeColors.set(c, i * 3));
    this.flatLabelColors = new Float32Array(MECH_SLOTS * 3);
    MECH_VARS.forEach((v, i) => this.flatLabelColors.set(cssRGB(v), i * 3));
    this.dirty = true;
  }

  /* ---------------------------------------------------------- geometry -- */

  /**
   * @param {{n:number, extent:[number,number], x:Float32Array, y:Float32Array,
   *          type:Uint8Array}} g
   */
  setGeometry(g) {
    const gl = this.gl;
    const n = g.n;
    const ex = g.extent[0] || 1, ey = g.extent[1] || 1;

    // Spatially shuffle once so that drawing the first K vertices is a uniform
    // random subsample of the whole section — that is what makes prefix-LOD
    // visually honest instead of lopping off a corner of the tissue.
    const order = new Uint32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
    let seed = 0x9e3779b9;
    for (let i = n - 1; i > 0; i--) {
      seed ^= seed << 13; seed >>>= 0;
      seed ^= seed >>> 17;
      seed ^= seed << 5; seed >>>= 0;
      const j = seed % (i + 1);
      const t = order[i]; order[i] = order[j]; order[j] = t;
    }

    const xs = new Float32Array(n), ys = new Float32Array(n), types = new Uint8Array(n);
    const packed = new Uint16Array(n * 2);
    const sx = 65535 / ex, sy = 65535 / ey;
    for (let k = 0; k < n; k++) {
      const i = order[k];
      const x = g.x[i], y = g.y[i];
      xs[k] = x; ys[k] = y; types[k] = g.type[i];
      let px = x * sx; px = px < 0 ? 0 : px > 65535 ? 65535 : px;
      let py = y * sy; py = py < 0 ? 0 : py > 65535 ? 65535 : py;
      packed[k * 2] = px; packed[k * 2 + 1] = py;
    }

    if (this.vao) gl.deleteVertexArray(this.vao);
    if (this.posBuf) gl.deleteBuffer(this.posBuf);
    if (this.typeBuf) gl.deleteBuffer(this.typeBuf);

    this.posBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.posBuf);
    gl.bufferData(gl.ARRAY_BUFFER, packed, gl.STATIC_DRAW);

    this.typeBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.typeBuf);
    gl.bufferData(gl.ARRAY_BUFFER, types, gl.STATIC_DRAW);

    this.vao = gl.createVertexArray();
    gl.bindVertexArray(this.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.posBuf);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.UNSIGNED_SHORT, true, 0, 0);   // normalised
    gl.bindBuffer(gl.ARRAY_BUFFER, this.typeBuf);
    gl.enableVertexAttribArray(1);
    gl.vertexAttribPointer(1, 1, gl.UNSIGNED_BYTE, false, 0, 0);   // 0..4 as float
    gl.bindVertexArray(null);

    const t0 = performance.now();
    this.index = new GridIndex(xs, ys, ex, ey, 2);
    const indexMs = performance.now() - t0;

    // Kept so a per-cell attribute payload arriving later (mechanism labels)
    // can be permuted into the same order as the position buffer without
    // re-shuffling the geometry.
    this.order = order;
    this.clearLabels();

    this.geom = { n, ex, ey, xs, ys, types };
    this.counts = new Int32Array(CELL_TYPES.length);
    for (let i = 0; i < n; i++) if (types[i] < CELL_TYPES.length) this.counts[types[i]]++;

    this.stats.total = n;
    this.stats.gpuBytes = packed.byteLength + types.byteLength;
    this.stats.indexBytes = this.index.bytes;
    this.stats.indexMs = indexMs;
    this.hovered = -1; this.selected = -1;
    this.lod = 1;
    this.dirty = true;
    return { indexMs };
  }

  clearGeometry() {
    this.geom = null; this.index = null; this.counts = null; this.order = null;
    this.clearLabels();
    this.hovered = -1; this.selected = -1;
    this.stats.total = 0; this.stats.drawn = 0;
    this.dirty = true;
  }

  /* --------------------------------------------- colour-by-attribute -- */

  /**
   * Install a per-cell u8 attribute (the mechanism-attribution label) as the
   * third vertex attribute, so colorMode 2 paints from it. The payload arrives
   * in the geometry payload's cell order and is permuted here into the
   * viewer's internal shuffled order — the same permutation the positions got,
   * so prefix-LOD stays a uniform subsample of every category.
   *
   * One extra byte per vertex. At 850k cells that is 850 kB of additional
   * vertex fetch against a 3.4 MB position buffer, and it is only fetched at
   * all while a geometry is bound; the branch that reads it is uniform-guarded,
   * so cell-type and field colouring cost exactly what they cost before.
   *
   * @param {Uint8Array|null} labels one code per cell, geometry order
   */
  setLabels(labels) {
    const gl = this.gl;
    if (!labels || !this.geom || !this.order) { this.clearLabels(); return; }
    const n = this.geom.n;
    if (labels.length !== n) {
      throw new Error(
        `mechanism payload has ${labels.length} labels but the geometry has ` +
        `${n} cells — refusing to render a misaligned attribution`);
    }
    const shuffled = new Uint8Array(n);
    for (let k = 0; k < n; k++) shuffled[k] = labels[this.order[k]];

    if (!this.labelBuf) this.labelBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.labelBuf);
    gl.bufferData(gl.ARRAY_BUFFER, shuffled, gl.STATIC_DRAW);
    gl.bindVertexArray(this.vao);
    gl.enableVertexAttribArray(2);
    gl.vertexAttribPointer(2, 1, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.bindVertexArray(null);

    this.labels = labels;
    this.labelCounts = new Int32Array(MECH_SLOTS);
    for (let i = 0; i < n; i++) {
      if (labels[i] < MECH_SLOTS) this.labelCounts[labels[i]]++;
    }
    this.labelVis.fill(1);
    this.stats.labelBytes = shuffled.byteLength;
    this.dirty = true;
  }

  clearLabels() {
    const gl = this.gl;
    this.labels = null;
    this.labelCounts = null;
    this.labelVis.fill(1);
    this.stats.labelBytes = 0;
    if (this.vao) {
      gl.bindVertexArray(this.vao);
      gl.disableVertexAttribArray(2);
      gl.vertexAttrib1f(2, 0);       // constant attribute; no buffer read
      gl.bindVertexArray(null);
    }
    if (this.colorMode === COLOR_MODE.LABEL) this.colorMode = COLOR_MODE.TYPE;
    this.dirty = true;
  }

  setLabelVisible(code, on) {
    if (code < 0 || code >= MECH_SLOTS) return;
    this.labelVis[code] = on ? 1 : 0;
    this.dirty = true;
  }

  /* ------------------------------------------------------------- field -- */

  /** @param {{name,rows,cols,dx,min,max,values:Float32Array,unit?:string}|null} f */
  setField(f) {
    const gl = this.gl;
    if (!f) { this.field = null; this.win = null; this.dirty = true; return; }
    if (!this.fieldTex) {
      this.fieldTex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
      const filt = this.floatLinear ? gl.LINEAR : gl.NEAREST;
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filt);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filt);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    }
    gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, f.cols, f.rows, 0, gl.RED, gl.FLOAT, f.values);
    this.field = f;
    this.resetWindow();
    this.dirty = true;
  }

  /* --------------------------------------------------------- window/level */

  /** Full data range of the loaded field, guarded against a degenerate span. */
  fullRange() {
    const f = this.field;
    if (!f) return [0, 1];
    let lo = f.min, hi = f.max;
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) return [0, 1];
    if (hi <= lo) hi = lo + 1e-9;
    return [lo, hi];
  }

  /** The window actually applied to the ramp: [lo, hi] in field units. */
  windowRange() {
    if (!this.win) return this.fullRange();
    const half = Math.max(1e-12, this.win.w) / 2;
    return [this.win.c - half, this.win.c + half];
  }

  setWindow(centre, width) {
    if (!Number.isFinite(centre) || !Number.isFinite(width)) return;
    this.win = { c: centre, w: Math.max(1e-12, width) };
    this.dirty = true;
  }

  /**
   * Replace only the PIXELS of the already-installed field, keeping its
   * dimensions, its texture object and (optionally) the operator's window.
   *
   * This is the time-lapse path. Re-calling setField() per frame would
   * reallocate a 1.15 MB R32F texture and reset the window on every scrub
   * step, which both costs a driver allocation per frame and makes the
   * animation flicker as the ramp re-normalises to each frame's own range.
   * texSubImage2D uploads into the existing storage instead.
   *
   * @param {Float32Array} values row-major, must match the installed rows*cols
   * @param {number} min  finite min to report for this frame
   * @param {number} max  finite max to report for this frame
   * @param {boolean} keepWindow leave the window/level alone (default true)
   */
  setFieldValues(values, min, max, keepWindow = true) {
    const gl = this.gl;
    const f = this.field;
    if (!f || !this.fieldTex) {
      throw new Error('setFieldValues: no field is installed; call setField first');
    }
    if (values.length !== f.rows * f.cols) {
      throw new Error(
        `setFieldValues: ${values.length} samples for a ${f.rows}x${f.cols} ` +
        'grid — frame does not match the installed field');
    }
    gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, f.cols, f.rows, gl.RED, gl.FLOAT, values);
    f.values = values;
    if (Number.isFinite(min)) f.min = min;
    if (Number.isFinite(max)) f.max = max;
    if (!keepWindow) this.resetWindow();
    this.dirty = true;
  }

  /** Restore the window to the full served range (TMSF min..max). */
  resetWindow() {
    const [lo, hi] = this.fullRange();
    this.win = { c: (lo + hi) / 2, w: hi - lo };
    this.dirty = true;
  }

  setRamp(name) {
    if (name === this.rampName) return;
    const gl = this.gl;
    gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, lut(name));
    this.rampName = name;
    this.dirty = true;
  }

  /** Sample the loaded field at a world position (bilinear, um-space). */
  sampleField(x, y) {
    const f = this.field;
    if (!f || !this.geom) return null;
    const u = (x / this.geom.ex) * (f.cols - 1);
    const v = (y / this.geom.ey) * (f.rows - 1);
    if (!(u >= 0 && v >= 0 && u <= f.cols - 1 && v <= f.rows - 1)) return null;
    const i0 = Math.floor(u), j0 = Math.floor(v);
    const i1 = Math.min(f.cols - 1, i0 + 1), j1 = Math.min(f.rows - 1, j0 + 1);
    const fu = u - i0, fv = v - j0;
    const a = f.values[j0 * f.cols + i0], b = f.values[j0 * f.cols + i1];
    const c = f.values[j1 * f.cols + i0], d = f.values[j1 * f.cols + i1];
    return (a * (1 - fu) + b * fu) * (1 - fv) + (c * (1 - fu) + d * fu) * fv;
  }

  /* -------------------------------------------------------------- view -- */

  get cssSize() {
    // No ResizeObserver support: degrade to the direct read every call had
    // before, rather than risk a stale cached size forever.
    if (!this._ro) {
      const r = this.canvas.getBoundingClientRect();
      return { w: Math.max(1, r.width), h: Math.max(1, r.height) };
    }
    return { w: this._cssW, h: this._cssH };
  }

  fit(pad = 0.04) {
    if (!this.geom) return;
    const { w, h } = this.cssSize;
    const z = Math.min(w / this.geom.ex, h / this.geom.ey) * (1 - pad);
    this.cam.zoom = z > 0 ? z : 1;
    this.cam.cx = this.geom.ex / 2;
    this.cam.cy = this.geom.ey / 2;
    this.dirty = true;
  }

  screenToWorld(px, py) {
    const { w, h } = this.cssSize;
    return {
      x: this.cam.cx + (px - w / 2) / this.cam.zoom,
      y: this.cam.cy + (py - h / 2) / this.cam.zoom,
    };
  }

  worldToScreen(x, y) {
    const { w, h } = this.cssSize;
    return {
      x: (x - this.cam.cx) * this.cam.zoom + w / 2,
      y: (y - this.cam.cy) * this.cam.zoom + h / 2,
    };
  }

  zoomAt(px, py, factor) {
    const before = this.screenToWorld(px, py);
    const minZ = this.geom ? Math.min(this.cssSize.w / this.geom.ex, this.cssSize.h / this.geom.ey) * 0.15 : 1e-4;
    this.cam.zoom = Math.max(minZ, Math.min(64, this.cam.zoom * factor));
    const after = this.screenToWorld(px, py);
    this.cam.cx += before.x - after.x;
    this.cam.cy += before.y - after.y;
    this._touch();
  }

  panByScreen(dx, dy) {
    this.cam.cx -= dx / this.cam.zoom;
    this.cam.cy -= dy / this.cam.zoom;
    this._touch();
  }

  _touch() {
    this.dirty = true;
    this.interacting = true;
    clearTimeout(this._interactTimer);
    this._interactTimer = setTimeout(() => {
      this.interacting = false;
      if (this.lodAuto && this.lod < 1) { this.lod = 1; this.dirty = true; }
      this.dirty = true;
    }, 200);
  }

  setVisible(typeIdx, on) { this.visible[typeIdx] = on ? 1 : 0; this.dirty = true; }
  setColorMode(m) { this.colorMode = m; this.dirty = true; }

  /** um per device pixel — used for the scale bar and pick radius. */
  get umPerPx() { return 1 / this.cam.zoom; }

  /** True cell diameter at the current zoom, in CSS px — no floor, no boost. */
  get trueDiamCssPx() { return this.cellDiamUm * this.cam.zoom; }

  /** Drawn cell diameter in CSS px: the legibility floor, then the operator's
   *  boost. The boost is applied AFTER the floor so that it stays fully
   *  authoritative — an operator who drags the size control down still gets
   *  smaller dots instead of silently hitting the clamp. */
  get pointSizeCssPx() {
    return Math.max(MIN_POINT_CSS_PX, this.trueDiamCssPx) * this.sizeBoost;
  }

  /**
   * How much larger than life the dots are currently drawn, 1 = to scale.
   * Covers both causes: the sub-pixel floor, and the operator's size boost.
   * The scale bar surfaces it; see MIN_POINT_CSS_PX for why the floor exists.
   */
  get pointOversize() {
    const d = this.trueDiamCssPx;
    return d > 0 ? this.pointSizeCssPx / d : 1;
  }

  get pointSizePx() {
    // Floor in CSS px, THEN scale to device px. Flooring in device px (as this
    // used to, via Math.max(1, ...)) makes the on-screen dot size depend on the
    // display: 1 device px is a full CSS pixel on a 1x panel and half of one on
    // a retina panel, so the same section looked twice as faint on the better
    // display. The floor is a legibility floor, so it belongs in CSS px.
    return Math.max(1, Math.min(28, this.pointSizeCssPx * this.dpr));
  }

  /* ----------------------------------------------------------- picking -- */

  pick(px, py) {
    if (!this.index) return -1;
    const w = this.screenToWorld(px, py);
    const r = Math.max(this.cellDiamUm * 0.6, 6 / this.cam.zoom);
    return this.index.nearest(w.x, w.y, r, this.geom.types, this.visible);
  }

  cellInfo(i, compositionRadiusUm = 60) {
    if (!this.geom || i < 0 || i >= this.geom.n) return null;
    const x = this.geom.xs[i], y = this.geom.ys[i];
    const comp = this.index.composition(x, y, compositionRadiusUm, this.geom.types, CELL_TYPES.length);
    return {
      index: i, x, y,
      type: this.geom.types[i],
      typeName: CELL_TYPES[this.geom.types[i]] || 'UNKNOWN',
      composition: comp,
      compositionRadiusUm,
      fieldValue: this.sampleField(x, y),
      // `labels` is in the payload's own order; `order[i]` maps the viewer's
      // internal index back to it.
      label: this.labels && this.order ? this.labels[this.order[i]] : null,
    };
  }

  /* ------------------------------------------------------------ render -- */

  start() {
    if (this._running) return;
    this._running = true;
    this._lastT = performance.now();
    requestAnimationFrame(this._loop);
  }

  stop() { this._running = false; }

  _loop(t) {
    if (!this._running) return;
    const dt = t - this._lastT;
    this._lastT = t;

    if (this.dirty || this.continuous) {
      const t0 = performance.now();
      this.render();
      const ms = performance.now() - t0;
      this.stats.cpuMs = ms;
      this.dirty = false;

      if (dt > 0 && dt < 500) {
        this._ema = this._ema ? this._ema * 0.85 + dt * 0.15 : dt;
        this.stats.frameMs = this._ema;
        this.stats.fps = 1000 / this._ema;
        this._frames.push(dt);
        if (this._frames.length > 240) this._frames.shift();
      }

      // Adaptive LOD: only while the camera is moving, only if we are late.
      if (this.lodAuto && this.geom && this.interacting) {
        if (this._ema > TARGET_MS && this.lod > 0.12) {
          this.lod = Math.max(0.12, this.lod * 0.82);
        } else if (this._ema < RELAX_MS && this.lod < 1) {
          this.lod = Math.min(1, this.lod * 1.14);
        }
      }
      this.stats.lod = this.lod;
      this.stats.mode = this.continuous ? 'continuous' : (this.interacting ? 'interacting' : 'on-demand');
      if (this.onFrame) this.onFrame(this.stats);
    }
    requestAnimationFrame(this._loop);
  }

  resize() {
    const { w, h } = this.cssSize;
    // Device-pixel budget. A points renderer gains almost nothing visually
    // from >1.5x oversampling — the sprites are 1-10 px — but fill and
    // compositing cost scales linearly with framebuffer area. On a 4K
    // workstation panel an uncapped dpr-2 buffer is ~8 Mpx and that alone
    // costs more than the 850k draw. Cap total area, then derive dpr.
    const raw = Math.min(window.devicePixelRatio || 1, 2);
    const area = Math.max(1, w * h);
    const dpr = Math.max(1, Math.min(raw, Math.sqrt(MAX_DRAW_PX / area)));
    const W = Math.max(1, Math.round(w * dpr)), H = Math.max(1, Math.round(h * dpr));
    if (this.canvas.width !== W || this.canvas.height !== H || this.dpr !== dpr) {
      this.canvas.width = W; this.canvas.height = H; this.dpr = dpr;
      this.dirty = true;
    }
    return { W, H, w, h };
  }

  render() {
    const gl = this.gl;
    const { W, H, w, h } = this.resize();
    gl.viewport(0, 0, W, H);
    gl.disable(gl.DEPTH_TEST);
    gl.clearColor(this.bg[0], this.bg[1], this.bg[2], 1);
    gl.clear(gl.COLOR_BUFFER_BIT);

    if (!this.geom) { this.stats.drawn = 0; return; }

    const z = this.cam.zoom;
    const scale = [2 * z / w, -2 * z / h];
    const offset = [-this.cam.cx * 2 * z / w, this.cam.cy * 2 * z / h];
    const extent = [this.geom.ex, this.geom.ey];

    const [winLo, winHi] = this.windowRange();

    /* --- field underlay ------------------------------------------------ */
    if (this.field && this.fieldVisible && this.fieldOpacity > 0.001) {
      const { p, u } = this.prog.field;
      gl.useProgram(p);
      gl.bindVertexArray(this.quadVao);
      gl.uniform2fv(u.u_extent, extent);
      gl.uniform2fv(u.u_scale, scale);
      gl.uniform2fv(u.u_offset, offset);
      gl.uniform2f(u.u_range, winLo, winHi);
      gl.uniform1f(u.u_opacity, this.fieldOpacity);
      gl.uniform1f(u.u_gamma, this.fieldGamma);
      gl.uniform1f(u.u_alphaRamp, this.fieldAlphaRamp ? 1 : 0);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
      gl.uniform1i(u.u_field, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
      gl.uniform1i(u.u_ramp, 1);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);   // premultiplied
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    }

    /* --- cells ---------------------------------------------------------- */
    const { p, u } = this.prog.points;
    gl.useProgram(p);
    gl.bindVertexArray(this.vao);
    gl.uniform2fv(u.u_extent, extent);
    gl.uniform2fv(u.u_scale, scale);
    gl.uniform2fv(u.u_offset, offset);
    gl.uniform3fv(u.u_typeColor, this.flatTypeColors);
    gl.uniform1fv(u.u_typeVis, this.visible);
    gl.uniform1i(u.u_colorMode, this.colorMode);
    gl.uniform1i(u.u_hasField, this.field ? 1 : 0);
    gl.uniform1i(u.u_hasLabel, this.labels ? 1 : 0);
    gl.uniform3fv(u.u_labelColor, this.flatLabelColors);
    gl.uniform1fv(u.u_labelVis, this.labelVis);
    gl.uniform2f(u.u_fieldRange, winLo, winHi);
    gl.uniform4f(u.u_force, 0, 0, 0, 0);
    gl.uniform1f(u.u_alpha, this.pointAlpha);
    gl.uniform1f(u.u_ring, 0);

    const psz = this.pointSizePx;
    gl.uniform1f(u.u_pointSize, psz);
    gl.uniform1f(u.u_round, this.roundPoints && psz > 2.5 ? 1 : 0);
    if (this.field) {
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
      gl.uniform1i(u.u_field, 0);
    }
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
    gl.uniform1i(u.u_ramp, 1);

    if (this.pointAlpha < 1) {
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    } else {
      gl.disable(gl.BLEND);
    }

    const count = Math.max(1, Math.min(this.geom.n, Math.ceil(this.geom.n * this.lod)));
    gl.drawArrays(gl.POINTS, 0, count);
    this.stats.drawn = count;

    /* --- hover / selection markers -------------------------------------- */
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.uniform1f(u.u_alpha, 1);
    const drawMark = (idx, color, size) => {
      if (idx < 0 || idx >= this.geom.n) return;
      gl.uniform1f(u.u_pointSize, size);
      gl.uniform1f(u.u_round, 1);
      gl.uniform1f(u.u_ring, 1);
      gl.uniform4f(u.u_force, color[0], color[1], color[2], 1);
      gl.uniform1fv(u.u_typeVis, ALL_TYPES_VISIBLE);
      // A pinned cell must stay visible even when its channel — or its
      // mechanism category — is toggled off; otherwise the marker silently
      // vanishes and the read-out describes something the operator cannot see.
      gl.uniform1fv(u.u_labelVis, ALL_LABELS_VISIBLE);
      gl.drawArrays(gl.POINTS, idx, 1);
    };
    if (this.selected >= 0) drawMark(this.selected, this.markColor, Math.max(16, psz + 10));
    if (this.hovered >= 0 && this.hovered !== this.selected) {
      drawMark(this.hovered, [1, 1, 1], Math.max(12, psz + 7));
    }
    gl.uniform1f(u.u_ring, 0);
    gl.uniform4f(u.u_force, 0, 0, 0, 0);

    /* --- section outline ------------------------------------------------ */
    {
      const { p: lp, u: lu } = this.prog.line;
      gl.useProgram(lp);
      gl.bindVertexArray(this.loopVao);
      gl.uniform2fv(lu.u_extent, extent);
      gl.uniform2fv(lu.u_scale, scale);
      gl.uniform2fv(lu.u_offset, offset);
      gl.uniform4f(lu.u_color, this.outline[0], this.outline[1], this.outline[2], 0.85);
      gl.drawArrays(gl.LINE_LOOP, 0, 4);
    }
    gl.bindVertexArray(null);
  }

  /* --------------------------------------------------------- benchmark -- */

  /** Force the GPU to complete the frame before returning.
   *  gl.finish() is NOT a fence in WebGL — it returns as soon as the command
   *  buffer has been handed to the GPU process, so timing around it measures
   *  command submission (microseconds) and nothing else. A 1-pixel readPixels
   *  is the portable way to actually block on completion. */
  _syncGpu() {
    const gl = this.gl;
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, this._syncPx);
  }

  /**
   * Two-phase measurement, because "fps" alone is ambiguous on a vsync-locked
   * display:
   *
   *   PRESENTED — a requestAnimationFrame sweep with full detail. This is what
   *   the operator actually sees, and it is capped at the panel refresh rate.
   *   This is the number the >= 30 fps requirement is about.
   *
   *   GPU-BOUND — a tight loop with a readPixels fence per iteration, which
   *   removes vsync and reports the true cost of one frame. This is the
   *   headroom number: it says how much room is left before the requirement
   *   is at risk on slower hardware or a larger viewport.
   */
  async benchmark({ frames = 180, syncFrames = 90, pan = true } = {}) {
    if (!this.geom) throw new Error('no geometry loaded');
    const prevLod = this.lod, prevAuto = this.lodAuto, prevCont = this.continuous;
    const wasRunning = this._running;
    const cam0 = { ...this.cam };
    // Suspend the normal render loop: leaving it running would draw the scene
    // twice per display frame and halve the number we are trying to measure.
    this._running = false;
    this.lod = 1; this.lodAuto = false; this.continuous = false;
    await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));

    const sweep = (i, total) => {
      if (!pan) return;
      const th = (i / total) * Math.PI * 2;
      this.cam.cx = cam0.cx + Math.cos(th) * this.geom.ex * 0.12;
      this.cam.cy = cam0.cy + Math.sin(th) * this.geom.ey * 0.12;
      this.cam.zoom = cam0.zoom * (1 + 0.35 * (0.5 - 0.5 * Math.cos(th)));
    };

    /* --- phase 1: presented frame rate -------------------------------- */
    const presented = [];
    let last = performance.now();
    let i = 0;
    await new Promise((resolve) => {
      const step = () => {
        sweep(i, frames);
        this.dirty = true;
        this.render();
        const now = performance.now();
        const dt = now - last; last = now;
        if (i > 8) presented.push(dt);       // discard warm-up
        i++;
        if (i < frames) requestAnimationFrame(step); else resolve();
      };
      requestAnimationFrame(step);
    });

    /* --- phase 2: GPU-bound cost -------------------------------------- */
    this.render(); this._syncGpu();          // prime
    const gpu = [];
    for (let k = 0; k < syncFrames; k++) {
      sweep(k, syncFrames);
      const t0 = performance.now();
      this.dirty = true;
      this.render();
      this._syncGpu();
      gpu.push(performance.now() - t0);
    }

    this.cam = cam0;
    this.lod = prevLod; this.lodAuto = prevAuto; this.continuous = prevCont;
    this.dirty = true;
    if (wasRunning) { this._running = false; this.start(); }

    const stat = (arr) => {
      const a = arr.slice().sort((x, y) => x - y);
      const q = (f) => a[Math.min(a.length - 1, Math.max(0, Math.round(f * (a.length - 1))))];
      const mean = a.reduce((s, v) => s + v, 0) / a.length;
      return { n: a.length, meanMs: mean, medianMs: q(0.5), p95Ms: q(0.95), minMs: a[0], maxMs: a[a.length - 1] };
    };
    const p = stat(presented), g = stat(gpu);
    return {
      n: this.geom.n,
      frames: p.n,
      // presented (what the requirement is measured against)
      meanMs: p.meanMs, medianMs: p.medianMs, p95Ms: p.p95Ms, minMs: p.minMs, maxMs: p.maxMs,
      fpsMean: 1000 / p.meanMs, fpsMedian: 1000 / p.medianMs, fpsP5: 1000 / p.p95Ms,
      // GPU-bound headroom
      gpuMedianMs: g.medianMs, gpuP95Ms: g.p95Ms, gpuMeanMs: g.meanMs,
      gpuFpsMedian: 1000 / g.medianMs,
      pointSizePx: this.pointSizePx,
      dpr: this.dpr,
      canvas: [this.canvas.width, this.canvas.height],
      drawn: this.stats.drawn,
    };
  }
}

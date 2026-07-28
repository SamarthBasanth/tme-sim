// picker.js — CPU spatial index for per-cell hover picking.
//
// Why not GPU colour picking: a readPixels() after every mousemove stalls the
// pipeline and costs more than the whole 850k draw. A uniform bin grid built
// once at load answers a hover in O(cells-per-bin) with zero GPU involvement,
// which is what keeps hover smooth while the viewer is still rendering.
//
// Build cost for 850k cells is one counting sort: two linear passes, ~10 ms,
// and 2 * 4 bytes/cell of index memory.

export class GridIndex {
  /**
   * @param {Float32Array} xs  cell x in um
   * @param {Float32Array} ys  cell y in um
   * @param {number} extentX
   * @param {number} extentY
   * @param {number} targetPerBin
   */
  constructor(xs, ys, extentX, extentY, targetPerBin = 2) {
    const n = xs.length;
    this.xs = xs; this.ys = ys; this.n = n;
    this.ex = extentX || 1; this.ey = extentY || 1;

    const wanted = Math.max(1, Math.ceil(n / targetPerBin));
    const aspect = this.ex / this.ey;
    let nx = Math.max(1, Math.round(Math.sqrt(wanted * aspect)));
    let ny = Math.max(1, Math.round(wanted / nx));
    nx = Math.min(nx, 4096); ny = Math.min(ny, 4096);
    this.nx = nx; this.ny = ny;
    this.bw = this.ex / nx;
    this.bh = this.ey / ny;

    const nbins = nx * ny;
    const counts = new Int32Array(nbins + 1);
    const bin = new Int32Array(n);

    for (let i = 0; i < n; i++) {
      let bx = (xs[i] / this.bw) | 0;
      let by = (ys[i] / this.bh) | 0;
      if (bx < 0) bx = 0; else if (bx >= nx) bx = nx - 1;
      if (by < 0) by = 0; else if (by >= ny) by = ny - 1;
      const b = by * nx + bx;
      bin[i] = b;
      counts[b + 1]++;
    }
    for (let b = 0; b < nbins; b++) counts[b + 1] += counts[b];

    const items = new Int32Array(n);
    const cursor = counts.slice(0, nbins);
    for (let i = 0; i < n; i++) items[cursor[bin[i]]++] = i;

    this.start = counts;     // length nbins+1
    this.items = items;
    this.bytes = (counts.byteLength + items.byteLength);
  }

  /**
   * Nearest cell to (x,y) within `radius` um, restricted to types whose flag
   * in `visible` is truthy. Returns -1 when nothing qualifies.
   */
  nearest(x, y, radius, types, visible) {
    const { nx, ny, bw, bh, start, items, xs, ys } = this;
    const r = Math.max(radius, 1e-6);
    let bx0 = Math.floor((x - r) / bw), bx1 = Math.floor((x + r) / bw);
    let by0 = Math.floor((y - r) / bh), by1 = Math.floor((y + r) / bh);
    if (bx1 < 0 || by1 < 0 || bx0 >= nx || by0 >= ny) return -1;
    if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0;
    if (bx1 >= nx) bx1 = nx - 1; if (by1 >= ny) by1 = ny - 1;

    // Guard against a pathological ring search when the user is zoomed out
    // far enough that the pick radius spans thousands of bins.
    const spanX = bx1 - bx0 + 1, spanY = by1 - by0 + 1;
    if (spanX * spanY > 40000) return -1;

    let best = -1, bestD2 = r * r;
    for (let by = by0; by <= by1; by++) {
      const row = by * nx;
      for (let bx = bx0; bx <= bx1; bx++) {
        const b = row + bx;
        const s = start[b], e = start[b + 1];
        for (let k = s; k < e; k++) {
          const i = items[k];
          if (visible && !visible[types[i]]) continue;
          const dx = xs[i] - x, dy = ys[i] - y;
          const d2 = dx * dx + dy * dy;
          if (d2 < bestD2) { bestD2 = d2; best = i; }
        }
      }
    }
    return best;
  }

  /**
   * Visit every cell inside an axis-aligned um rectangle. Used by the ROI
   * tool. Bin-bounded like `nearest`, so the cost is proportional to the
   * region, not to the section — an ROI over 1% of the tissue reads ~1% of
   * the cells rather than scanning all 850k.
   */
  forEachInRect(x0, y0, x1, y1, cb) {
    const { nx, ny, bw, bh, start, items, xs, ys } = this;
    const ax = Math.min(x0, x1), bx = Math.max(x0, x1);
    const ay = Math.min(y0, y1), by = Math.max(y0, y1);
    let bx0 = Math.floor(ax / bw), bx1 = Math.floor(bx / bw);
    let by0 = Math.floor(ay / bh), by1 = Math.floor(by / bh);
    if (bx1 < 0 || by1 < 0 || bx0 >= nx || by0 >= ny) return 0;
    if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0;
    if (bx1 >= nx) bx1 = nx - 1; if (by1 >= ny) by1 = ny - 1;
    let n = 0;
    for (let byi = by0; byi <= by1; byi++) {
      const row = byi * nx;
      for (let bxi = bx0; bxi <= bx1; bxi++) {
        const s = start[row + bxi], e = start[row + bxi + 1];
        for (let k = s; k < e; k++) {
          const i = items[k];
          const x = xs[i], y = ys[i];
          if (x < ax || x > bx || y < ay || y > by) continue;
          cb(i);
          n++;
        }
      }
    }
    return n;
  }

  /**
   * Local composition: counts per cell type within `radius` um of (x,y),
   * plus the total. This is the "what is around this cell" read-out.
   */
  composition(x, y, radius, types, nTypes = 5) {
    const { nx, ny, bw, bh, start, items, xs, ys } = this;
    const out = new Int32Array(nTypes);
    let total = 0;
    const r2 = radius * radius;
    let bx0 = Math.floor((x - radius) / bw), bx1 = Math.floor((x + radius) / bw);
    let by0 = Math.floor((y - radius) / bh), by1 = Math.floor((y + radius) / bh);
    if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0;
    if (bx1 >= nx) bx1 = nx - 1; if (by1 >= ny) by1 = ny - 1;
    if ((bx1 - bx0 + 1) * (by1 - by0 + 1) > 40000) return { counts: out, total: -1 };

    for (let by = by0; by <= by1; by++) {
      const row = by * nx;
      for (let bx = bx0; bx <= bx1; bx++) {
        const s = start[row + bx], e = start[row + bx + 1];
        for (let k = s; k < e; k++) {
          const i = items[k];
          const dx = xs[i] - x, dy = ys[i] - y;
          if (dx * dx + dy * dy > r2) continue;
          const t = types[i];
          if (t < nTypes) out[t]++;
          total++;
        }
      }
    }
    return { counts: out, total };
  }
}

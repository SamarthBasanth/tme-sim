/**
 * A specimen label for the hero's tissue render.
 *
 * WHY THIS EXISTS. The hero shows 40,887 real measured cells, but with no
 * legend and no dimensions a viewer cannot tell the colours mean anything —
 * it reads as decorative confetti, and the page's claim that the tissue is
 * real lands as an assertion rather than something you can check. A figure
 * caption is how this data is presented everywhere else in the field, and it
 * turns the claim into something legible.
 *
 * Every value here is real and traceable: counts and percentages are the
 * actual cell-type breakdown of the exported section, and the hex values are
 * the same hues `tmesim.viz.explainer.CELLTYPE_COLORS` uses for every figure
 * in `results/`, so the page and the papers agree.
 *
 * DELIBERATELY NOT a graphical scale bar. The camera moves, so a drawn bar
 * would be wrong at most scroll positions, and a wrong scale bar is worse
 * than none. The physical extent is stated as text instead.
 */
const CELL_TYPES = [
  { label: "Tumor", pct: "58.2", hex: "#d65d5c" },
  { label: "Immune", pct: "19.5", hex: "#5c93d6" },
  { label: "Stromal", pct: "11.8", hex: "#5cd65c" },
  { label: "Endothelial", pct: "10.5", hex: "#7668ca" },
] as const;

export default function SpecimenCaption() {
  return (
    <figcaption className="specimen" aria-label="Specimen details for the tissue shown">
      <div className="specimen__id">
        Xenium GBM · 319-gene panel
      </div>

      <dl className="specimen__facts">
        <div>
          <dt>Cells</dt>
          <dd>40,887</dd>
        </div>
        <div>
          <dt>Section</dt>
          <dd>8.09 × 5.07 mm</dd>
        </div>
        <div>
          <dt>Depth</dt>
          <dd>not measured</dd>
        </div>
      </dl>

      <ul className="specimen__legend">
        {CELL_TYPES.map((t) => (
          <li key={t.label}>
            <span
              className="specimen__swatch"
              style={{ backgroundColor: t.hex }}
              aria-hidden="true"
            />
            {t.label}
            <span className="specimen__pct">{t.pct}%</span>
          </li>
        ))}
      </ul>

      <p className="specimen__note">
        10x Genomics, CC BY 4.0 · positions and types as measured
      </p>
    </figcaption>
  );
}

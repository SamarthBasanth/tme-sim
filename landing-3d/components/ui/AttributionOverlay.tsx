"use client";

import { forwardRef } from "react";

interface AttributionOverlayProps {
  initiallyHidden?: boolean;
}

/**
 * The per-cell mechanism attribution for one run, shown alongside section 03.
 *
 * WHY THIS EXISTS. Section 03 makes the product's central argument -- that a
 * compound can fail because it never arrives or because it arrives too weak,
 * and those need different fixes -- and it was making it in pure prose, with
 * the actual decomposition buried mid-sentence as "delivery 0.172, dose 0.825,
 * resistance 0.003". That is the one output most worth *showing*: three
 * numbers that sum to one and answer a question a compartmental model cannot
 * even pose. As running text it read as a claim; as a bar it reads as a
 * reading, which is what it is.
 *
 * It also gives the number ONE home. The triple previously lived only inside a
 * copy string, which is how it survived a physics correction that changed it
 * -- commit c6b6921 rewrote this page's ranking claim for exactly that change
 * and missed the number fourteen lines above. Provenance (run id, and which
 * delivery law produced it) now travels with the value.
 *
 * PROVENANCE, and the reason the caption is not decoration. These are the real
 * shares from run r_00067 (app/runs/r_00067/result.json), computed under the
 * delivery law the engine currently SHIPS -- `PDEConfig.transvascular_exchange`
 * defaults to False. The corrected gradient-driven law gives 0.423 / 0.577 /
 * 0.000 on the same run; which law ships is an open decision. When it is made,
 * this component is the single place the number changes.
 *
 * Deliberately NOT shown: absolute kill fraction. It is not a defensible
 * output and it does not go on this page in any form. Shares of a decomposition
 * are a different quantity -- they are relative by construction.
 */
const BANDS = [
  {
    key: "delivery",
    label: "Never arrived",
    sub: "delivery-limited",
    value: 0.172,
    className: "attribution__seg--delivery",
  },
  {
    key: "dose",
    label: "Arrived too weak",
    sub: "dose-limited",
    value: 0.825,
    className: "attribution__seg--dose",
  },
  {
    key: "resistance",
    label: "Resistant",
    sub: "cell state",
    value: 0.003,
    className: "attribution__seg--resistance",
  },
] as const;

const AttributionOverlay = forwardRef<HTMLDivElement, AttributionOverlayProps>(
  function AttributionOverlay({ initiallyHidden = true }, ref) {
    return (
      <div
        ref={ref}
        className="attribution"
        style={initiallyHidden ? { opacity: 0 } : undefined}
      >
        <div className="attribution__header">
          <span>Why they survived</span>
          <span className="attribution__tag">Run r_00067</span>
        </div>

        {/* One bar, three segments, widths straight from the values -- no
            rounding to "nice" numbers, and the resistance sliver stays a
            sliver because that is the finding. */}
        <div
          className="attribution__bar"
          role="img"
          aria-label="Delivery-limited 17.2 percent, dose-limited 82.5 percent, resistant 0.3 percent"
        >
          {BANDS.map((b) => (
            <span
              key={b.key}
              className={`attribution__seg ${b.className}`}
              style={{ width: `${b.value * 100}%` }}
            />
          ))}
        </div>

        <dl className="attribution__rows">
          {BANDS.map((b) => (
            <div key={b.key} className="attribution__row">
              <dt>
                <span className={`attribution__swatch ${b.className}`} aria-hidden="true" />
                {b.label}
                <span className="attribution__sub">{b.sub}</span>
              </dt>
              <dd>{b.value.toFixed(3)}</dd>
            </div>
          ))}
        </dl>

        <p className="attribution__footnote">
          Shares of surviving tumor cells, one run. Under the shipped delivery
          law; the corrected law gives 0.423 / 0.577 / 0.000.
        </p>
      </div>
    );
  },
);

export default AttributionOverlay;

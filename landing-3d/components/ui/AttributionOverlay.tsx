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
 * shares from run r_00084 (app/runs/r_00084/result.json), computed under the
 * delivery law the engine currently SHIPS -- `PDEConfig.transvascular_exchange`
 * defaults to False.
 *
 * RE-POINTED FROM r_00067 ON 2026-08-22, AND THE REASON IS THE ILLUSTRATIVE-
 * EXAMPLE GUARDRAIL IN lib/copy.ts, NOT A NUMBER GOING STALE. r_00067 is
 * PALBOCICLIB (0.172 / 0.825 / 0.003) and the overlay never said so. It sits
 * directly beside beat 03, which now states a measured result about a
 * DIFFERENT compound on the flagship section, and a reader has no way to tell
 * the bar is not that sentence's illustration -- they will read 82.5% where the
 * sentence says 99.8%, and conclude one of the two is wrong. The example must
 * be at least as well supported as the claim it decorates, and the only way to
 * satisfy that here is for the bar to BE the claim's own run.
 *
 * WHAT THE RE-POINT COST, RECORDED SO IT IS A KNOWN OMISSION. The old footnote
 * carried "the corrected law gives 0.423 / 0.577 / 0.000 on the same run". No
 * corrected-law counterpart has been measured for r_00084, and the r_00067 pair
 * cannot be carried over -- it is a different compound, and it appears in no
 * results/ file, so it fails the citation rule at the head of lib/copy.ts
 * regardless. The disclosure itself did NOT go: beat 06's ledger now carries
 * the delivery law as its first row, with the external anchor (12.8x high /
 * 10.5x low against human microdialysis) that this footnote never had. That is
 * a strictly stronger statement of the same limit, and it is sourced.
 *
 * THE NEAR-EMPTY BANDS ARE THE FINDING AND MUST NOT BE ROUNDED AWAY OR
 * RE-BALANCED FOR LEGIBILITY. A bar that is 99.8% one colour is what "the drug
 * arrives almost everywhere and is almost everywhere too weak" looks like. The
 * resistance 0.000 is not a small number either: with the whole population
 * below IC50 the exposure bands bind first by the documented precedence rule,
 * so RESISTANT is exactly zero cells rather than nearly zero. If a future edit
 * wants a more interesting bar it must find a run whose tissue is more
 * interesting, not a rendering that makes this one look like it.
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
    value: 0.002,
    className: "attribution__seg--delivery",
  },
  {
    key: "dose",
    label: "Arrived too weak",
    sub: "dose-limited",
    value: 0.998,
    className: "attribution__seg--dose",
  },
  {
    key: "resistance",
    label: "Resistant",
    sub: "cell state",
    value: 0.0,
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
          <span className="attribution__tag">r_00084 · temozolomide</span>
        </div>

        {/* One bar, three segments, widths straight from the values -- no
            rounding to "nice" numbers, and the two near-empty bands render at
            their real width (one hairline, one nothing) because that is the
            finding. See the note above before "fixing" the look of this. */}
        <div
          className="attribution__bar"
          role="img"
          aria-label="Delivery-limited 0.2 percent, dose-limited 99.8 percent, resistant 0 percent"
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
          Shares of surviving tumor cells, one compound, one section, one run.
          Under the shipped delivery law — which 06 says can invert it.
        </p>
      </div>
    );
  },
);

export default AttributionOverlay;

"use client";

import { forwardRef } from "react";

interface MetricOverlayProps {
  /**
   * The cinematic experience starts this panel invisible and drives its
   * opacity/transform imperatively every scroll frame (via the forwarded
   * ref), so it must mount hidden or it would flash at full opacity before
   * the first ScrollTrigger update. The reduced-motion fallback has no such
   * driver, so it renders this at full opacity from the start.
   */
  initiallyHidden?: boolean;
}

/**
 * Compact instrument-style readout shown alongside the spatial-attribution
 * section. Its own opacity/transform are set imperatively by
 * ScrollExperience (the same ref-driven pattern used for narrative text)
 * rather than via React state, so it stays in lockstep with the 3D
 * animation without triggering re-renders on every scroll frame.
 *
 * The panel's scroll window (metricPanelOpacity, 3/6 -> 4/6) puts it under
 * section 04, so it shows that section's number and nothing else: relative
 * drug exposure retained by tumor cells in the macrophage-densest vs the
 * sparsest quintile of one real section, and the spread between them.
 *
 * The spread is only meaningful against its nulls, so they are shown as rows
 * rather than prose. A well-mixed compartment returns 1.0 by construction —
 * fair, but it is not what bounds the claim. Two spatially-aware nulls do:
 * relocating the same macrophages while preserving their clustering returns
 * 1.66, and a 400 µm block shuffle returns 1.85, which puts the part
 * attributable to where these macrophages actually ARE at 1.08-1.21x, not
 * 2.0x. Showing 2.0 against 1.0 alone overstates the finding by roughly
 * ten-fold.
 *
 * WHAT CHANGED ON 2026-09-20, AND WHY THE PANEL AND NOT ONLY THE COPY. Three
 * defects, all of them the panel flattering the page rather than the page
 * flattering itself.
 *
 * (1) 2.0x HAD TYPOGRAPHIC PRIMACY. It sat on the summary row, under a rule,
 *     in the emphasis colour, while the corrected 1.08-1.21x appeared only in
 *     a ledger in a DIFFERENT BEAT. results/claim_surface_audit_2026-08-30.md
 *     sec. 2.4 graded that SUPERSEDED-soft: every number on the panel was real
 *     and the footnote was honest, and the emphasis was not. A reader skims the
 *     bold number, so the bold number is now the attributable one.
 *
 *     2.0x IS NOT PRINTED ANYWHERE ON THE PANEL NOW, AND NOTHING WAS HIDDEN BY
 *     THAT. The two bars ARE the 2.0x — 0.43 against 0.86, drawn at their real
 *     widths — so setting the same contrast again as a number was not a second
 *     fact, it was a second emphasis, which is exactly what the row had to
 *     lose. The raw ratio is still stated in beat 04's own copy, where it
 *     belongs, next to what takes it away.
 *
 * (2) ONE NULL WAS SHOWN WHERE TWO ARE PRIMARY, and the shown one was the
 *     flattering half: 1.6581 yields the larger excess (1.208x) and 1.8522 the
 *     smaller (1.081x). With only the first on the panel, the 1.08 end of the
 *     interval could not be derived from anything a visitor could see. Both are
 *     PRIMARY in results/immune_sink_impact.md sec. 5.3; neither is a
 *     sensitivity. They travel as one row at their measured precision rather
 *     than as "~1.7", because the rounding is what hid the arithmetic, and as
 *     one row because the panel's height is spoken for — it is reserved out of
 *     beat 04's own column (see the `padding-bottom` band in globals.css), so
 *     every line here is a line the copy beside it does not get.
 *
 * (3) THE TAG SAID "SAMPLE RUN", WHICH IS THE ONE THING IT IS NOT. No run this
 *     product performs produces these numbers: `tmesim/molecular/immune_sink.py`
 *     is imported by its own validation module and by nothing else, so the term
 *     is off by construction. Calling the panel a sample run asserted a product
 *     capability that does not exist, on the strength of a measurement that does.
 *     The tag now says what the panel is.
 *
 * Do NOT re-add "monotone across all five quintiles" as corroboration. It
 * was there and it was removed on measurement: every control reproduces the
 * monotone fall (label permutation 0.7031->0.5002, uniform-in-tissue
 * 0.6942->0.5780), so it discriminates nothing. Details and the Visium
 * invalidity result: results/immune_sink_impact.md.
 *
 * These are illustrative interface content in the sense that a product
 * screenshot is: real values from one run, labeled "sample run", never to
 * be phrased as a predicted patient outcome. Do not add an absolute kill
 * fraction here; it is not a defensible output.
 */
const MetricOverlay = forwardRef<HTMLDivElement, MetricOverlayProps>(function MetricOverlay(
  { initiallyHidden = true },
  ref,
) {
  return (
    <div ref={ref} className="metric-overlay" style={initiallyHidden ? { opacity: 0 } : undefined}>
      <div className="metric-overlay__header">
        <span>Relative drug exposure</span>
        {/* NOT "Sample run" — see (3) in the docstring. This is a validation
            measurement of a term the product does not run. */}
        <span className="metric-overlay__tag">Term ships off</span>
      </div>
      <dl className="metric-overlay__grid">
        <div className="metric-overlay__row">
          <dt>Macrophage-dense</dt>
          <dd>
            <span className="metric-overlay__bar">
              <span className="metric-overlay__bar-fill" style={{ width: "42.91%" }} />
            </span>
            <span className="metric-overlay__word">~0.43</span>
          </dd>
        </div>
        <div className="metric-overlay__row">
          <dt>Macrophage-sparse</dt>
          <dd>
            <span className="metric-overlay__bar">
              <span className="metric-overlay__bar-fill" style={{ width: "85.92%" }} />
            </span>
            <span className="metric-overlay__word">~0.86</span>
          </dd>
        </div>
        {/* The summary row is the ATTRIBUTABLE excess, not the raw ratio.
            The raw ratio has not been dropped — it moved down to sit with the
            two nulls that bound it, which is the only place it means anything. */}
        <div className="metric-overlay__row metric-overlay__row--total">
          <dt>Attributable</dt>
          <dd>1.08–1.21x</dd>
        </div>
      </dl>
      {/* The reference points sit IN the readout rather than in the footnote,
          because the summary row means nothing without them and a footnote is
          where a number goes to be skipped. 1.0 is the floor a well-mixed
          compartment returns by construction; 1.66 and 1.85 are the two
          spatially-aware nulls that actually bound the claim, and dividing the
          bars' own contrast by them is where the summary row above comes
          from. */}
      <dl className="metric-overlay__refs">
        <div>
          <dt>Well-mixed model</dt>
          <dd>1.0x</dd>
        </div>
        <div>
          <dt>Spatially-aware nulls</dt>
          <dd>1.66–1.85x</dd>
        </div>
      </dl>
      {/* "Not measurable on spot data" was imprecise and is now known to be
          wrong in the literal sense: a spot-resolution field reproduces
          essentially the whole effect (1.98 against the observed 2.00). What
          spot data cannot do is LOCALISE an immune cell — its per-spot immune
          content is a deconvolution estimate — so the number can be computed
          there but not attributed to where the macrophages are, which is the
          entire claim. The distinction is the honest one. */}
      <p className="metric-overlay__footnote">
        One section. Needs measured cell positions; spot platforms estimate
        them.
      </p>
    </div>
  );
});

export default MetricOverlay;

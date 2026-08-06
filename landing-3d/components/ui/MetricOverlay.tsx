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
 * sparsest quintile of one real section, and the spread between them. The
 * spread is the point, because a well-mixed compartmental model returns
 * exactly 1.0 there by construction.
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
        <span className="metric-overlay__tag">Sample run</span>
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
        <div className="metric-overlay__row">
          <dt>Spread</dt>
          <dd>2.0x</dd>
        </div>
      </dl>
      <p className="metric-overlay__footnote">
        One run, one tissue section, monotone across all five quintiles. A well-mixed model returns 1.0x here.
      </p>
    </div>
  );
});

export default MetricOverlay;

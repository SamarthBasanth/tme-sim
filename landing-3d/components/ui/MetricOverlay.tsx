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
 * Compact clinical-style readout for the "AI Analysis" scene. Its own
 * opacity/transform are set imperatively by ScrollExperience (the same
 * ref-driven pattern used for narrative text) rather than via React state,
 * so it stays perfectly in lockstep with the 3D scan-plane animation
 * without triggering re-renders on every scroll frame.
 *
 * The values below are illustrative interface content only — the same way
 * a device manufacturer's own product photography shows a plausible sample
 * reading. They are explicitly labeled "sample case" rather than phrased as
 * a real, validated measurement.
 */
const MetricOverlay = forwardRef<HTMLDivElement, MetricOverlayProps>(function MetricOverlay(
  { initiallyHidden = true },
  ref,
) {
  return (
    <div ref={ref} className="metric-overlay" style={initiallyHidden ? { opacity: 0 } : undefined}>
      <div className="metric-overlay__header">
        <span>Analysis preview</span>
        <span className="metric-overlay__tag">Sample case</span>
      </div>
      <dl className="metric-overlay__grid">
        <div className="metric-overlay__row">
          <dt>Est. volume</dt>
          <dd>14.2 cm³</dd>
        </div>
        <div className="metric-overlay__row">
          <dt>Model confidence</dt>
          <dd>
            <span className="metric-overlay__bar">
              <span className="metric-overlay__bar-fill" style={{ width: "78%" }} />
            </span>
            <span className="metric-overlay__word">High</span>
          </dd>
        </div>
        <div className="metric-overlay__row">
          <dt>Growth trend</dt>
          <dd>Stable</dd>
        </div>
      </dl>
      <p className="metric-overlay__footnote">
        Illustrative output for a sample case — not a diagnosis.
      </p>
    </div>
  );
});

export default MetricOverlay;

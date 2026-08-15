"use client";

import type { SectionLedger as SectionLedgerData } from "@/lib/copy";

/**
 * An in-flow instrument ledger: a named quantity, the verdict the project is
 * registered as holding on it, and the measurement that put it there.
 *
 * WHY IT IS IN FLOW AND NOT A FLOATING PANEL. The page already has two
 * floating readouts (MetricOverlay, AttributionOverlay) and a 3D-anchored
 * callout, and all three have collided with narrative text at some scroll
 * position — the callout at BOTH viewport sizes, because its screen position
 * is a camera projection that sweeps ~130px (desktop) / ~190px (mobile)
 * across its window. A disclosure that can be covered by another element is
 * not a disclosure, so the rows that carry the project's registered limits
 * are laid out in normal flow inside the beat that owns them, where nothing
 * can move on top of them.
 *
 * WHY IT LOOKS LIKE THIS. It borrows SpecimenCaption's register deliberately
 * — mono, hairline rules, no card, no radius, no shadow — because that is the
 * page's most successful existing device and because these rows are the same
 * kind of object: a figure caption in a paper, not a UI component. The state
 * column carries NO colour. Red/amber/green would render a registered limit
 * as an exception to an otherwise working product, and on this page the
 * limits are the product surface, so they are set in the same ink as
 * everything else and separated by weight and tracking alone.
 */
export default function SectionLedger({ caption, rows }: SectionLedgerData) {
  return (
    <div className="ledger">
      <p className="ledger__caption">{caption}</p>
      <dl className="ledger__rows">
        {rows.map((row) => (
          <div key={row.label} className="ledger__row">
            <dt className="ledger__label">{row.label}</dt>
            {/* The verdict sits on the label's baseline, right-set, so the
                column of states can be read straight down without reading
                the prose — the way a spec sheet is actually used. */}
            <dd className="ledger__state">{row.state}</dd>
            <dd className="ledger__detail">{row.detail}</dd>
          </div>
        ))}
      </dl>
    </div>
  );
}

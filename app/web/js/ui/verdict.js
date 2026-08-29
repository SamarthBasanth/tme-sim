// verdict.js — the caveats that qualify the delivery/dose verdict, authored
// ONCE and rendered wherever that verdict appears.
//
// WHY THIS FILE EXISTS
// --------------------
// The verdict is stated twice on the right rail, by design: js/ui/finding.js
// leads with it in counts of cells, and js/ui/reach.js restates it as the
// anchor for the delivery evidence underneath. Both panels therefore have to
// carry the limits that change what it means — a limit that lives a thousand
// words from the number it qualifies does not get read, which is the argument
// finding.js already makes for keeping them inline.
//
// What was wrong was not the repetition but that the two copies were WRITTEN
// TWICE. `joint_limited` was paraphrased independently in each panel, and the
// two paraphrases had already drifted: finding.js carried the cell count
// (19,603) and reach.js did not, so the same caveat was materially weaker on
// one of the two screens that states the claim it qualifies. A caveat that is
// authored in two places is a caveat that will be corrected in one.
//
// So the words live here, once, and both call sites render the same node. That
// is deduplication of the AUTHORING, not of the DISCLOSURE: the caveat still
// appears beside every number it qualifies, which is the rule this repo does
// not bend. Nothing here is ever collapsed, clamped or hidden behind a
// disclosure control — see css/components.css on the caveat list.
//
// The backend also ships its own sentence for this in `joint_limited.note`.
// Neither panel used it. It is now printed VERBATIM underneath the composed
// wording rather than replacing it: the composed sentence names the count and
// the consequence in a bench scientist's words, the backend's is the authority
// the number was produced under, and dropping either loses something. If the
// backend rewords its note the change shows up on screen instead of drifting
// silently away from a paraphrase nobody re-reads.

import { h } from '../util/dom.js';
import { int, pct } from '../util/fmt.js';
import { isSynthetic } from './provenance.js';

/**
 * "Too weak" is not "more would work".
 *
 * @param {object|null} run a RunResult
 * @param {{lead?: string}} [opts] optional replacement for the bold lead, for a
 *   panel whose surrounding wording makes a different opening read better. The
 *   BODY is not parameterised — that is the part that must not vary.
 * @returns {Node|null} null when the payload does not carry the block, in which
 *   case nothing is substituted for it.
 */
export function jointLimitedCaveat(run, opts) {
  const mech = (run && run.mechanism) || null;
  const jl = (mech && mech.joint_limited) || null;
  if (!jl || !Number.isFinite(jl.fraction_of_exposure_limited)) return null;
  const lead = (opts && opts.lead) || '"TOO WEAK" IS NOT "MORE WOULD WORK". ';

  /* The COUNT leads where the payload has one. A share of a share — "92.7% of
   * the exposure-limited ones" — is the shape of number a reader re-anchors to
   * the wrong denominator, and this project has already had that happen to a
   * figure off this very panel. */
  const count = Number.isFinite(jl.count) ? jl.count : null;
  const quantity = count !== null
    ? `${int(count)} of those survivors — `
      + `${pct(jl.fraction_of_exposure_limited, 1)}% of the exposure-limited ones — are`
    : `${pct(jl.fraction_of_exposure_limited, 1)}% of the exposure-limited survivals are`;

  return h('div.caveat', h('span',
    h('b', lead),
    `${quantity} also at or below the modelled resistance threshold. They are `
    + 'counted under their exposure band because below the half-effect '
    + 'concentration the exposure is the binding constraint, but raising the '
    + 'dose alone would not have killed them. Read the verdict as "exposure is '
    + 'the first thing in the way", never as "more drug would fix it".',
    jl.note
      ? h('div', {
        style: {
          fontSize: 'var(--fs-disclosure)', color: 'var(--fg-2)',
          lineHeight: 'var(--lh-body)', margin: 'var(--sp-2) 0 0',
        },
        title: isSynthetic()
          ? 'js/mock.js\'s wording for this limit. There is no backend in this '
            + 'session; against a live one this line is the backend\'s own '
            + 'sentence, verbatim.'
          : 'The backend\'s own wording for this limit, verbatim.',
      }, (isSynthetic() ? 'js/mock.js · ' : 'backend · ') + String(jl.note))
      : null));
}

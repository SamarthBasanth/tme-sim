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

  /* THE LEAD IS SET AT --fs-lg, AND THAT IS THE HOUSE RULE, NOT A PREFERENCE.
   *
   * A caveat is never rendered smaller than the number it qualifies.
   * css/components.css enforces it for `.caveat.crit` and `.caveat.sev-2` by
   * promoting their first `<b>` to --fs-lg, and its own comment names the case
   * it was written for: "the single most load-bearing sentence in the READ
   * rail, sitting directly under a 15px `.m-val`".
   *
   * This caveat is a plain `.caveat`, so it got none of that — and it is the
   * one that qualifies the LARGEST number on either of the two panels that
   * render it: a --fs-xl (19 px) verdict share in js/ui/reach.js and a --fs-lg
   * (15 px) survivor count in js/ui/finding.js, against an 11.5 px lead in
   * both. It stays severity 1, because it qualifies a reading rather than
   * forbidding a quotation; what changes is only that its subject line is no
   * longer small print beside the figure it exists to correct.
   *
   * Set from JS because neither call site owns a stylesheet, and because
   * authoring it here is the entire point of this module — a size fixed in one
   * of the two panels would re-open the drift this file was written to close. */
  const LEAD = {
    fontSize: 'var(--fs-lg)', lineHeight: 'var(--lh-tight)',
    display: 'inline-block', marginBottom: 'var(--sp-2)',
  };

  return h('div.caveat', h('span',
    h('b', { style: LEAD }, lead),
    `${quantity} also at or below the modelled resistance threshold. They are `
    + 'counted under their exposure band because below the half-effect '
    + 'concentration the exposure is the binding constraint, but raising the '
    + 'dose alone would not have killed them. Read the verdict as "exposure is '
    + 'the first thing in the way", never as "more drug would fix it".',
    /* The backend's own sentence, at the SAME size as the caveat it sits
     * inside (`.caveat` is --fs-sm; see the measured argument in
     * components.css). It was set at --fs-disclosure, which made the one line
     * carrying the authority the number was produced under the smallest text
     * in the block — a limit quoted verbatim, set below the paraphrase of it.
     * The mono face and --fg-2 keep it subordinate; subordinate is not the
     * same as smaller. */
    jl.note
      ? h('div', {
        style: {
          fontFamily: 'var(--font-mono)', fontSize: 'var(--fs-sm)',
          color: 'var(--fg-2)', lineHeight: 'var(--lh-body)',
          margin: 'var(--sp-3) 0 0', overflowWrap: 'anywhere',
        },
        title: 'The backend\'s own wording for this limit, verbatim.',
      }, 'backend · ' + String(jl.note))
      : null));
}

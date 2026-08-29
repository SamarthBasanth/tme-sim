// finding.js — THE FIRST SCREEN. What did this run find, in cells, in words a
// bench scientist has never had to learn.
//
// WHY THIS EXISTS
// ---------------
// This panel was written after driving the Studio as someone who had never
// seen it. The product's actual claim — on the flagship section, of the tumour
// cells that survived, essentially all of them had drug ARRIVE and stay too
// weak, and essentially none of them were missed by delivery — was on screen,
// correct, and effectively invisible. It sat five blocks down, behind four
// metrics, under a heading reading "DELIVERY OR DOSE · the discrimination",
// and stated as two percentages with no denominator a reader could hold.
//
// The failure was not that the claim was hidden. It was that reaching it
// required learning a vocabulary first: `E_ref`, `F(0.1)`, "reach ladder",
// "under-reached fraction", "the rung you are using". A first-time reader
// cannot tell an invented term from a term of art they are embarrassed not to
// know, so they slow down, and the finding is behind the slowdown.
//
// So this block leads the rail and does three things, in this order:
//
//   1. states the finding in COUNTS of cells, not fractions. "21,094 of the
//      21,140 tumour cells that survived" is checkable by a reader who has
//      never met this product; "99.8% DOSE-LIMITED" is a term of art wearing a
//      number. The percentage is still shown, beside the count it came from.
//   2. names what would REFUTE it. That is the entire reason this claim
//      replaced a compound ranking — a ranking was only checkable against an
//      outcome nobody has, and this is checkable on a bench. A falsifiable
//      claim that never states its falsification test is being sold as an
//      unfalsifiable one.
//   3. puts the two limits that most change the reading IMMEDIATELY under the
//      sentence they limit, not in the caveats block a thousand words below.
//      Both come from the backend verbatim (`joint_limited.note`,
//      `exposure_axis_note`); neither is composed here.
//
// WHAT THIS BLOCK IS NOT ALLOWED TO DO
// ------------------------------------
// * It does not lead with, headline, or restate `kill_fraction`. Absolute kill
//   is not a defensible output of this project on any section that cannot
//   measure hypoxia, which includes the flagship. The kill count appears here
//   only as the denominator arithmetic — "of N tumour cells, M survived" —
//   because the survival count is meaningless without it, and it is marked as
//   the bounded quantity it is.
// * It does not predict a response, rank a compound, or describe a patient.
// * It invents no number. Every figure is read off `run.mechanism`, which
//   reconciles bit-exactly with the engine's own labels. If the payload does
//   not carry a quantity, this block says so and prints nothing in its place.

import { h, clear } from '../util/dom.js';
import { int, pct, num, NULL_GLYPH } from '../util/fmt.js';
import { doseInputBlock } from './doseinput.js';
import { MECH } from '../mechanism.js';
import { jointLimitedCaveat } from './verdict.js';
import { isSynthetic } from './provenance.js';

/** The limiter buckets, in the words a reader who has not read this repo can
 *  check. Keys are `mechanism.survival_breakdown` keys — see
 *  tmesim/sim/mechanism.py::LABEL_LIMITER, which is the authority. */
const SURVIVAL_PHRASE = {
  delivery: {
    count: 'SURVIVED_UNREACHED',
    lead: 'never had drug arrive in a pharmacologically meaningful amount',
    verdict: 'a DELIVERY result',
    gist: 'On this section, for this compound, the drug did not get to the '
      + 'cells that survived.',
  },
  dose: {
    count: 'SURVIVED_UNDERDOSED',
    lead: 'had drug arrive, and it stayed below the half-effect concentration',
    verdict: 'a POTENCY result, not a delivery result',
    gist: 'On this section, for this compound, the drug arrives almost '
      + 'everywhere and is almost everywhere too weak.',
  },
  resistance: {
    count: 'SURVIVED_RESISTANT',
    lead: 'had adequate exposure, and the modelled cell state suppressed the '
      + 'kill anyway',
    verdict: 'a CELL-STATE result',
    gist: 'On this section, for this compound, exposure was adequate and the '
      + 'modelled state of the cell is what stands in the way.',
  },
  chance: {
    count: 'SURVIVED_STOCHASTIC',
    lead: 'had adequate exposure and no modelled resistance, and drew above '
      + 'their kill probability',
    verdict: 'a result dominated by the model’s own coin flips',
    gist: 'Exposure and modelled biology were both adequate; these cells '
      + 'survived the sampling. Read that as a property of the model, not of '
      + 'the tissue.',
  },
  unattributed: {
    count: 'UNATTRIBUTED',
    lead: 'are not explained by any of the model’s own bands',
    verdict: 'NOT A RESULT — the model cannot account for these cells',
    gist: 'A large unattributed share is itself the finding: the attribution '
      + 'has failed to explain the survivors, and nothing below should be read '
      + 'as an explanation of them.',
  },
};

/** The mechanism label each survival bucket colours the tissue by. Keyed on the
 *  SURVIVAL_PHRASE `count` name so the mapping cannot drift from the labels the
 *  engine actually emits (tmesim/sim/mechanism.py::LABEL_LIMITER). */
const COUNT_TO_MECH = {
  SURVIVED_UNREACHED: MECH.SURVIVED_UNREACHED,
  SURVIVED_UNDERDOSED: MECH.SURVIVED_UNDERDOSED,
  SURVIVED_RESISTANT: MECH.SURVIVED_RESISTANT,
  SURVIVED_STOCHASTIC: MECH.SURVIVED_STOCHASTIC,
  UNATTRIBUTED: MECH.UNATTRIBUTED,
};

/** One count / share / sentence row. */
function row(n, share, text, emphasis, isolate) {
  const node = h('div', {
    style: {
      display: 'grid',
      gridTemplateColumns: 'minmax(0, 5.5em) minmax(0, 4em) minmax(0, 1fr)',
      columnGap: 'var(--sp-4)',
      alignItems: 'baseline',
      padding: 'var(--sp-2) 0',
      borderTop: '1px solid var(--line-1)',
    },
  },
  h('span', {
    style: {
      fontFamily: 'var(--font-mono)', textAlign: 'right',
      fontSize: emphasis ? 'var(--fs-lg)' : 'var(--fs-sm)',
      color: emphasis ? 'var(--fg-0)' : 'var(--fg-1)',
      fontWeight: emphasis ? '600' : '400',
    },
  }, Number.isFinite(n) ? int(n) : NULL_GLYPH),
  h('span', {
    style: {
      fontFamily: 'var(--font-mono)', textAlign: 'right',
      fontSize: 'var(--fs-xs)', color: 'var(--fg-2)',
    },
  }, Number.isFinite(share) ? pct(share, 1) + '%' : NULL_GLYPH),
  h('span', {
    style: {
      fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
      lineHeight: 'var(--lh-body)',
    },
  }, text));

  // THE CLAIM IS A CONTROL. Until 2026-08-27 this panel stated the product's
  // central finding — "21,094 of 21,140 surviving tumour cells had drug arrive
  // and stay too weak" — as text, with no affordance of any kind: zero buttons,
  // zero handlers, zero pointer cursors in the whole file. Meanwhile the viewer
  // could already colour every cell by exactly that label, and isolate one
  // bucket, behind a `colour by: mech` control in the top bar that a first-time
  // reader has no reason to find.
  //
  // So the number and the picture of the number were one click apart, and
  // nothing suggested the click. Clicking a row now colours the tissue by
  // mechanism and shows ONLY that bucket, which turns the sentence into a map
  // of where those cells are — the thing a reader actually wants next and the
  // only way to see that the finding is spatially structured rather than
  // uniform.
  if (isolate) {
    node.tabIndex = 0;
    node.setAttribute('role', 'button');
    node.style.cursor = 'pointer';
    node.title = 'Show only these cells in the tissue view';
    node.addEventListener('click', isolate);
    node.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); isolate(); }
    });
    node.classList.add('finding-row--actionable');
  }
  return node;
}

function para(...kids) {
  return h('p', {
    style: {
      margin: '0 0 var(--sp-4)', fontSize: 'var(--fs-disclosure)',
      lineHeight: 'var(--lh-body)', color: 'var(--fg-1)',
    },
  }, ...kids);
}

export function createFindingPanel({ onIsolate } = {}) {
  const el = h('div');

  function empty(msg) {
    clear(el);
    const p = para(msg);
    p.classList.add('dim');
    el.appendChild(p);
  }

  /** @param {object|null} run a RunResult, or null. */
  function render(run) {
    if (!run) {
      empty('No run is loaded. Submit the run configured on the left — the '
        + 'defaults are a complete, real experiment — and this panel says what '
        + 'it found.');
      return;
    }
    if (run.status !== 'done') {
      empty(`Run ${run.run_id} is ${run.status}. This panel says what a `
        + 'finished run found; nothing is shown for one that has not finished.');
      return;
    }
    const mech = run.mechanism || null;
    const bd = (mech && mech.survival_breakdown) || null;
    const counts = (mech && mech.counts) || null;
    if (!mech || !bd || !counts) {
      // THREE situations land here and a single "no data" would flatten them:
      // a drug-free reconstruction, a run whose backend omitted the
      // attribution, and -- the common one on a hosted instance -- a run no
      // backend produced at all.
      //
      // "That is a backend omission" was what every visitor got whose backend
      // was asleep. There is no backend in stub mode to have omitted anything:
      // the shell already says STUB DATA in the header and js/mock.js is what
      // made the run, so naming that cause sent the reader looking for a
      // server-side failure that never happened.
      const noDrug = !(run.spec && run.spec.drug);
      empty(noDrug
        ? 'This run carried no compound, so there is no delivery-versus-'
          + 'potency finding to state: it reconstructed the tissue and stopped. '
          + 'Pick a drug in the RunSpec on the left and submit again.'
        : isSynthetic()
          ? `Run ${run.run_id} carries no per-cell attribution, and it was not `
            + 'produced by a backend: no server answered, so js/mock.js '
            + 'fabricated this run locally. Nothing here is a result. When the '
            + 'backend is reachable this panel states what the run found.'
          : `Run ${run.run_id} carried no per-cell mechanism attribution, so `
            + 'this panel has nothing to read. That is a backend omission '
            + 'rather than a result, and no number is substituted for it.');
      return;
    }
    clear(el);

    const nTumor = mech.n_tumor;
    const nSurv = mech.n_survived;
    const nKilled = counts.KILLED;

    /* ---- the sentence a reader can check ------------------------------- */
    const drug = (run.spec && run.spec.drug) || 'this compound';
    const section = (run.spec && run.spec.section) || 'this section';
    el.appendChild(h('div', {
      style: {
        fontSize: 'var(--fs-micro)', color: 'var(--fg-2)',
        fontFamily: 'var(--font-mono)', letterSpacing: 'var(--tr-wide)',
        textTransform: 'uppercase', margin: '0 0 var(--sp-3)',
      },
    }, `${drug} · ${section}`));

    el.appendChild(para(
      'Of ', h('b', int(nTumor)), ' tumour cells measured in this section, ',
      h('b', int(nSurv)), ' survived this run. Every cell below is one of those '
      + 'survivors, labelled by what the model says stood in the way.'));

    /* ---- the rows, ordered by size, in counts -------------------------- */
    const entries = Object.entries(bd)
      .filter(([k, v]) => Number.isFinite(v) && v > 0 && SURVIVAL_PHRASE[k])
      .sort((a, b) => b[1] - a[1]);
    entries.forEach(([k, v], i) => {
      const p = SURVIVAL_PHRASE[k];
      const code = COUNT_TO_MECH[p.count];
      const act = (onIsolate && code !== undefined && counts[p.count])
        ? () => onIsolate(code)
        : null;
      el.appendChild(row(counts[p.count], v, p.lead, i === 0, act));
    });
    // The rows above are shares of SURVIVORS and sum to 1 by construction. The
    // backend says so in its own words; quoting it keeps the denominator from
    // being inferred wrongly, which is the single most common way a number off
    // this panel gets misquoted (the shares are NOT of all tumour cells).
    if (mech.survival_breakdown_note) {
      el.appendChild(h('div', {
        style: {
          fontSize: 'var(--fs-micro)', color: 'var(--fg-2)',
          lineHeight: 'var(--lh-body)', borderTop: '1px solid var(--line-1)',
          padding: 'var(--sp-2) 0 0', margin: '0 0 var(--sp-4)',
        },
      }, String(mech.survival_breakdown_note)));
    }

    /* ---- the verdict, in one line -------------------------------------- */
    //
    // THE ONE SENTENCE ON THIS PANEL THAT SURVIVES A SCREENSHOT. Every other
    // marker of a stubbed session -- the striped frame, the magenta viewport
    // tag, the per-value `synthetic` chips -- is chrome around the panel or
    // several paragraphs below this line. This sentence reads as the product's
    // claim wherever it is quoted from, so in stub mode the disclaimer is ON
    // it, in the same block, and the rule is coloured to match.
    const [topKey] = entries[0] || [null];
    if (topKey) {
      const p = SURVIVAL_PHRASE[topKey];
      const stub = isSynthetic();
      el.appendChild(h('div', {
        style: {
          borderLeft: `3px solid var(${stub ? '--synthetic' : '--accent'})`,
          padding: '0 0 0 var(--sp-4)',
          margin: '0 0 var(--sp-4)',
        },
      }, para(
        stub
          ? h('b', { style: { color: 'var(--synthetic)' } },
            'SYNTHETIC — NOT A RESULT. ')
          : null,
        h('b', `This is ${p.verdict}. `), p.gist,
        stub
          ? h('span', { style: { color: 'var(--synthetic)' } },
            ' No backend answered, so js/mock.js fabricated the labels this '
            + 'sentence is counted from. It is a demonstration of the panel, '
            + 'not a finding about any tissue.')
          : null)));
    }

    /* ---- HOW MUCH OF THAT SENTENCE WAS AN INPUT? ------------------------
     * Added 2026-08-22, and it goes DIRECTLY under the verdict because it is
     * about that sentence and nothing else.
     *
     * The shipped default RunSpec carries plasma_peak = 50 uM against a
     * temozolomide IC50 of 606 uM, so the ratio the whole potency half of this
     * panel turns on is 0.083 — fixed in the spec, before any tissue is read —
     * and therapeutic_coverage_1x comes back 0.0. The panel then reported the
     * compound as dose-limited without once saying what dose it was given. A
     * reader's first question is the right one: you told it the peak was 12x
     * under IC50 and it told you the dose was too low; what did the tissue
     * contribute?
     *
     * js/ui/reach.js has answered the structurally identical question for
     * bevacizumab since the panel was written — bbb_penetrance = 0 gives an
     * identically zero source term and it says THIS IS AN INPUT, NOT A
     * FINDING. That branch fires only on bbb == 0. This is the same treatment
     * on the dose axis, and it deliberately picks NO threshold: it prints the
     * ratio. js/ui/doseinput.js carries the full argument for why a cut chosen
     * now would be post-hoc selection. */
    el.append(...doseInputBlock(run));

    /* ---- THE VERDICT'S OWN SENSITIVITY, stated beside it -----------------
     * Added 2026-08-22, and it is the most important caveat on this panel
     * because it is about THIS SENTENCE rather than about a number below it.
     *
     * The delivery/dose discrimination is the product's headline claim. It was
     * measured against the one term the transport law is known to be missing —
     * a transvascular exchange resistance, which ships OFF — and on a shipped
     * compound at shipped parameters IT INVERTS. Not a per-cell nibble: the
     * dominant limiter flips, so the sentence a reader would carry away
     * reverses. A headline that can reverse under a known-missing term, on a
     * page that does not say so, is exactly the failure this product exists to
     * prevent.
     *
     * It is DISCLOSED rather than fixed, deliberately. The split is built and
     * tested (tmesim/treatment/transvascular.py) and is NOT switched on,
     * because making it honest needs a per-compound wall permeability P [um/s]
     * and a tissue surface-to-volume ratio [1/um] that this repository does not
     * have. Inverting the shipped law implies P = 0.0206-0.102 um/s depending
     * on a vessel radius nobody here has measured. Turning it on would trade a
     * law known to be incomplete for one that is unvalidated, which is not
     * progress — it would only move the uncertainty somewhere less visible. */
    const jl = mech.joint_limited;
    const vstab = (run.extra_metrics || {}).verdict_stability || null;
    if (jl) {
      /* THREE STATES, because two rounded one of them off.
       *
       * The census measured TWO reversals per compound and they are NOT the
       * same reversal: this panel's own delivery/dose sentence (5 of 12) and
       * the four-way dominant limiter printed under it (7 of 12). This block
       * used to carry one boolean and the second rate, so IRINOTECAN,
       * PALBOCICLIB and RIBOCICLIB -- limiter flips, verdict holds -- were
       * each told "THIS VERDICT REVERSES FOR <name>", in a sentence that
       * cited by name the file measuring it holding.
       *
       * The tier now arrives from the server as verdict_stability.verdict,
       * where the compound lists are pinned to
       * results/verdict_inversion_census_data.json by
       * tests/test_verdict_stability.py. Absent means NOT SCORED and renders
       * as the panel rate, never as reassurance. */
      const vs = vstab && vstab.verdict ? vstab.verdict
        /* A record written before the tiering carries `stable` but no
         * `verdict`, and if `stable` is false its note may assert a reversal
         * that was never measured. Where it is unambiguous (stable) it is
         * honoured; where it is not, the run falls back to the panel rate and
         * its note is withheld below. A stored claim is not rewritten, but a
         * wrong one is not re-asserted either. */
        : (vstab && vstab.stable === true ? 'holds' : null);
      const who = vstab ? String(vstab.compound).toUpperCase() : '';
      el.appendChild(h('div.caveat.crit', h('span',
        h('b',
          vs === 'reverses'
            ? `THIS VERDICT REVERSES FOR ${who} — THE COMPOUND IN THIS RUN. `
            : vs === 'holds_with_caveat'
              ? `THIS VERDICT HOLDS FOR ${who} — BUT THE FOUR-WAY LIMITER `
                + 'UNDER IT DOES NOT. '
              : vs === 'holds'
                ? `THIS VERDICT HOLDS FOR ${who} — one of the 4 of 12 for `
                  + 'which both measures hold. '
                : 'THIS VERDICT REVERSES FOR 5 OF THE 12 SHIPPED COMPOUNDS. '),
        'The transport law has no transvascular exchange resistance — drug '
        + 'crosses the vessel wall without one, a term the engine ships '
        + 'switched OFF. Switching it on at the shipped parameters, nothing '
        + 'scaled, reverses THIS delivery/dose verdict for 5 of the 12 scored '
        + 'compounds (42%), and flips the four-way dominant limiter printed '
        + 'below it for 7 of 12 (58%). Those are two counts of two different '
        + 'claims, and only the first one is the sentence above. Bevacizumab '
        + 'is excluded, its source is identically zero. Up to 57% of tumour '
        + 'cells change their per-cell label. ',
        vs && vstab && vstab.note
          ? vstab.note + ' Measured on ' + vstab.measured_on + '; see '
            + vstab.source + '.'
          : 'Nothing here says which side this run’s compound falls on: '
            + 'it was not one of the twelve scored, and absence from that '
            + 'census is not evidence of stability.',
        h('br'), h('br'),
        h('b', 'WHICH ONES CAN BE TRUSTED, AND THE RULE IS MECHANICAL. '),
        'Entry elasticity is exactly 1, so the missing wall resistance does '
        + 'one thing — it cuts tissue exposure — and how much it cuts decides '
        + 'what moves. The delivery/dose verdict reverses exactly where '
        + 'exposure collapses to 0.05 to 0.14 of the original: CARMUSTINE, '
        + 'CISPLATIN, LOMUSTINE, PROCARBAZINE, TEMOZOLOMIDE. It survives in '
        + 'the middle band, 0.20 to 0.63, but the four-way limiter flips there: '
        + 'IRINOTECAN, PALBOCICLIB, RIBOCICLIB. Both hold only where most of '
        + 'the exposure survives, 0.69 to 0.93: DOXORUBICIN, ETOPOSIDE, '
        + 'METHOTREXATE, VINCRISTINE. The three bands do not overlap. Read a '
        + 'near-even split as undecided, not as a finding.',
        h('br'), h('br'),
        'Neither arm is the truth. The term is not switched on because doing '
        + 'so needs a per-compound wall permeability and a tissue '
        + 'surface-to-volume ratio that are not measured here — the implied '
        + 'value spans 0.0206 to 0.102 um/s entirely on a vessel radius nobody '
        + 'has measured, which would replace an incomplete law with an '
        + 'unvalidated one. What the census establishes is not that the other '
        + 'answer is right, but that for 5 of the 12 scored compounds this '
        + 'answer is NOT DETERMINED BY THE DATA, and for 3 more the limiter '
        + 'under it is not: two defensible treatments of the same missing '
        + 'physics disagree. There is also a physical bound this law '
        + 'does not respect: tissue exposure rises exactly in proportion to '
        + 'entry without limit, so nothing stops it exceeding the plasma level '
        + 'it came from. See results/verdict_inversion_census.md and '
        + 'results/permeability_split.md.')));
    }

    /* ---- WHAT WOULD REFUTE IT -------------------------------------------
     * The reason this claim replaced a compound ranking. A ranking was only
     * ever checkable against a treatment outcome these sections do not carry;
     * this is checkable with a measurement on the same block of tissue. A
     * falsifiable claim that never names its own falsification test is being
     * presented as an unfalsifiable one, so the test is stated on the same
     * screen as the claim rather than in a document.
     *
     * The test named is a MEASUREMENT, not a prediction: it says what
     * observation would contradict the map, and it does not say what a patient
     * would do. */
    const fold = mech.exposure_shortfall_fold;
    el.appendChild(h('div.caveat.info', h('span',
      h('b', 'WHAT WOULD REFUTE THIS. '),
      `Measure ${drug} directly in this same tissue — imaging mass `
      + 'spectrometry, or microdialysis on a matched block. The claim above '
      + 'says drug is present essentially everywhere in the tumour and is '
      + 'essentially everywhere below the concentration that produces half its '
      + 'effect. Finding tumour regions with no drug in them would refute the '
      + 'delivery half directly. ',
      Number.isFinite(fold)
        ? `Finding concentrations at or above ${num(fold, 2)}x what this run `
          + 'puts in the median tumour cell would refute the potency half. '
        : '',
      'This is the reason this panel exists in place of a compound ranking: a '
      + 'ranking was only ever checkable against a treatment outcome these '
      + 'sections do not carry.')));

    /* ---- the two limits that most change the reading, INLINE ------------
     * Both are the backend's own strings. They are repeated in the Caveats
     * block below with everything else; they are ALSO here because each one
     * changes what the sentence directly above it means, and a limit that
     * lives a thousand words away from the number it limits does not get
     * read. Nothing is removed from Caveats to pay for this. */
    // AUTHORED IN js/ui/verdict.js, not here, and rendered identically by the
    // reach panel. It used to be written out twice — once here with the cell
    // count and once in reach.js WITHOUT it — so the same caveat was materially
    // weaker on one of the two screens that states the claim it qualifies. The
    // caveat still appears beside every number it qualifies; what stopped being
    // duplicated is the wording, not the disclosure.
    const jlCaveat = jointLimitedCaveat(run);
    if (jlCaveat) el.appendChild(jlCaveat);
    if (mech.exposure_axis_note) {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'THESE TWO LABELS ARE ONE AXIS, NOT TWO MECHANISMS. '),
        String(mech.exposure_axis_note))));
    }
    // The band edges are a CHOICE, and the top row's label depends on them.
    // 0.1 is where "never arrived" is drawn; a different edge moves cells
    // between the two rows above without anything physical changing.
    const th = mech.thresholds || null;
    if (th && Number.isFinite(th.unreached_ratio)) {
      el.appendChild(h('div.caveat', h('span',
        h('b', 'THE LINE BETWEEN THE TWO ROWS IS A CONVENTION. '),
        `A survivor is called "never arrived" below ${num(th.unreached_ratio, 2)}x `
        + `the half-effect concentration and "arrived but too weak" between `
        + `there and ${num(th.adequate_ratio, 2)}x. Nothing in the tissue or the `
        + 'compound marks that edge out — it is this project\'s reporting '
        + 'convention, and moving it moves cells between the two rows above '
        + 'with no physics changing. The full ladder of thresholds is plotted '
        + 'under Delivery or dose, below, so the choice stays visible.')));
    }

    /* ---- the arithmetic, and its own bound ------------------------------
     * The kill count is HERE and not at the top, and it is stated as a
     * denominator rather than as a result. Absolute kill is not defensible on
     * a section that cannot measure hypoxia: it moves by up to 0.3252 absolute
     * (1.77x) purely with the choice of state route, with tissue, drug,
     * transport and RNG held identical (results/cell_state_impact.md sec.5).
     * That statement lives in the run's own caveats and in STANDING_LIMITS; it
     * is pointed at from here because this is the one place on this panel
     * where a kill number is printed, and a reader must not be able to read
     * the number without meeting the bound. */
    const rec = mech.reconciliation || null;
    el.appendChild(h('div.caveat.crit', h('span',
      h('b', 'THE KILL COUNT IS THE DENOMINATOR HERE, NOT A RESULT. '),
      `${int(nKilled)} of ${int(nTumor)} tumour cells are labelled killed, which `
      + 'is what leaves the survivor count above. Do not quote that as a kill '
      + 'rate. It is an upper bound with no intrinsic resistance modelled, and '
      + 'on a section that cannot measure hypoxia it moves by up to 0.3252 '
      + 'absolute (1.77x) purely with which state route is chosen, everything '
      + 'else held identical. The shares above are ratios among survivors and '
      + 'do not inherit that movement. See Caveats, below, in full.',
      isSynthetic()
        ? ' These labels were fabricated by js/mock.js. Both sides of the '
          + 'reconciliation below them are the same tally of that one '
          + 'fabricated array, so it is an identity and audits nothing.'
        : rec && rec.bit_exact
        ? ' The labels above reconcile bit-exactly with the engine\'s own '
          + 'kill_fraction, so they are a restatement of its arithmetic rather '
          + 'than a second opinion about it.'
        : '')));
  }

  return { el, render };
}

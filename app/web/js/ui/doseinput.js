// doseinput.js — HOW MUCH OF THE VERDICT WAS AN INPUT?
//
// THE DEFECT THIS FIXES
// ---------------------
// The rail's headline verdict is DELIVERY-LIMITED or DOSE-LIMITED. On the
// shipped default RunSpec its SIGN is handed to the model in the spec, and
// until 2026-08-22 no panel said so:
//
//   plasma_peak       50 uM   (a field in the RunSpec, on the left of this app)
//   temozolomide IC50 606 uM  (GDSC, median of 51 glioma cell lines)
//   ratio             0.083
//
// therapeutic_coverage_1x comes back 0.0, the mechanism attribution comes back
// 99.8% DOSE-LIMITED, and the panel reports the compound as dose-limited. The
// first question a pharmacologist asks is exactly the right one: *you told it
// the peak was 12x under IC50 and it told you the dose was too low — what did
// the tissue contribute?*
//
// js/ui/reach.js already answers the structurally identical question for
// bevacizumab, where bbb_penetrance = 0 gives an identically zero source term
// and the panel says THIS IS AN INPUT, NOT A FINDING. That branch fires only on
// bbb == 0. This module extends the same treatment to the dose axis, and it is
// used by BOTH panels that state the verdict so the wording exists once.
//
// WHY THERE IS NO THRESHOLD IN HERE
// ---------------------------------
// The obvious implementation is "if plasma_peak/IC50 < 0.5, warn". It is also
// the exact mistake this product removed a headline for: a cut chosen NOW is
// chosen knowing what the shipped data does, which is post-hoc selection
// wearing a convention's clothes (see reach.js on the withdrawal of the F(0.1)
// headline). So this module chooses nothing. It PRINTS THE RATIO and lets the
// reader hold it.
//
// The one comparison it does make is against 1.0, and 1.0 is not a chosen
// threshold: it is the definitional point where the plasma peak entered for the
// run equals the concentration that produces half the compound's effect. Below
// it, a statement follows by arithmetic rather than by convention — see
// `plasmaBoundHolds` below for exactly which one, and under which delivery law.
//
// EVERY NUMBER COMES FROM THE PAYLOAD
// -----------------------------------
// Nothing here is a constant and nothing is a fixture:
//
//   plasma_peak_uM, plasma_mean_uM,
//   median_E_over_ic50, p10_E_over_ic50, E_p50_uM   run.extra_metrics
//   ic50_uM                                          run.drug_params
//   therapeutic_coverage_1x                          run.metrics
//   ic50_source / _confidence / _n_cell_lines        run.provenance
//   delivery_law.transvascular_exchange / .law       run.provenance
//
//   plasma_peak_over_ic50, tissue_median_over_plasma_peak,
//   tissue_median_over_plasma_mean                   run.extra_metrics
//
// The three ratios were briefly divided in the client and are now supplied by
// `app/server/jobs.py`, so they carry provenance, can be asserted against and
// land in result.json. Records written before those fields existed fall back to
// the division and SAY they did — see `derived` in `doseFraming()`. If a field
// is absent the block that needs it is not rendered and nothing is substituted
// for it.
//
// WHICH RATIO IS THE BOUND
// ------------------------
// `tissue_median_over_plasma_PEAK`, never the mean. C <= Cp is a statement
// about the highest concentration the blood ever reached; the mean is the
// intuitive comparison and it is contaminated by the dosing schedule, since a
// longer window lowers the mean and inflates the ratio for a reason that has
// nothing to do with the bound. Measured: 3.71 (peak) / 5.96 (mean) on the
// shipped 6 h default — both over 1, that run really does violate it — against
// 0.836 / 0.848 on a short window, BELOW 1 on both readings while the law is
// exactly as unbounded. Reasoning against the mean would report the bound
// satisfied on short windows and violated on long ones, when the only thing
// that changed was the schedule.

import { h } from '../util/dom.js';
import { auto, int, num, pct, NULL_GLYPH } from '../util/fmt.js';

const fmtFrac = (v) => (Number.isFinite(v) ? pct(v, 1) + '%' : NULL_GLYPH);

/* A concentration, at reading precision. `auto` renders an entered 50 as
 * "50.00", which puts three digits of false precision on a number the operator
 * typed; an integer entry prints as an integer and everything else keeps
 * `auto`'s four significant figures. */
const conc = (v) => (Number.isInteger(v) ? int(v) : auto(v, 3));

/**
 * Read the dose-versus-potency framing off a RunResult.
 *
 * @param {object|null} run
 * @returns {null|{
 *   peak_uM: number, ic50_uM: number, ratio: number, foldUnder: number,
 *   mean_uM: number|null, medianOverIc50: number|null, p10OverIc50: number|null,
 *   medianE_uM: number|null, tissueOverPlasmaMean: number|null,
 *   tissueOverPlasmaPeak: number|null, coverage1x: number|null,
 *   ic50Source: string|null, ic50N: number|null, ic50Confidence: string|null,
 *   law: string|null, plasmaBoundHolds: boolean|null,
 * }}
 *
 * Returns null — and therefore renders nothing — unless BOTH halves of the
 * ratio are present and positive. A half-known ratio is not a weaker version of
 * this disclosure, it is a different and unfounded one.
 */
export function doseFraming(run) {
  if (!run) return null;
  const em = run.extra_metrics || {};
  const dp = run.drug_params || {};
  const prov = run.provenance || {};
  const peak = em.plasma_peak_uM;
  const ic50 = dp.ic50_uM;
  if (!Number.isFinite(peak) || !Number.isFinite(ic50) || peak <= 0 || ic50 <= 0) {
    return null;
  }

  const mean = Number.isFinite(em.plasma_mean_uM) ? em.plasma_mean_uM : null;
  const medE = Number.isFinite(em.E_p50_uM) ? em.E_p50_uM : null;

  /* WHICH SOURCE TERM PRODUCED THE MAP, and therefore whether the plasma level
   * bounds it. Read off the run's own provenance, never assumed:
   *
   *   S = PS*(Cp - C)  gradient exchange   C_ss <= Cp BY CONSTRUCTION
   *   S = PS*Cp        one-way leak        no path back to blood, unbounded
   *
   * app/server/delivery_law.py is the authority for both sentences and ships
   * the leak as the default. Absent (an older record) is a THIRD state and is
   * reported as unknown rather than folded into either — that is the same rule
   * reach.js applies to a missing `window_mode`. */
  const dl = prov.delivery_law || null;
  const exch = dl ? dl.transvascular_exchange : undefined;
  const plasmaBoundHolds = typeof exch === 'boolean' ? exch : null;

  /* READ, NOT DIVIDED. All three ratios are supplied by the backend as of
   * 2026-08-22 (app/server/jobs.py, extra_metrics) so that they carry
   * provenance, can be asserted against, and survive into result.json. A
   * number computed in a browser has none of those properties.
   *
   * THE FALLBACK IS A THIRD STATE AND IS LABELLED ONE. There are ~90 run
   * records in this repo written before those fields existed. Refusing to
   * disclose on them would hide the ratio precisely where a reader is most
   * likely to be re-reading an old result; silently dividing would present a
   * browser arithmetic as a recorded quantity. So the division is done, and
   * `derived` says it was — the same rule the reach panel applies to a missing
   * `window_mode`. */
  const supplied = Number.isFinite(em.plasma_peak_over_ic50);
  const ratio = supplied ? em.plasma_peak_over_ic50 : peak / ic50;
  const overPeak = Number.isFinite(em.tissue_median_over_plasma_peak)
    ? em.tissue_median_over_plasma_peak
    : (medE !== null ? medE / peak : null);
  const overMean = Number.isFinite(em.tissue_median_over_plasma_mean)
    ? em.tissue_median_over_plasma_mean
    : (medE !== null && mean !== null && mean > 0 ? medE / mean : null);

  return {
    peak_uM: peak,
    ic50_uM: ic50,
    ratio,
    foldUnder: 1 / ratio,
    derived: !supplied,
    mean_uM: mean,
    medianOverIc50: Number.isFinite(em.median_E_over_ic50) ? em.median_E_over_ic50 : null,
    p10OverIc50: Number.isFinite(em.p10_E_over_ic50) ? em.p10_E_over_ic50 : null,
    medianE_uM: medE,
    tissueOverPlasmaMean: overMean,
    tissueOverPlasmaPeak: overPeak,
    coverage1x: run.metrics && Number.isFinite(run.metrics.therapeutic_coverage_1x)
      ? run.metrics.therapeutic_coverage_1x : null,
    ic50Source: prov.ic50_source ? String(prov.ic50_source) : null,
    ic50N: Number.isFinite(prov.ic50_n_cell_lines) ? prov.ic50_n_cell_lines : null,
    ic50Confidence: prov.ic50_confidence ? String(prov.ic50_confidence) : null,
    law: dl && dl.law ? String(dl.law) : null,
    plasmaBoundHolds,
  };
}

/** The ratio, printed with both of its terms, so the reader can do the division
 *  themselves rather than take it on trust. */
function ratioText(f) {
  return `plasma_peak / IC50 = ${conc(f.peak_uM)} uM / ${conc(f.ic50_uM)} uM `
    + `= ${num(f.ratio, 3)}`;
}

/** Where the IC50 came from, as a NOUN PHRASE that can sit in any sentence.
 *
 *  It is the half of the ratio that is NOT an input, and a reader has to be
 *  able to see that without leaving the panel — otherwise "this is a ratio
 *  between two things you chose" is the obvious and wrong reading.
 *
 *  The `n` is appended only when the source string does not already carry one:
 *  the shipped provenance reads "GDSC measured (n=51)", and appending to it
 *  produced "GDSC measured (n=51), n = 51 cell lines" on screen. */
function ic50Provenance(f) {
  if (!f.ic50Source) return 'a constant this record carries no provenance for';
  const hasN = /\bn\s*=/.test(f.ic50Source);
  const n = !hasN && f.ic50N ? `, n = ${f.ic50N} cell lines` : '';
  return `${f.ic50Source}${n}`;
}

/**
 * THE SAME RATIO, BEFORE THE RUN — for the dose field in the setup form.
 *
 * WHY IT LIVES IN THIS MODULE. Everything above is written for a RunResult, and
 * a RunResult only exists after the operator has already spent the run. The
 * ratio it prints is fixed the moment the dose is typed, which means the panel
 * that explains it is one screen too late: `plasma_peak` is a field of the
 * RunSpec on the LEFT of this app, and until this function existed nothing
 * beside that field said what it was a fraction of. `js/ui/spec.js` calls this
 * under its dose row so the wording exists exactly once — the same rule the
 * header of this file states for the two result panels.
 *
 * WHAT IT MAY NOT DO. It may not predict the verdict. What a sub-IC50 ratio
 * IMPLIES is read off the run's delivery law (`lawClause` above), and no run
 * has happened yet, so the implication is named as pending rather than guessed:
 * under the gradient-exchange source the verdict would be arithmetic, under the
 * shipped one-way leak it would not, and this function knows neither. It prints
 * the two numbers, their quotient, and where the half that is not an input came
 * from. Nothing else.
 *
 * NO THRESHOLD, for the reason at the top of this file. The only comparison is
 * against 1.0, which is definitional rather than chosen.
 *
 * @param {{peak_uM:number, ic50_uM:number, ic50Source?:string|null,
 *          drug?:string|null, userCompound?:boolean}} d
 * @returns {Node|null} null when either half is missing — a half-known ratio is
 *   not a weaker version of this disclosure, it is a different and unfounded
 *   one, which is the same rule `doseFraming()` applies.
 */
export function plasmaPeakNote(d) {
  const peak = d && Number(d.peak_uM);
  const ic50 = d && Number(d.ic50_uM);
  if (!Number.isFinite(peak) || !Number.isFinite(ic50) || peak <= 0 || ic50 <= 0) {
    return null;
  }
  const ratio = peak / ic50;
  const under = ratio < 1;
  const src = d.ic50Source
    ? String(d.ic50Source)
    : 'a constant this catalogue carries no provenance for';
  return h('div', {
    style: { lineHeight: 'var(--lh-body)', color: 'var(--fg-2)' },
  },
  h('span.mono.strong', { style: { color: 'var(--fg-0)' } },
    `plasma_peak / IC50 = ${conc(peak)} uM / ${conc(ic50)} uM = ${num(ratio, 3)}`),
  ' — the peak plasma concentration entered above is ',
  h('b', under ? `${num(1 / ratio, 1)}x under` : `${num(ratio, 2)}x`),
  ' the concentration that produces half this compound\'s effect',
  ` (IC50 ${conc(ic50)} uM, ${src}`,
  d.userCompound
    ? ' — a compound you compiled, so this half of the ratio is your own '
      + 'number or a flagged placeholder'
    : '',
  '), which is not a RunSpec field. The ratio is fixed here, before any tissue '
  + 'is read. What follows from it depends on the delivery law the run uses '
  + 'and is stated on the result, not predicted here.');
}

/** The mono ratio strip — the number, beside the words. Same shape as the
 *  delivery-spread strip in reach.js so the two read as one family. */
function ratioStrip(f, big) {
  return h('div', {
    style: {
      display: 'flex', alignItems: 'baseline', gap: 'var(--sp-4)',
      margin: '0 0 var(--sp-3)',
    },
  },
  h('span', {
    style: {
      fontSize: big ? 'var(--fs-lg)' : 'var(--fs-sm)',
      fontFamily: 'var(--font-mono)', color: 'var(--fg-0)',
    },
  }, num(f.ratio, 3)),
  h('span', {
    style: {
      fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
      lineHeight: 'var(--lh-body)',
    },
  },
  h('b', 'PLASMA PEAK / IC50'), ' — the peak plasma concentration entered for '
    + `this run (${conc(f.peak_uM)} uM, a RunSpec field) over the concentration `
    + 'that produces half this compound\'s effect '
    + `(${conc(f.ic50_uM)} uM, ${ic50Provenance(f)}, which is not a RunSpec `
    + 'field). This ratio is fixed before any tissue is read.',
  /* WHERE THE RATIO CAME FROM. The backend supplies it as
   * `extra_metrics.plasma_peak_over_ic50`, so on a current run this is a
   * RECORDED quantity that is in result.json and can be asserted against. On a
   * record written before that field existed it is a division done in this
   * browser, which has none of those properties — and that difference is
   * stated rather than smoothed over. */
  f.derived
    ? h('span', { style: { color: 'var(--fg-2)' } },
      ' Derived in this browser: this record predates '
      + 'extra_metrics.plasma_peak_over_ic50, so the ratio is not in its '
      + 'result.json and carries no provenance of its own.')
    : null));
}

/**
 * The law clause: does this run's source term bound tissue by plasma?
 *
 * Three branches, because there are three states, and the difference between
 * them is the difference between "the verdict was arithmetic" and "the verdict
 * survived a law that could have contradicted it".
 */
function lawClause(f) {
  const cov = f.coverage1x !== null ? fmtFrac(f.coverage1x) : null;

  if (f.plasmaBoundHolds === true) {
    return [
      'What follows from a ratio below 1 depends entirely on the source term, '
      + 'so it is read off this run\'s own. This run\'s delivery law is the '
      + `gradient exchange ${f.law || 'S = PS*(Cp - C)'}, which is bounded by `
      + 'the plasma level BY CONSTRUCTION (app/server/delivery_law.py). That '
      + 'bound is therefore in force here: no point of this map could reach '
      + 'IC50 whatever the vessels did'
      + (cov ? `, and therapeutic_coverage_1x = ${cov} follows by arithmetic ` : ', ')
      + (cov
        ? 'for every tissue and every vessel geometry — it carries no '
          + 'information about either. '
        : 'at this ratio, for every tissue and every vessel geometry. '),
      h('b', 'The sign of a dose-limited verdict was settled by the RunSpec, '
        + 'not by this tissue'), '.',
    ];
  }

  if (f.plasmaBoundHolds === false) {
    /* THE BOUND IS AGAINST THE PEAK, and the mean is context at most.
     *
     * C <= Cp is a statement about the HIGHEST concentration the blood ever
     * reached: absent an accumulation or binding mechanism, which this engine
     * does not model, tissue cannot exceed it. The plasma MEAN is the intuitive
     * comparison and it is the wrong one — it is contaminated by the schedule.
     * Plasma decays across the window, so a longer window lowers the mean and
     * inflates the ratio for a reason that has nothing to do with the bound.
     *
     * Measured, and this is why it matters rather than being a preference: on
     * the shipped 6 h default the two read ~3.71 (peak) and ~5.96 (mean), both
     * above 1, so that run does violate the bound. On a short window they read
     * 0.836 and 0.848 — BELOW 1 on both readings, while the law is exactly as
     * unbounded as it is on the flagship. A panel that reasoned against the
     * mean would report the bound satisfied on short windows and violated on
     * long ones, when all that changed was the dosing schedule. */
    const amp = [];
    if (f.tissueOverPlasmaPeak !== null && f.medianE_uM !== null) {
      const violates = f.tissueOverPlasmaPeak > 1;
      amp.push(
        ' The bound that law is missing is C <= Cp — tissue cannot exceed the '
        + 'HIGHEST concentration the blood ever reached, absent an accumulation '
        + 'or binding mechanism this engine does not model. On this run the '
        + `median tumour cell's window-averaged concentration is `
        + `${conc(f.medianE_uM)} uM against an entered plasma peak of `
        + `${conc(f.peak_uM)} uM: `,
        h('b', `${num(f.tissueOverPlasmaPeak, 2)}x the peak`),
        violates
          ? '. That is the bound broken, on screen, by this run'
          : '. That sits under the bound on this run, which is NOT evidence the '
            + 'law respects it — it is unbounded either way, and a short window '
            + 'can put this ratio below 1 while nothing about the law has '
            + 'changed',
        cov ? ` — and it still returned therapeutic_coverage_1x = ${cov}.` : '.',
        violates
          ? ' Read the excess as the source term\'s known unboundedness, NOT as '
            + 'this tissue concentrating the compound: nothing in this model '
            + 'does active uptake.'
          : '',
        f.tissueOverPlasmaMean !== null && f.mean_uM !== null
          ? ` For context only, the same concentration is `
            + `${num(f.tissueOverPlasmaMean, 2)}x the plasma MEAN over the `
            + `window (${conc(f.mean_uM)} uM). Do not read that as the bound: `
            + 'the mean falls as the window lengthens, so it inflates with the '
            + 'schedule rather than with the violation.'
          : '',
      );
    }
    return [
      'What follows from a ratio below 1 depends entirely on the source term, '
      + 'so it is read off this run\'s own. This run\'s delivery law is the '
      + `one-way leak ${f.law || 'S = PS*Cp'}, the shipped default, which has `
      + 'no path from interstitium back to blood and whose steady state is '
      + 'therefore UNBOUNDED ABOVE (app/server/delivery_law.py). ',
      h('b', 'So the arithmetic did not force the answer here'),
      ' — this map was free to put cells over IC50 at a ratio below 1.',
      ...amp.flat(),
      ' Take that as WEAK reassurance rather than strong: the law it survived '
      + 'is one this repo has measured against human microdialysis and found '
      + 'wrong in both directions (results/breast_kp_validation.md).',
    ];
  }

  return [
    'This record does not say which delivery law produced it, so what follows '
    + 'from a ratio below 1 CANNOT be read off it. Under the gradient-exchange '
    + 'source, tissue is bounded by plasma and the verdict is arithmetic; under '
    + 'the one-way leak, which is unbounded above, it is not. Neither is being '
    + 'assumed, and no coverage figure here should be read without settling '
    + 'which law produced it.',
  ];
}

/** What survives the input — stated as its own list so it is not read as an
 *  apology. The spread is a within-map ratio and carries no dose; that is the
 *  same argument reach.js makes for the delivery spread, pointed at rather than
 *  restated. */
function whatSurvived(f) {
  const bits = [];
  if (f.p10OverIc50 !== null && f.medianOverIc50 !== null) {
    bits.push(
      `the SPREAD across the tissue — the tenth-percentile tumour cell sits at `
      + `${num(f.p10OverIc50, 3)} of IC50 and the median at `
      + `${num(f.medianOverIc50, 3)} — which is a ratio taken inside one map and `
      + 'so carries neither the dose nor the potency (the argument is under '
      + 'Delivery or dose)',
    );
  }
  if (f.plasmaBoundHolds === false && f.coverage1x !== null) {
    bits.push(
      `that nothing crossed IC50 anyway — coverage at 1x came back `
      + `${fmtFrac(f.coverage1x)} — under a source term that has no ceiling and `
      + (f.tissueOverPlasmaPeak !== null && f.tissueOverPlasmaPeak > 1
        ? 'was demonstrably over the plasma peak on this run'
        : 'was free to exceed the plasma peak'),
    );
  }
  if (!bits.length) return null;
  return [
    h('b', 'WHAT THE MAP ADDED, SEPARATELY. '),
    'Two things on this run are not in the ratio above: ' + bits.join('; ') + '. ',
    'Neither of those is "the entered dose was low". Everything about the '
    + 'LEVEL of this verdict is downstream of a ratio that was fixed before the '
    + 'tissue was read; what is not downstream of it is how unevenly the '
    + 'compound is distributed across this section.',
  ];
}

/**
 * The disclosure block.
 *
 * @param {object|null} run
 * @param {{compact?: boolean, pointTo?: string}} [opts]
 *   `compact` drops the separated "what the map added" list and the coverage
 *   arithmetic, keeping the ratio, the input statement and the law clause —
 *   for the second panel on the same rail, where the full version is above.
 *   Nothing that is only in the compact form is dropped from the full one.
 * @returns {Node[]} zero or more nodes; append them all.
 */
export function doseInputBlock(run, opts) {
  const f = doseFraming(run);
  if (!f) return [];
  const { compact = false, pointTo = null } = opts || {};
  const under = f.ratio < 1;

  const out = [ratioStrip(f, !compact)];

  if (!under) {
    /* AT OR ABOVE 1. The complement is stated just as plainly, and it is the
     * branch that makes this block informative rather than merely cautionary:
     * here a dose-limited verdict is NOT forced, so it means more. Saying
     * nothing in this branch would make the block read as a standing apology
     * for the shipped default rather than as a measurement of one run. */
    out.push(h('div.caveat.info', h('span',
      h('b', `THE ENTERED DOSE DOES NOT FORCE THIS VERDICT: ${ratioText(f)}. `),
      'The peak plasma concentration entered for this run is at or above the '
      + 'concentration that produces half this compound\'s effect, so nothing '
      + 'in the RunSpec prevents a cell from reaching IC50 and no verdict on '
      + 'this rail is arithmetic at this ratio. Any shortfall reported here is '
      + 'therefore something the transport map produced. The ratio is stated '
      + 'anyway, in both directions, so that its absence is never the reason it '
      + 'went unchecked.')));
    return out;
  }

  /* BELOW 1. */
  const body = [
    h('b', `PART OF THIS VERDICT IS AN INPUT: ${ratioText(f)}. `),
    'The plasma peak is a field of the RunSpec on the left of this app. The '
    + `IC50 is not — it is ${ic50Provenance(f)}. Their ratio is below 1: the `
    + `entered peak is ${num(f.foldUnder, 1)}x under the concentration that `
    + 'produces half this compound\'s effect, and it was fixed before this '
    + 'tissue was read. ',
    ...lawClause(f),
  ];
  if (pointTo) body.push(' ' + pointTo);
  out.push(h('div.caveat.crit', h('span', ...body)));

  if (!compact) {
    const surv = whatSurvived(f);
    if (surv) out.push(h('div.caveat', h('span', ...surv)));
  }
  return out;
}

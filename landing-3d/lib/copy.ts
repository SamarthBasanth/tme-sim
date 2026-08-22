/**
 * Narrative copy, shared between the cinematic scroll experience and the
 * prefers-reduced-motion fallback so the two never drift out of sync.
 *
 * Positioning guardrails (these are product constraints, not style
 * preferences). Codon is a mechanistic drug-transport simulator for
 * research users. Copy here must NOT claim or imply: clinical response or
 * patient outcome prediction, any diagnostic / treatment-planning /
 * surgical use, AI or machine learning of any kind, "digital twin", or
 * personalised medicine for an individual patient. Absolute kill fraction
 * is not a defensible output. Every number that appears below must be
 * traceable to a real run or a real dataset count.
 *
 * DATA PROVENANCE, AND WHOSE PROPERTY A PRIVACY CLAIM IS ABOUT. Copy here
 * must not invite, imply or accommodate identifiable patient data: Codon is
 * scoped to de-identified tissue (docs/DEPLOYMENT.md §0), and there is no
 * tissue upload to describe — sections are registered by directory path, so
 * the exposed surface is the path, id and label, and the label becomes the
 * <h1> of every exported report. Never write "upload".
 *
 * The rule that the rest of this paragraph exists to enforce: ANY
 * data-locality or privacy claim must be traceable to a verified property of
 * THE SOFTWARE, never to the customer's network. Ours is the only one of the
 * two we can inspect, and a promise about a perimeter we cannot see is not a
 * promise, it is a guess printed in the imperative.
 *
 * The claim that failed this test, kept here so it cannot drift back: beat 06
 * read "proprietary tissue and compound data never cross your firewall". That
 * is a statement about the customer's perimeter, and it was false in a
 * specific way — the API shipped wildcard CORS with no default auth, so any
 * page the operator browsed could read the compound library, and the data did
 * have a path out. Note that the failure was not caught by reading the
 * sentence; it was caught by auditing the software the sentence described.
 *
 * What is sayable instead, because it was checked rather than assumed: the
 * solve runs as a single local process and sends no tissue or compound data
 * anywhere (`app/server/` imports no HTTP client, no socket library and no
 * analytics or crash-reporting SDK; this landing page self-hosts its fonts
 * and makes zero third-party requests at runtime). Do NOT upgrade this to the
 * broader "makes no outbound network requests" — that is FALSE of the package
 * as a whole and falsifiable by grep: tmesim/data/{acquire,spatial_ingest}.py
 * and tmesim/data/real/mgmt_spatial.py fetch public reference datasets from
 * GDC, 10x, Sanger, the Allen Institute and NCBI GEO. Those are deliberate,
 * operator-invoked downloads of public data that carry nothing of the
 * customer's, which is exactly why the honest claim is about what the
 * software SENDS and not about how many sockets it opens.
 *
 * A SCORE MUST NAME WHAT THE MODEL WAS GIVEN, AND MUST NOT BE PRINTED AS A
 * POINT WHEN WHAT EXISTS IS AN INTERVAL. Two rules, one paragraph, because the
 * same number broke both and they fail the same way — by looking like more
 * than they are.
 *
 * (i) NAME THE INPUTS. Every performance figure has two halves and the copy has
 * historically printed one: the data it was tested on, and the inputs it was
 * handed to make the prediction. A number earned on descriptors, fields or
 * parameters supplied by the source of the measurement is a claim about THE
 * PHYSICS. A number earned on inputs the product computes for itself is a
 * claim about THE PRODUCT. They are different claims, they can differ by a
 * factor, and a reader has no way to tell which they are looking at unless the
 * sentence says. "Blind" does not say. It answers whether the targets were
 * seen, and is silent on where the inputs came from. So answer it in the copy,
 * in the same breath as the number: what was the model GIVEN?
 *
 * (ii) AND CHECK THE INPUTS WERE THE RIGHT QUANTITY. Naming them is not enough,
 * because the supplier's column can be a different quantity from the one the
 * model consumes — same name, different convention — and then the score is not
 * a determinate number at all. When recomputing it on the right quantity leaves
 * a RANGE, the range is the finding and the range is what ships. Do NOT print
 * the old point estimate with a caveat attached, and do NOT substitute a new
 * point estimate: quote the interval, say what would narrow it, and say plainly
 * that an interval that wide validates nothing. A page that publishes its own
 * refutations does not get to round one off.
 *
 * The claim that failed both halves, kept here so it cannot drift back: beat 05
 * read "R² = 0.718 blind on 71 measured compounds". Every word was defensible
 * alone — the law was never fitted to those diffusivities, so "blind" is
 * earned. But (i) 0.718 was produced by handing the law Pruijn's own measured
 * logP₇.₄ and hydrogen-bond-donor columns as inputs
 * (experiments/two_route_transport.py:246) while the shipping engine computes
 * both itself from a structure (tmesim/molecular/compile_drug.py:190), so the
 * sentence licensed "hand it a compound, get 0.718"; and (ii) the donor column
 * it was handed counts donor HYDROGENS while the law's coefficient is defined
 * per donor HEAVY ATOM, so 0.718 is not a value the model's own convention
 * produces. Nothing in the sentence was false. Note the failure mode both
 * times: it was not caught by reading the sentence. It was caught by reading
 * the run.
 * See results/descriptor_falsification.md, results/hbd_convention_audit.md.
 *
 * AN INVARIANCE CLAIM MUST NAME THE PERTURBATIONS THAT WERE ACTUALLY SWEPT.
 * This rule is written down because guardrail 1 below went through it twice.
 * "Invariant to A, B and C" is read by everyone who does not check which
 * letters are listed as "invariant" — the list is heard as an illustration of
 * a general property rather than as its full extent. So the perturbation set
 * is not a detail that belongs in a footnote or in a linked results file: it
 * belongs in the same sentence as the word, and a sentence that cannot carry
 * it does not get to use the word. The corollary is the load-bearing half: an
 * axis that has NOT been swept must never be implied to be inside the list by
 * the sentence's silence about it. Tissue capacity was outside the list for
 * as long as the claim existed, and when it was finally swept it broke the
 * claim — which is exactly the outcome an unnamed axis is most likely to
 * produce, since nobody chose it in the first place.
 *
 * AN ILLUSTRATIVE EXAMPLE MUST BE CHECKED AGAINST THE SAME DEFENCES AS THE
 * HEADLINE VERDICT. Every row on this page has two halves — the registered
 * verdict, and the vivid particular that makes it land — and only the first
 * has ever been audited, on the reasoning that the example is not the claim.
 * That is backwards. The example is the part a reader retains, so it is the
 * part most likely to be repeated, and it must be at least as well supported
 * as the sentence it is decorating. Run the standing objections at it before
 * it ships, not only at the verdict.
 *
 * The example that failed this test, kept here so the failure mode cannot
 * drift back: beat 05's ranking row illustrated a capacity-driven ranking
 * failure with a named compound swap (`vinorelbine` -> `doxorubicin` on the
 * breast panel). The standing objection to the whole result was that the
 * capacity rule partly re-applies a per-compound ordering the engine already
 * carries. Run against the VERDICT that objection fails, and fails on its own
 * terms: the arms that break the order are uniform, they have no per-compound
 * content, and nothing can be double-counted out of them. Run against the
 * EXAMPLE it lands squarely — the correlation is +0.780 on that panel
 * specifically, and residualising restores `vinorelbine`. So the verdict
 * survived scrutiny that its own illustration did not, and every reader who
 * carried away the compound pair rather than the threshold carried away the
 * weaker half of the row. The objection was not even new: it was already on
 * file (capacity_ranking_robustness.md sec. 6.2 raised it and explicitly left
 * it open). Nobody had thought to point it at the example.
 *
 * TWO OPERATIONAL RULES FALL OUT OF THAT, AND THE SECOND IS THE ONE THAT GETS
 * SKIPPED. First: prefer the illustration the objection cannot reach. A
 * threshold from an arm with no per-compound content (36.3, 5.28) is not
 * merely safer than a pair of drug names, it is a sharper fact — it says the
 * order breaks BELOW the constant the measurement itself recommends, which no
 * compound pair can say. Second, and this is the half a defensive page drops:
 * a correction that only ever runs downhill is not honesty, it is flinching.
 * When a retraction frees room, the room goes to the disclosure that costs us
 * most, not back into prose — here it bought the rat-brain-slice provenance
 * of the swept range, which makes that range look LESS transferable to the
 * human panels the row is about. See results/uptake_proxy_double_count.md.
 *
 * Two guardrails narrowed by measurement, both of which read as weaker than
 * the obvious marketing line and must not drift back to it:
 *
 * 1. RANKING IS NOT UNCONDITIONALLY INVARIANT, AND AS OF 2026-08-21 IT IS NOT
 *    ROBUST EITHER. This guardrail has been narrowed twice and the copy must
 *    not drift back to EITHER of the wider versions.
 *    WAS (still true, still say it): ordering is invariant to common-mode
 *    perturbations — cell-state route, oxygen route, alpha, grid — and NOT to
 *    drug-specific transport parameters; correcting the transvascular delivery
 *    law changes the top-1 GBM pick.
 *    See results/transvascular_exchange_revalidation.md.
 *    IS: a pre-registered sweep of the one term that was never in that list
 *    has now been run, and the FULL ORDERING does not survive it.
 *    `s_eff`, the tissue's capacity to store a compound, is 1.000000 for every
 *    compound in the shipping engine. Measured rat brain-slice `Vu,brain`
 *    spans 3.2-4184 mL/g (n = 33) — a spread of 1307.5x — with a
 *    GMFE-minimising single constant of 39.
 *    THE DEFENCE THAT FAILS, AND IT FAILS ON ITS OWN TERMS: "a constant error
 *    cancels out of a ranking" is exactly the UNIFORM case, and under a
 *    uniform capacity the published order breaks at R = 36.3 on the 10-drug
 *    GBM panel and R = 5.28 on the 10-drug breast panel — BOTH below the best
 *    available constant of 39, and the breast figure only 1.7x above the
 *    bottom of the measured range. The mechanism is not a coincidence and is
 *    worth knowing before editing this: attenuation depends on k*T/R, so a
 *    uniform R still reorders a panel whose drugs differ in k by 27x.
 *    THE NAMED TOP-PICK SWAP WAS RETRACTED ON 2026-08-21 AND MUST NOT COME
 *    BACK IN ANY FORM. It read: "AND THE TOP PICK MOVES TOO, on a NAMED
 *    assignment rather than a sampling tail: on the breast panel this repo's
 *    own structural K_b rule (two_route.binding_capacity_from_logp) replaces
 *    `vinorelbine` with `doxorubicin`, as does its access-limited variant."
 *    WHY IT WENT: the engine ALREADY carries a per-compound quantity that its
 *    own docstring calls a capacity ordering — `validation.calibration.
 *    uptake_proxy`, which `alpha` multiplies into every compiled drug's
 *    `binding` — and on the BREAST panel specifically it is strongly
 *    correlated with the structural K_b rule: Spearman +0.780 (p = 0.008),
 *    Pearson on logs +0.741 (p = 0.014), R^2 0.549. So applying that rule on
 *    top of the shipping engine partly applied the same per-compound ordering
 *    twice. Strip the component the proxy already predicts and `vinorelbine`
 *    comes back as breast top-1 (panel spread 874.7x -> 141.1x), and the GBM
 *    10-drug FULL published order comes back too. Both named findings are
 *    substantially double-counted and are discounted. Note the earlier
 *    Spearman of +0.333 (p = 0.35) quoted in capacity_ranking_robustness.md
 *    sec. 6.2 was the GBM panel ONLY, and quoting it alone understated the
 *    overlap exactly where the swap happened.
 *    WHY THE VERDICT IS UNCHANGED, AND THIS IS THE HALF THAT MUST BE CARRIED
 *    WITH EQUAL FORCE — the retraction reaches the illustration and nothing
 *    else:
 *      (i) THE UNIFORM ARMS HAVE NO PER-COMPOUND CONTENT AT ALL, so no overlap
 *          of any strength can shrink them. 36.3 and 5.28 stand exactly as
 *          written above, and they are what FULL ORDER NOT ROBUST rests on.
 *     (ii) PUSHED TO ITS LIMIT THE DEFENCE MAKES THE RESULT WORSE. Give every
 *          compound a capacity exactly proportional to `uptake_proxy` — the
 *          maximal double-count, where none of the per-compound shape can be
 *          new physics — and the published order breaks at a geometric-mean
 *          capacity of 14.7 on GBM 10-drug (against 36.3 uniform) and 10 on
 *          GBM 5-drug (against 435 uniform). Only breast is more benign (7.57
 *          un-clipped against 5.28). The mechanism: R_i proportional to k_i
 *          makes k*T/R the same for every compound, so attenuation FLATTENS
 *          from its shipping range of 0.9075-0.9965, and flattening reorders.
 *    (iii) THE SAMPLED ARMS BARELY MOVE and nothing crosses any bar.
 *          P(full published order) 0.0047 -> 0.0055, 0.2002 -> 0.2397, 0.0006
 *          -> 0.0000; P(top pick) 0.730 -> 0.765, 0.584 -> 0.566, 0.797 ->
 *          0.792, all three still failing the decision-grade bar. This is the
 *          capacity-side top-1 evidence that SURVIVES, and it is weaker than
 *          the named swap was — a sampling result, not a named assignment —
 *          which is why beat 05 leans on the delivery-law top-1 change
 *          instead and does not print these three numbers.
 *    ALSO ESTABLISHED THERE, AND IT IS NOT A POINT IN THE ENGINE'S FAVOUR:
 *    `uptake_proxy` reaches the solver at one place, `physics/diffusion.py:315
 *    k = clearance + binding`, as an irreversible first-order SINK with no
 *    bound pool. Measured on the shipped solver, x10 on `binding` moves the
 *    steady state to 0.1089x while a capacity of 10 leaves it at
 *    0.999999999998x. They are different operators, so the double-count does
 *    not net out; the engine is using a capacity ORDERING to set a quantity
 *    (penetration length) a real capacity cannot touch, while the real
 *    capacity slot `s_eff` stays pinned at 1.0.
 *    NOT CLAIMABLE FROM ANY OF IT: that `uptake_proxy` agrees, or disagrees,
 *    with measured tissue capacity. Exactly 1 of the 35 rat brain-slice
 *    compounds is on either shipping panel (`docetaxel`) and the frozen
 *    dataset has no SMILES column, so n = 1 and there is no test. The two
 *    access-limited variants were NOT re-tested.
 *    See results/uptake_proxy_double_count.md sec. 0, 3.1, 3.2, 3.3, 4.1.
 *    WHAT CAN NO LONGER BE WRITTEN IN ANY FORM: "the compound ranking is
 *    robust"; "the absolute values may be off but the ranking holds"; "our top
 *    pick is robust" without the qualification below. The middle one is the
 *    dangerous one — it is the sentence a reader constructs for themselves out
 *    of beats 02 and 05 if the copy does not close it off.
 *    WHAT MAY BE WRITTEN, sec. 9's condition, and do NOT invent a weaker or a
 *    stronger version of it: "the predicted ordering is stable against a
 *    uniform tissue-capacity term up to R ~ 5.28 for the full order and
 *    R ~ 435 for the top-1 pick. The engine assumes R = 1. Measured
 *    `Vu,brain` spans 3.2-4184 with a best constant of 39." The stable regime
 *    does not contain the measured range, and for the full order it does not
 *    contain even 39.
 *    WHAT SURVIVES, AND MUST BE CARRIED WITH EQUAL FORCE — the page is not
 *    permitted to be more pessimistic than the evidence either:
 *      (a) STEADY-STATE AND EQUILIBRIUM-REFERENCED OUTPUTS ARE EXACTLY
 *          CAPACITY-INVARIANT, measured to 1e-9 (sec. 2b): penetration length,
 *          interstitial:plasma ratio, the shape of a converged gradient. R
 *          divides D, v, S and k alike, so C_ss = S/k cannot move. This
 *          failure mode does not reach them AT ALL, it is specific to
 *          finite-window AUC, and it is now the strongest demonstrated
 *          statement on the page. It lives in beat 02, whose own
 *          shape-versus-scale argument it is a measured instance of.
 *      (b) On the five-drug GBM sub-panel that results/ranking_robustness.md
 *          actually scores, EVERY named assignment — uniform R = 39, the
 *          structural rule, the access-limited rule — reproduces the published
 *          order exactly.
 *      (c) The GBM 10-drug top-1 (`vincristine`) never changes anywhere on the
 *          swept range under a uniform capacity.
 *    THREE FIGURES THAT WERE CIRCULATED WRONG AND MUST NOT COME BACK: the
 *    measured spread is 1307.5x, NOT 631x; both shipping panels run a 6-HOUR
 *    window (21,600 s), not 12 h — PDEConfig's 43,200 s DEFAULT is not what
 *    either panel constructs, and a shorter window makes capacity matter MORE;
 *    and 0.730 / 0.584 / 0.797 are P(the TOP PICK holds) under per-compound
 *    sampling, NOT P(the full order holds), which is 0.005 / 0.200 / 0.001.
 *    All three of the top-pick figures fail this project's own decision-grade
 *    bar (leader's 90 % interval entirely above 0.90).
 *    NOT claimable from any of this: that a capacity term BELONGS in the
 *    engine. results/intracellular_retardation.md refuted the only available
 *    per-compound rule on measured Pruijn HT-29 rank agreement (77.4 % ->
 *    51.6 %). The finding is not "the model is missing a term", it is "the
 *    ordering is conditional on a term whose value is unknown, and setting it
 *    to 1.0 is a load-bearing choice rather than a neutral default".
 *    See results/capacity_ranking_robustness.md sec. 9 and sec. 10.
 * 1b. THE DIFFUSIVITY LAW NOW USES STRUCTURE, AND THE COPY MUST NOT OVERSTATE
 *    HOW WELL. This guardrail CHANGED on 2026-08-13 and the previous version is
 *    kept here because the copy must not drift back to EITHER extreme.
 *    WAS: the only molecular input to D_eff was molecular weight; logP, HBD,
 *    HBA, TPSA and charge were computed and discarded; against 71 compounds
 *    measured in one cell line the measurement spanned 29.7x where the engine
 *    spanned 2.4x (R^2 = 0.022), and at identical MW the engine was exactly
 *    1.00x against a measured 3.7x. "descriptors -> transport" was therefore
 *    measured and REFUTED.
 *    IS: the shipping model is two-route (paracellular + transcellular) and
 *    consumes logP and hydrogen-bond donor count. Unfitted, it scores an R^2
 *    that lies somewhere in 0.498-0.839 (see 1c -- 0.718 is one point inside
 *    that interval, not the answer), 99 % within 3x as published, and orders 24
 *    of 31 exact-mass pairs correctly against the old model's 0 of 31. So
 *    "a compound's structure drives its transport" is STILL SAYABLE, and note
 *    WHY it survives: the ENTIRE interval, low end included, is more than an
 *    order of magnitude above the 0.022 the MW-only model scored. The interval
 *    is too wide to validate the model and far too high to refute it. Say it in
 *    that form, about the physics, and no further.
 *    THE LIMIT THAT REPLACES THE OLD ONE, AND IT IS NOT OPTIONAL: the model
 *    reproduces 8.3x of a measured 29.7x spread -- about a QUARTER -- so it
 *    UNDER-SEPARATES similar-mass analogues rather than resolving them, and the
 *    whole test is ONE cell line. (That ratio comes from the same run as 0.718
 *    and inherits its input defect; it is quotable because it is an UNDER-claim
 *    -- across the searched space the spread stays ~6.4-8.3x, i.e. ~21-28 %.)
 *    Do NOT write "validated", "accurate", "gets transport right", or any
 *    absolute diffusivity figure. Say that structure reaches transport, say
 *    what the model was GIVEN, quote the INTERVAL rather than 0.718 alone, and
 *    say it recovers about a quarter of the measured spread.
 *    "n = 71" is the DATASET; "n = 68" is what 0.718 was computed over.
 *    See results/two_route_promotion.md, results/diffusivity_falsification.md.
 * 1c. 0.718 IS NOT A DETERMINATE NUMBER AND MUST NEVER SHIP AS A POINT
 *    ESTIMATE. THIS IS THE STRONGEST OF THE THREE AND IT IS THE ONE MOST
 *    LIKELY TO BE ROUNDED OFF BY A FUTURE EDIT.
 *    WHAT HAPPENED: 0.718 was scored by feeding the law Pruijn's `HD` column,
 *    which counts donor HYDROGENS (settled from the paper's own Methods
 *    sentence, whose checkable other half reproduces 5/5). The law's
 *    HBD_LOGP_DECREMENT = 1.34 is defined per donor HEAVY ATOM -- the project's
 *    own frozen pre-registration says so in as many words
 *    (results/sweep/tworoute_promotion_prereg.json), and the measured Abraham A
 *    values agree (heavy-atom convention RMSE 0.90 and unbiased; hydrogen
 *    convention 1.43 with a +0.75 systematic over-penalty).
 *    THE DIRECTION MATTERS AND IS EASY TO GET BACKWARDS: THE ENGINE IS RIGHT,
 *    THE VALIDATION WAS WRONG. Do NOT write this as a bug in the product, and
 *    do NOT let anyone "fix" the engine to match -- switching to NHOHCount
 *    would INTRODUCE the dimensional error, and would move calibrate_alpha by
 *    32 % (22.2235 -> 15.1009) for nothing.
 *    THE INTERVAL, WHICH IS THE ONLY HONEST FORM OF THE NUMBER: recomputed on
 *    the engine's own convention, 45 of the 71 rows pin EXACTLY with no
 *    structure needed (42 because HD <= 1 forces the two conventions to agree,
 *    plus 5 with resolvable structures, 2 overlapping). The other 26 are free,
 *    and searching that 6.7e7-configuration space puts R^2 anywhere in
 *    0.4978-0.8388. That interval CONTAINS 0.718, so the published figure is
 *    NEITHER REFUTED NOR CONFIRMED, and an interval 0.34 wide validates
 *    nothing. Quote the interval. Do NOT quote 0.718 with a footnote, and do
 *    NOT substitute a replacement point estimate -- including the "likeliest"
 *    all-free-rows-carry-an-NH2 case at ~0.638, which is a scaffold argument
 *    and is not measured. R^2 is NOT MONOTONE in the per-row shift, so no
 *    single scenario is a bound and the older "-1 -> 0.569" bracket is not one
 *    either; it silently shifts the 42 rows that provably cannot shift.
 *    UNRESOLVED, NOT UNKNOWABLE, AND THE COPY SHOULD READ THAT WAY: the 26
 *    structures exist. The CSV's `ref_footnote` column points at the synthesis
 *    papers and the corresponding author is a named route
 *    (data_cache/pruijn2008/PROVENANCE.md). Not by OCR-ing the Table S1
 *    drawings -- that is inventing inputs.
 *    See results/hbd_convention_audit.md.
 * 1d. SEPARATELY, THE DESCRIPTOR PATH A USER ACTUALLY EXERCISES HAS NEVER BEEN
 *    SCORED AT ALL. This is a DIFFERENT defect from 1c and about a DIFFERENT
 *    descriptor. Even set the donor convention aside: 0.718 was also handed
 *    Pruijn's MEASURED logP7.4 column, while the engine computes Crippen
 *    MolLogP from a structure (tmesim/molecular/compile_drug.py:190). No score
 *    exists for that path, because 66 of the 71 compounds are Auckland SN code
 *    numbers with no structure resolvable in PubChem. What was measured, on the
 *    five that resolve: RMS 1.62 log units of logP error, +3.28 on the one
 *    ionisable base, one-signed, because ionisation is modelled NOWHERE on the
 *    transport path (Crippen is a neutral-species logP; the paper's column is a
 *    pH-7.4 distribution coefficient).
 *    THE FIGURE THAT MUST NOT BE PRINTED AS A RESULT: propagating that error
 *    projects R^2 to roughly 0.22-0.50, central ~0.27. It is a projection off a
 *    FIVE-COMPOUND estimate, it is the weakest number in that file, and the
 *    weakness IS the finding -- so do NOT write 0.27, or any point value from
 *    that band. Silence beats false precision: "has never been scored" is the
 *    whole of what is known. (Its low edge additionally assumed the engine's
 *    donor count was the wrong quantity, which 1c refutes.) Equally, do NOT let
 *    the copy imply the physics collapses: even that projection sits well above
 *    0.022, so the transcellular physics is MASKED by descriptor error, not
 *    erased by it, and the failures are separable.
 *    See results/descriptor_falsification.md.
 * 2. THE SPATIAL SPREAD MUST CARRY ITS HONEST NULL. The 2.0x is real and
 *    reproducible, but a well-mixed model's 1.0 is not the comparator that
 *    bounds it: a null that preserves the macrophages' own clustering
 *    returns ~1.7, so the part attributable to where these macrophages
 *    actually are is ~1.08-1.21x (1.081 against the 400 um block shuffle,
 *    1.208 against the clumpiness-matched permutation). Quote both numbers
 *    or neither. Carry the confound too: distance to the tissue edge alone,
 *    with no immune input whatsoever, already reproduces 82% of the raw
 *    2.0x, and it is the most deflating single fact about this claim. Do
 *    NOT cite monotonicity across quintiles as corroboration — every control
 *    reproduces it, so it carries no information. The statistic is invalid
 *    on spot data, where immune positions are synthesised rather than
 *    measured. See results/immune_sink_impact.md.
 */

/**
 * One row of an in-flow instrument ledger: a named quantity, the STATE the
 * project is registered as holding on it, and the measurement that put it
 * there.
 *
 * This shape exists because two beats needed the same device and were both
 * expressing it as prose. `state` is deliberately a short mono string rather
 * than an enum with a colour: the page encodes limit severity typographically
 * (mono, uppercase, hairline rule) and NOT with a red/amber/green vocabulary,
 * because a caveat rendered as a warning chip reads as an exception to a
 * working product, and here the caveats ARE the product surface. No new hue
 * enters the palette for this.
 */
export interface LedgerRow {
  label: string;
  /** Short registered verdict, e.g. "SOLVED", "NOT CALIBRATED". */
  state: string;
  detail: string;
}

export interface SectionLedger {
  caption: string;
  rows: readonly LedgerRow[];
}

export interface SectionCopy {
  kicker: string;
  eyebrow: string;
  headline: string;
  body: string;
  align: "left" | "right" | "center";
  /**
   * Optional structured readout rendered under `body`, in flow, in both the
   * cinematic and reduced-motion experiences. In flow rather than floating:
   * the two floating panels already collide with narrative text at some
   * scroll positions (see globals.css), and a disclosure that can be covered
   * is not a disclosure.
   */
  ledger?: SectionLedger;
}

export const SECTION_COPY: readonly SectionCopy[] = [
  {
    kicker: "01",
    eyebrow: "Real measured tissue",
    headline: "See Where the Drug Actually Goes.",
    body: "The section behind this text is a real one. Codon rebuilds tumor architecture from spatial transcriptomics, then runs transport physics on the result to work out where a drug actually reaches.",
    align: "left",
  },
  {
    kicker: "02",
    eyebrow: "Transport physics",
    headline: "The shape is solved. The scale is not.",
    /*
     * THE LAST SENTENCE IS THE PAGE'S STRONGEST DEMONSTRATED CLAIM AND IT WAS
     * ADDED ON 2026-08-21, THE DAY THE RANKING CLAIM WAS NARROWED. It belongs
     * here rather than in beat 05 for a reason that is not layout: this beat's
     * whole argument is shape-versus-scale, and capacity-invariance is a
     * measured instance of it — the one place where "the shape survives what
     * the scale does not" is an equality rather than an argument. R divides
     * diffusivity, velocity, source and sink alike, so the converged field
     * C_ss = S/k cannot move with it; run past equilibrium at R = 1, 39 and
     * 4184 the mean concentration agrees to 1e-9
     * (results/capacity_ranking_robustness.md sec. 2b).
     *
     * The second half of the sentence is not a hedge and must travel with the
     * first. The same term that cannot touch a converged field DOES move what
     * a finite exposure window accumulates, and the compound ordering is built
     * on exactly that finite window — which is why beat 05's first ledger row
     * now reads NOT ROBUST. Stating the invariance without the exception would
     * hand a reader the sentence guardrail 1 forbids ("the levels are off but
     * the ordering holds"), assembled out of two beats neither of which said
     * it. Do not split these two clauses across beats.
     */
    body: "Convection, diffusion and reaction, solved on that geometry. Interstitial pressure builds where vessels leak and drainage fails and pushes back against the drug coming in; oxygen and hypoxia solve on the same grid, and response follows Hill kinetics from the concentration each cell sees. What that buys is the shape of the field, not its magnitude — and the difference is measured, not asserted. One piece of it is exact rather than argued: how much drug the tissue can store is a number the engine assumes rather than measures, and it cannot move a converged field at all — checked to one part in a billion. It does move what a finite exposure window accumulates, which is where the compound ordering in 05 sits.",
    align: "left",
    /*
     * This beat used to be a four-item feature list with no number in it, and
     * its only instrumentation was a floating callout that named the quantity
     * and nothing else. Every value below is from results/ifp_validation.md,
     * which registers three rows against wick-in-needle measurements.
     *
     * The 41x FAIL is quoted WITH the under-resolution row, because that row
     * was registered specifically so the other two "cannot be quoted without
     * it" (ifp_validation.md:375). Do not drop it to shorten this block.
     *
     * The GBM row is NOT TESTABLE, not a pass: the published intracranial
     * value (2.0 +/- 2.5 mmHg SD, PMID 9062403) is statistically
     * indistinguishable from the same group's normal-tissue value, so it
     * cannot separate this engine from one with no tumour in it. That is a
     * far more useful thing to say than "validated on GBM", and it is the
     * only thing the anchor supports.
     */
    ledger: {
      caption: "Interstitial pressure · registered against wick-in-needle",
      rows: [
        {
          label: "Field shape",
          state: "SOLVED",
          detail:
            "On the measured geometry, not a template: where pressure rises, and the gradient it puts against incoming drug.",
        },
        {
          label: "Absolute scale",
          state: "FAIL · 41× LOW",
          detail:
            "Against primary breast carcinoma at 15 ± 9 mmHg (PMID 1423283) — the one IFP anchor that discriminates, and the engine fails it.",
        },
        {
          label: "Grid resolution",
          state: "UNDER-RESOLVED",
          detail:
            "Refining the grid moves the solve further from the measurement, not closer. Registered so the two rows above are never quoted alone.",
        },
        {
          label: "Glioblastoma",
          state: "NOT TESTABLE",
          detail:
            "The published intracranial value is statistically indistinguishable from normal brain, so it cannot separate this engine from one with no tumour in it.",
        },
      ],
    },
  },
  {
    kicker: "03",
    eyebrow: "Delivery or potency",
    headline: "Which one killed the compound?",
    /*
     * The middle sentence was REPLACED on 2026-08-14 and the previous version
     * is named here so it cannot come back by accident. It read: "Weight is
     * the only property reaching transport, so similar-mass analogues cannot
     * be separated." That was true of the model this page shipped when the
     * sentence was written, and it is now false — see the 1b guardrail above
     * and results/two_route_promotion.md.
     *
     * The replacement is deliberately NOT the triumphant form. "Structure
     * reaches transport" is sayable; "transport is solved" is not. The
     * quarter-of-the-spread limit travels in the SAME sentence as the claim,
     * because that is the clause a reader skips if it is moved to a footnote,
     * and beat 05 repeating it is not a reason to drop it here.
     *
     * The word "blind" was REMOVED on 2026-08-21 and must not come back on its
     * own. It was true and it was the wrong true thing: it certified that the
     * targets were unseen while saying nothing about the inputs, and the score
     * was earned on the paper's measured descriptor columns and not on the ones
     * the engine derives from a structure. The clause that replaces it names
     * the inputs AND says the engine's own path is unscored, because either
     * half alone still reads as "hand it a compound, get 0.718" — the first
     * sounds like a detail, and the second sounds like a to-do. See guardrail
     * 1c, and beat 05's "Close analogues" row, which registers the same
     * condition against the number itself.
     *
     * THIS BEAT DELIBERATELY QUOTES NO SCORE, and that is now load-bearing
     * rather than incidental. Under the model's own donor convention there is
     * no single score to quote — only the interval in guardrail 1c — so the
     * claim here is made in the one form the evidence supports: qualitative,
     * conditioned on its inputs, and bounded by the spread limit two sentences
     * later. If a future edit wants a number in this beat it must bring the
     * whole interval, and this beat does not have the room for it (below).
     *
     * ON LENGTH, MEASURED. This beat is vertically CENTRED, so every line it
     * gains pushes its kicker up half a line into the fixed header — at 390x844
     * the longer first draft of this clause collided with the CODON mark by
     * 16px across four scroll positions. The ceiling here is about 16 lines of
     * body (~770 characters); the version below is inside it with the
     * correction intact. If this beat must say more, the layout gives it the
     * room — do not buy the room by dropping the input condition.
     */
    body: "A compound compiles from its own structure: two routes through the cell layer, paracellular and transcellular, fed by lipophilicity and hydrogen-bond donors as well as mass, with potency from measured GDSC2 screens. Structure reaches transport now — scored unfitted against 68 measured compounds, but on descriptor columns supplied by that same measurement; the ones the engine derives from a structure have never been scored. Even so it recovers only about a quarter of the measured spread, so close analogues are under-separated rather than resolved. Each run then decomposes survival cell by cell: a molecule that never arrives and one that arrives too weak look identical in a survival curve, and need opposite fixes.",
    align: "right",
  },
  {
    kicker: "04",
    eyebrow: "Spatial attribution",
    headline: "The mean is not the finding. The spread is.",
    body: "In one section, tumor cells in the macrophage-densest quintile hold about 0.43 of drug exposure. In the sparsest, about 0.86 — a 2.0x spread inside a single tumor, where a well-mixed compartmental model returns exactly 1.0 by construction. That comparison flatters us. Distance to the tissue edge alone, with no immune input, reproduces 82% of that ratio; relocating the same macrophages at random while keeping their real clustering still returns about 1.7. What these particular macrophages add is 1.08-1.21x — smaller, and measured. A compartment expresses neither.",
    align: "left",
  },
  {
    kicker: "05",
    eyebrow: "Stated limits",
    headline: "The limits ship with it.",
    /*
     * THIS BEAT WAS A 1,180-CHARACTER PARAGRAPH AND IS NOW A LEDGER. Nothing
     * was cut. Every clause of the paragraph is below as its own row with its
     * own number; the restructure exists because a page whose argument is
     * "we publish what breaks" was burying five distinct limits in one block
     * of prose, where the fifth is unreachable. A limit a reader does not
     * finish is not a disclosure.
     *
     * If a future pass shortens this, shorten the DETAIL text. Do not remove
     * a row: each one is a separate registered result, and dropping any of
     * them re-opens a claim the project has already measured shut.
     *
     * WHY "CLOSE ANALOGUES" CARRIES AN INTERVAL IN ITS VERDICT COLUMN. It read
     * "UNDER-SEPARATED" and now reads "0.498–0.839, NOT 0.718", which is the
     * same device the "Spatial spread" row below already uses for the same
     * reason: the honest range, set against the headline number a reader would
     * otherwise carry away. It is deliberate that this is the STATE and not
     * only the detail — the state column is the part that gets read straight
     * down without the prose, and a skimmer who reads "UNDER-SEPARATED" and
     * the number 0.718 has been told the smaller of the two problems. The
     * under-separation verdict is not lost; it was already stated in full in
     * the detail text and still is. See guardrail 1c: 0.718 is one point
     * inside that interval, the interval cannot be narrowed without 26
     * molecular structures we do not currently have, and NO point estimate —
     * not 0.718, not a replacement — may ship in its place.
     *
     * ON LENGTH, MEASURED, BECAUSE THIS BEAT IS THE PAGE'S TIGHTEST. Beat 05
     * is top-aligned inside a pinned viewport that clips at overflow:hidden.
     * At 390x844 its bottom edge is 808.7px against the ~805px this file's
     * responsive block records as the budget for surviving its own 26px fade
     * transition (worst case 832.9 of 844, measured at 240 scroll positions
     * per viewport, zero clipped and zero out-of-viewport at every opacity).
     * That margin is SPENT. The next row or clause added here needs the layout
     * to give it room — it must not be bought back out of a caveat.
     *
     * HOW THE CAPACITY RESULT WAS PAID FOR ON 2026-08-21, since the rule above
     * says where the room may NOT come from. It came from this lead paragraph,
     * which went from three lines to one, and no caveat was touched: the two
     * lines it gave up carried the CLAIM ("ranking on the same tissue is the
     * output to trust — far more robust than any absolute number"), which the
     * capacity sweep refuted. Retracting an overclaim is the one edit that
     * frees vertical space without costing disclosure, and it was the only one
     * available here. Measured at 390x844: a lead line is 20.4px and a ledger
     * detail line is 13.83px, so the two lines surrendered (-40.8px) bought
     * three (+41.5px) and the beat lands 0.7px from where it started, 808.0 ->
     * 808.7. There is no second paragraph to spend, and the next edit that
     * needs a line has nothing left to sell.
     *
     * WHAT DID NOT FIT, RECORDED SO IT IS A KNOWN OMISSION AND NOT AN OVERSIGHT
     * — AND IT IS NOW PAID FOR, 2026-08-21. The 3.2-4184 range in the first row
     * is measured in RAT BRAIN SLICES (`Vu,brain`, n = 33), and applying it to a
     * human breast panel is itself an assumption — results/
     * capacity_ranking_robustness.md sec. 8 says so about its own docetaxel
     * value. Naming that provenance cost one line the row did not have. It was
     * never a self-serving omission: an unqualified "measured values run 3.2 to
     * 4184" makes the swept range look MORE authoritative and therefore the
     * model look WORSE, so the direction was against us, which is exactly why
     * this was registered as the first thing any spare room should buy.
     * Retracting the vinorelbine/doxorubicin example freed it, and the row now
     * reads "measured rat brain-slice values run 3.2 to 4184 (n = 33), assumed
     * to carry to human tissue". Measured at 390x844 the detail is 6 lines
     * before and after (83.0px, section bottom unchanged at 808.7) and at
     * 1512x950 it went 5 lines to 4. The margin is still SPENT: the next
     * addition here has nothing left to sell.
     */
    // The lead is now one line and it is about the METHOD, not about which
    // output to trust. That is a deliberate demotion: after the capacity
    // sweep, every claim on this page about the physics carries a named range
    // or a named exception, and the only thing that holds unconditionally is
    // that the project tests its own claims and prints the ones that break.
    // Saying so directly is more honest than promoting the next-strongest
    // physics claim into the vacancy and having it narrowed in turn.
    body: "We test what we ship, and publish what breaks.",
    align: "right",
    ledger: {
      caption: "Carried by every run · none of this is buried",
      rows: [
        /*
         * THIS ROW READ "CONDITIONAL" UNTIL 2026-08-21 AND THE UPGRADE IS THE
         * POINT, SO IT IS RECORDED HERE. "CONDITIONAL" was true and was
         * understated in one specific way: it named four perturbations the
         * ordering survives and one it does not, while the beat above it still
         * said ordering was the output to trust. A pre-registered sweep of
         * tissue capacity — an axis that had never been in the list at all —
         * then broke the full published order under a UNIFORM value, which is
         * precisely the case the "a constant error cancels out of a ranking"
         * defence covers. So the verdict column now carries the retraction
         * rather than a qualifier, the same way "Close analogues" and "Spatial
         * spread" carry theirs.
         *
         * "FULL ORDER" IS LOAD-BEARING AND IS NOT A SOFTENER. It is the
         * accurate scope: the top pick is stable against a uniform capacity
         * (to R ~ 435, and on the GBM panel across the whole measured range),
         * and on the five-drug GBM sub-panel every named assignment reproduces
         * the published order exactly. Writing a flat "NOT ROBUST" would be
         * more pessimistic than the evidence, which is its own kind of
         * dishonesty on a page whose argument is that it reports what it
         * measures.
         *
         * THE EXAMPLE CHANGED ON 2026-08-21; THE VERDICT DID NOT. The detail
         * used to close its middle sentence with "our own structural rule swaps
         * breast's top pick, vinorelbine for doxorubicin". That swap turned out
         * to be substantially double-counted — on the BREAST panel the
         * structural rule and the per-compound sink the engine already applies
         * are strongly correlated (Spearman +0.780, p = 0.008), and removing
         * what the engine already accounts for puts `vinorelbine` back on top.
         * The named pair is therefore GONE from the page and must not return in
         * any form, qualified or not: it was the vivid half of the row, so a
         * reader carried it away in place of the thresholds, and it was the
         * half that did not survive. See the ILLUSTRATIVE EXAMPLE guardrail at
         * the head of this file, which this row is the case study for.
         *
         * WHY THE ROW IS NOT WEAKER FOR IT, and why nothing else moved. The
         * example was never what the verdict rested on. 36.3 and 5.28 come from
         * a UNIFORM capacity — one value for every compound, so no per-compound
         * ordering to double-count — and no correlation of any strength can
         * reach an arm with no per-compound content. Both are still below the
         * best-fitting constant of 39, which is a sharper fact than a pair of
         * drug names: it says the failure happens at the very value the
         * measurement recommends. The line the retraction freed on desktop went
         * to the rat-brain-slice provenance this file had registered below as
         * the first thing any spare room should buy — a disclosure that makes
         * the swept range look LESS transferable, i.e. the correction did not
         * only run downhill.
         *
         * WHAT THE ROW NOW LEANS ON FOR THE TOP PICK, since the named breast
         * swap is gone and "full order" alone would let a reader keep the top-1
         * claim whole: the final clause, "correcting the delivery law moves the
         * glioblastoma top pick", which is a different result on a different
         * axis (results/transvascular_exchange_revalidation.md) and is untouched
         * by any of this. That clause is now load-bearing — do NOT drop it to
         * shorten this row. Guardrail 1 records the capacity-side top-1
         * evidence that also survives (the sampled arms) and why it did not fit
         * here.
         *
         * See guardrail 1 at the head of this file for sec. 9's condition
         * verbatim, for the three figures that were circulated wrong, and for
         * why none of this licenses adding a capacity term to the engine.
         */
        {
          label: "Ranking on one tissue",
          state: "FULL ORDER NOT ROBUST",
          detail:
            "Tissue capacity is 1.0 for every compound; measured rat brain-slice values run 3.2 to 4184 (n = 33), assumed to carry to human tissue. One uniform value breaks the full order at 36.3 on glioblastoma and 5.28 on breast, both below the best constant of 39. Invariant to cell-state route, oxygen route, alpha and grid; not to drug-specific transport, and correcting the delivery law moves the glioblastoma top pick.",
        },
        {
          label: "Close analogues",
          state: "0.498–0.839, NOT 0.718",
          detail:
            "0.718 was scored on the paper's donor counts, in a convention the model does not use. Redone on the engine's, 45 of 71 rows pin exactly and 26 have no published structure — until those are obtained, R² sits anywhere in that interval. It contains 0.718: nothing refuted, nothing confirmed, and a band that wide settles nothing. It still recovers only ~24% of the measured spread, from one cell line, under-separating a congeneric series ~4-fold. Its predecessor scored 0.022, below the whole interval; we published that refutation, then its replacement, and now this.",
        },
        {
          label: "Spatial spread",
          state: "1.08–1.21×, NOT 2.0×",
          detail:
            "A well-mixed model's 1.0 is not what bounds the raw 2.0×. A null preserving the macrophages' own clustering returns ~1.7, and tissue-edge distance alone — no immune input — reproduces 82% of it.",
        },
        {
          label: "Spot platforms",
          state: "NOT REPORTED",
          detail:
            "The spread needs measured cell positions. Spot data estimates them, so the number computes but cannot be attributed — which is the entire claim.",
        },
        {
          label: "Absolute kill fraction",
          state: "NOT DEFENSIBLE",
          detail:
            "Not on any grid, not on any tissue. It is not an output, and it is not on this page.",
        },
        {
          label: "Patient response",
          state: "OUT OF SCOPE",
          detail:
            "Not predicted. Not for diagnostic or treatment-planning use. Null results published alongside positive ones.",
        },
      ],
    },
  },
  {
    kicker: "06",
    eyebrow: "Get started",
    headline: "Run it on your own tissue.",
    // The previous sentence promised "proprietary tissue and compound data
    // never cross your firewall". That described the CUSTOMER'S PERIMETER,
    // which is not ours to make promises about, and it was concretely false
    // while the API shipped wildcard CORS with no default auth. What replaces
    // it claims only verified properties of the SOFTWARE. See the data
    // provenance guardrail at the head of this file before editing this line.
    //
    // "GET BACK A RANKING" WAS RE-READ AGAINST GUARDRAIL 1 ON 2026-08-21 AND
    // DELIBERATELY LEFT STANDING. It describes what the pipeline emits, which
    // is unchanged; it does not say the ordering is robust, and the clause that
    // follows attaches the list — now headed by FULL ORDER NOT ROBUST — to it.
    // The line to watch is "the same list, on your tissue, with your numbers in
    // it": that is a promise the limits TRAVEL with the output, so it is the
    // sentence that makes this beat honest rather than the one that needs
    // fixing. Do not promote the ranking to "the output to trust" here to
    // recover the phrase beat 05 gave up.
    body: "Codon runs where your tissue already is: the solve is a single process on your own machine, and it sends no tissue or compound data anywhere. De-identified tissue only. Bring a section and a candidate series, and get back a ranking, a per-cell mechanism attribution of what killed each compound, and the list above attached to both — the same list, on your tissue, with your numbers in it.",
    align: "center",
  },
] as const;

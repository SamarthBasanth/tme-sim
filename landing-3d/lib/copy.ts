/**
 * Narrative copy, shared between the cinematic scroll experience and the
 * prefers-reduced-motion fallback so the two never drift out of sync.
 *
 * Positioning guardrails (these are product constraints, not style
 * preferences). Vivum is a mechanistic drug-transport simulator for
 * research users. Copy here must NOT claim or imply: clinical response or
 * patient outcome prediction, any diagnostic / treatment-planning /
 * surgical use, AI or machine learning of any kind, "digital twin", or
 * personalised medicine for an individual patient. Absolute kill fraction
 * is not a defensible output. Every number that appears below must be
 * traceable to a real run or a real dataset count.
 *
 * DATA PROVENANCE, AND WHOSE PROPERTY A PRIVACY CLAIM IS ABOUT. Copy here
 * must not invite, imply or accommodate identifiable patient data: Vivum is
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
 * THE PRODUCT NO LONGER RANKS COMPOUNDS, AND THE PROHIBITION IS NOW ABSOLUTE
 * RATHER THAN CONDITIONAL. As of 23e51f6 (2026-08-22) ranking is gone from the
 * CLI, the API, the Studio, the protocol language and the reports — the FEATURE
 * was removed, not just the claim, on the reasoning that as long as the product
 * printed a ranked list that is what people would use, and that is the output
 * that was wrong. Guardrail 1 below used to police HOW an ordering could be
 * described. There is now no permitted description. No ordering, ranked list,
 * shortlist, "top pick", "which compound wins", pairwise call or cross-compound
 * comparison may appear on this page in any form, hedged or not, positive or
 * negative — including as a limit. A ledger row registering that an ordering is
 * not robust is not a disclosure once the ordering does not exist; it is a
 * fossil, and it invites a reader to believe the feature is still there and
 * merely caveated. Guardrail 1's body is kept below UNCHANGED as the record of
 * why the claim died, and it is now history rather than an editing rule: read
 * it before touching the capacity row, do not restore any sentence from it.
 *
 * WHAT REPLACED IT, AND WHY IT IS A DIFFERENT KIND OF CLAIM. One compound, one
 * real section, one question: is the failure a DELIVERY problem or a POTENCY
 * problem, answered per cell. On the flagship specimen the per-cell attribution
 * returns 99.8 % DOSE-limited (results/synthesised_positions.md sec. 5.2), and
 * the measured falloff decays only 1.00 -> 0.50 over 145 um at dx = 12 um (same
 * sentence, same file) — the drug gets almost everywhere and is almost
 * everywhere sub-therapeutic. The reason this is worth the page's lead position
 * is not that it is impressive; it is that it is REFUTABLE. "This compound
 * reached the tissue and was still short of IC50" is checkable against a bench
 * measurement on the same section. "Cisplatin outranks temozolomide" was only
 * ever checkable against an outcome nobody has. The page's honest through-line
 * is therefore not "we publish what breaks" on its own — that is the support,
 * not the thesis. It is that the claim was made SMALL ENOUGH TO CHECK, and the
 * four narrowings that got it there are the evidence that it will be checked.
 *
 * A MECHANISM VERDICT MUST NAME THE LAW IT WAS PRODUCED UNDER. This is the
 * newest guardrail and it is the one that can INVERT the page's lead rather
 * than merely bound it. Every exposure map and every per-cell verdict this
 * product emits is produced under `PDEConfig.transvascular_exchange = False`,
 * the one-way vascular leak `S = PS*Cp`, which has no path from interstitium
 * back to blood and therefore a steady state unbounded above the plasma level —
 * something passive delivery cannot do. Against human microdialysis it misses
 * in OPPOSITE directions (5-fluorouracil 12.8x too high, methotrexate 10.5x too
 * low; results/breast_kp_validation.md), so no single constant absorbs both,
 * and correcting it is not a rescale: it re-weights delivery by permeability
 * and changes the SHAPE of the field, which is the one thing a within-map ratio
 * does not cancel (results/transvascular_exchange_revalidation.md sec. 4.1).
 * app/server/delivery_law.py exists to carry that sentence with every run, and
 * the flagship run's own provenance block says the map "MUST NOT be read as
 * law-independent". The page's lead IS a mechanism verdict, so the page owes
 * the same disclosure the product already makes. It is on beat 06.
 *
 * A DELIVERY VERDICT HAS NEVER BEEN DEMONSTRATED NON-DEGENERATELY, AND THE
 * PAGE MUST NOT IMPLY THREE WORKING SETTINGS. It is tempting to say the
 * discriminator has three settings — dose, delivery, resistance — because all
 * three have been produced. Two of them were produced by the model. The third
 * was not: the only DELIVERY verdict obtained comes from a compound carrying
 * `bbb_penetrance = 0`, so the source term is identically zero by
 * parameterisation and the model is correctly propagating something it was
 * TOLD, not discovering a limit (results/synthesised_positions.md sec. 5.2,
 * and the "THIS IS AN INPUT, NOT A FINDING" branch in app/web/js/ui/reach.js).
 * None of the four measured-position sections is known to be poorly perfused.
 * Say two demonstrated settings and one the product can currently only be
 * handed — the qualified version is the one a sceptical reader trusts, and the
 * unqualified version is refuted by a file in this repo.
 *
 * A NUMBER THE PRODUCT COMPUTED ONCE IS NOT A CITATION. Every figure on this
 * page must trace to a results/ file, and "the panel printed it and someone
 * repeated it" does not qualify — it is exactly how a figure drifts a notch
 * more generous than its source. Three figures were proposed for this page and
 * REFUSED on that ground alone, all of them real live measurements: the median
 * tumour cell's exposure shortfall (`mechanism.exposure_shortfall_fold`, ~3.3x
 * on the flagship), the falloff shape's movement under a capacity term (0.3772
 * absolute at R = 39, penetration length collapsing 279.7 -> 76.7 um), and the
 * bit-identity of the map at `bbb_penetrance` 0.1 and 1.0. The first was left
 * out entirely; the other two ship as QUALITATIVE statements, which lose
 * nothing, because "the map cannot see the barrier" and "the shape moves where
 * the ladder does not" are the whole content. If a results/ file later carries
 * them, the numbers may be added — not before.
 *
 * results/mechanism_attribution.md IS STALE AND MUST NOT BE QUOTED FOR NUMBERS.
 * It was last regenerated on 2026-08-08, before the area-conserving vessel
 * source and before two-route transport went on, and both move exposure for
 * every drug. Its sec. 3.2 reports the flagship section as 97.10 % DELIVERY-
 * limited where the current engine reports 99.8 % DOSE-limited — opposite
 * conclusions about what limits the same compound on the same section. The
 * METHOD it describes (the categories, the bit-exact reconciliation, the
 * precedence rule) is unaffected and still citable. The measurements are not.
 *
 * Two guardrails narrowed by measurement, both of which read as weaker than
 * the obvious marketing line and must not drift back to it:
 *
 * 1. HISTORY, NOT A RULE, AS OF 2026-08-22 — SEE THE ABSOLUTE PROHIBITION
 *    ABOVE. RANKING IS NOT UNCONDITIONALLY INVARIANT, AND AS OF 2026-08-21 IT
 *    IS NOT ROBUST EITHER. This guardrail has been narrowed twice and the copy
 *    must not drift back to EITHER of the wider versions.
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
 *    ADDED 2026-09-20, AND IT IS THE HALF THAT WAS MISSING FOR FIVE WEEKS:
 *    THE TERM IS NOT IN THE SHIPPING PRODUCT. `tmesim/molecular/immune_sink.py`
 *    is imported by its own validation module and by nothing else, so the
 *    mechanism these numbers measure is off BY CONSTRUCTION and no run this
 *    product performs contains it. It did not improve the measured human-glioma
 *    anchor, which is why. Any copy about this beat must say so in the same
 *    breath as the number: the finding is a DEMONSTRATION that geometry carries
 *    what a compartment cannot express, and it is not a feature. Nothing on
 *    this page may present it as something the product computes for you, and
 *    the instrument panel beside it may not be labelled as an ordinary run.
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
    body: "The section behind this text is a real one. Vivum rebuilds tumor architecture from spatial transcriptomics, then runs transport physics on the result to work out where a drug actually reaches.",
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
     * a finite exposure window accumulates, and the product's window is six
     * hours — so every number this page quotes is read from the regime where
     * the invariance does NOT hold. Stating the invariance without the
     * exception would hand a reader the sentence the guardrails forbid ("the
     * levels are off but the shape holds"), assembled out of two beats neither
     * of which said it. Do not split these two clauses across beats.
     *
     * THE CLOSING CLAUSE WAS RE-POINTED ON 2026-08-22 AND THE OLD WORDING MUST
     * NOT COME BACK. It read "which is where the compound ordering in 05 sits".
     * The product no longer produces an ordering, so that clause named a thing
     * a reader could not find and implied the feature still existed; and beat
     * 05's first row no longer reads NOT ROBUST because it is no longer about
     * ranking. The exception itself is unchanged and is stated at the same
     * strength — what changed is which output it is pointed at.
     *
     * THE REPLACEMENT DELIBERATELY NAMES NO OUTPUT, and that is not vagueness —
     * it is the clause earning its keep within hours of being written. A
     * cutover of the reach map to a CONVERGED STEADY STATE ("the reach this
     * geometry permits" rather than "the reach a 6 h exposure achieved") landed
     * in 4cab273 while this pass was being measured. The earlier wording, "which
     * is the window every number here is read from", would have been silently
     * falsified by it: the MAP is now a steady object and is exactly
     * capacity-invariant, with the ladder measured to move 0.000e+00 at every R
     * (results/steady_state_reach.md sec. 4). "The one regime the invariance
     * does not cover" survived unchanged, because it names a regime rather than
     * an output.
     *
     * WHAT DID HAVE TO CHANGE, AND THE SPLIT IS NOW THE INTERESTING PART. Beat
     * 05's capacity row had said the falloff shape is not protected the way the
     * equilibrium quantities are. That was true of the 6 h march and is false
     * of the steady map, so it is gone. The map and the verdict are now in
     * DIFFERENT regimes and the row says which: the map is converged and the
     * term cannot reach it; the per-cell verdict in beat 03 is still read off a
     * finite exposure window, and there the term still moves things. Do not
     * collapse those two halves back together in either direction.
     *
     * ONE TRAP RECORDED, because it is the kind that gets quoted at us rather
     * than by us: with convection ON — the shipped setting — the published
     * capacity sweep is not a faithful capacity change at all (it multiplies
     * the Peclet number by R as well), so it produces a large residual against
     * a map that is exactly invariant. That is a statement about the
     * instrument, not about the map, and no version of it belongs in the copy.
     */
    body: "Convection, diffusion and reaction, solved on that geometry. Interstitial pressure builds where vessels leak and drainage fails and pushes back against the drug coming in; oxygen and hypoxia solve on the same grid, and response follows Hill kinetics from the concentration each cell sees. What that buys is the shape of the field, not its magnitude — and the difference is measured, not asserted. One piece of it is exact rather than argued: how much drug the tissue can store is a number the engine assumes rather than measures, and it cannot move a converged field at all — checked to one part in a billion. It does move what a finite exposure window accumulates — the one regime the invariance does not cover.",
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
     * ON LENGTH, MEASURED, AND THE CEILING IS TIGHTER THAN THE OLD NOTE SAID.
     * This beat is vertically CENTRED, so every line it gains pushes its kicker
     * up half a line into the fixed header. The previous note here put the
     * ceiling at "about 16 lines of body (~770 characters)"; re-measured
     * against a production `next start` build at 390x844 over 240 scroll
     * positions, the beat's worst-case top across the fade is 65.5px and the
     * site header's bottom edge is 70.8px, so the collision is ALREADY at zero
     * margin at the previous body length of 721 characters. Treat 721 as the
     * ceiling, not 770. The body below is 737 and was measured not to add a
     * wrapped line at either viewport; anything materially longer will.
     *
     * WHAT THE 2026-08-22 REWRITE DID, AND WHAT IT REFUSED TO PAY WITH. The
     * closing sentence used to be the SETUP with no payoff — "a molecule that
     * never arrives and one that arrives too weak look identical in a survival
     * curve, and need opposite fixes" — which stated why the discrimination
     * matters and never said what it found. It now carries the finding. The
     * room came from the compilation prose in the first two sentences
     * ("paracellular and transcellular" as a gloss on "two routes"; "now";
     * "Even so"), NOT from either protected clause: the input condition and
     * the quarter-of-the-spread limit are both still here, in the same
     * sentences as the claims they bound.
     *
     * THE LAST TWO SENTENCES ARE THE PAGE'S THESIS AND THEY TRAVEL TOGETHER.
     * The first is a positive, falsifiable result; the second is what makes it
     * different in kind from everything this page has retracted. Do not keep
     * the finding and drop the refutability clause — without it the sentence is
     * just a smaller claim, and the point is that it is a CHECKABLE one.
     *
     * THE COMPARATOR WENT ON 2026-09-20 AND THE CLAUSE STAYED. It read "A bench
     * assay on that section could refute it; a ranking never could." That was
     * the only occurrence of the word "ranking" in rendered copy on this page,
     * and results/claim_surface_audit_2026-08-30.md section 2.4 recorded it as
     * a judgment call rather than settling it. It is not a judgment call
     * against the rule as written: the prohibition at the head of this file
     * covers a cross-compound ordering "in any form, hedged or not, POSITIVE OR
     * NEGATIVE — including as a limit", and a sentence whose rhetorical work is
     * done by contrast with an ordering is an ordering used as a limit. It also
     * does the specific damage the prohibition exists to prevent: a reader who
     * meets "a ranking never could" on a product page has been told the product
     * has one, and that the page is merely being modest about it.
     *
     * The requirement above survives the edit intact — the refutability clause
     * is still here and still attached to the finding. What it no longer does
     * is buy its force from the withdrawn feature. Four characters shorter,
     * measured, so the beat's registered height ceiling is not touched.
     *
     * WHAT THE VERDICT SENTENCE OWES, AND WHERE THE DEBT IS PAID. "The drug
     * arrives almost everywhere" is the model's own reading and is the exact
     * direction in which the engine is measured to be WRONG — its one
     * out-of-sample perivascular gradient is 2.2-3.0x too shallow on the map,
     * which flatters reach. (That figure was ~14x here until 2026-09-20; see
     * the row itself on beat 06 for why the larger number was wrong and why
     * correcting it DOWNWARD was still the honest edit.) "Under-dosed" also reads as "raise the dose", which the product
     * explicitly warns against because most of those cells are resistance-
     * flagged as well. Neither caveat fits in this beat (see the ceiling
     * above), so both are ledger rows on beat 06 and both name this sentence.
     * If this beat is ever given room, those two are what it should buy first.
     */
    body: "A compound compiles from its own structure: two routes through the cell layer, fed by lipophilicity and hydrogen-bond donors as well as mass, potency from measured GDSC2 screens. Structure reaches transport — scored unfitted on 68 measured compounds, but on descriptor columns supplied by that measurement, and the engine's own path has never been scored. It recovers about a quarter of the measured spread, so close analogues are under-separated, not resolved. Each run then decomposes survival cell by cell. On the flagship section 99.8% of surviving tumour cells are under-dosed, not unreached: the drug arrives almost everywhere and is almost everywhere too weak. A bench assay on that section could refute it; that is the point.",
    align: "right",
  },
  {
    kicker: "04",
    eyebrow: "Spatial attribution",
    headline: "The mean is not the finding. The spread is.",
    /*
     * TWO CORRECTIONS LANDED HERE ON 2026-09-20 AND BOTH RUN AGAINST US. This
     * beat had the most room on the page and was the last one audited, which
     * is the usual order in which a page goes wrong.
     *
     * (1) THE SECOND NULL IS NOW NAMED, BECAUSE WITHOUT IT THE BEAT DID NOT
     * RECONCILE WITH ITSELF. The sentence said one null "returns about 1.7" and
     * then that the attributable excess is 1.08-1.21x. 2.0 / 1.7 is 1.18, so a
     * reader who did the division could not reach 1.08 from anything on the
     * page — the low end comes from a SECOND primary null, the 400 µm block
     * shuffle at 1.8522, which was not mentioned. Quoting only the
     * clumpiness-matched null at 1.6581 also happens to be the flattering
     * choice: it is the one that yields the LARGER excess. Both nulls are
     * PRIMARY in results/immune_sink_impact.md sec. 5.3 and neither is a
     * sensitivity; omitting one was recorded as a defect in
     * results/CLAIM_INVENTORY.md sec. 3 — "quoting only the weaker null
     * flatters the result" — and it went uncorrected for five weeks.
     *
     * (2) THE TERM IS NOT IN THE PRODUCT, AND THE BEAT NOW SAYS SO. This is the
     * more serious of the two. `tmesim/molecular/immune_sink.py` is imported by
     * exactly one file — its own validation module — so no run this product
     * performs contains the mechanism these numbers measure. It is off BY
     * CONSTRUCTION rather than by default (immune_sink_impact.md: "no shipped
     * code path imports tmesim.molecular.immune_sink"), and the reason is
     * itself a published failure: it did not improve the measured human-glioma
     * anchor. A measurement that is sourced is not automatically a CAPABILITY,
     * and this beat sits on a product page under a panel labelled with a run —
     * every cue on the page said the product does this and it does not.
     *
     * WHAT THE BEAT STILL GETS TO SAY, because the correction is not a
     * withdrawal: immune_sink_impact.md sec. 6 endorses the finding in almost
     * these words — "It remains a genuine demonstration that geometry carries
     * information a compartmental model cannot express; it is not a 2.0x
     * demonstration." The demonstration survives; the size was already
     * corrected here; what was missing was that it is a demonstration rather
     * than a feature.
     */
    body: "In one section, tumor cells in the macrophage-densest quintile hold about 0.43 of drug exposure. In the sparsest, about 0.86 — a 2.0x spread inside a single tumor, where a well-mixed compartmental model returns exactly 1.0 by construction. That flatters us. Tissue-edge distance alone, with no immune input, reproduces 82% of it; nulls preserving the macrophages' clustering return 1.66 and 1.85. What these macrophages add is 1.08-1.21x — smaller, and measured. A compartment expresses neither. No shipped run carries the immune drug-loss term this needs: it did not improve the measured human-glioma anchor.",
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
         * 2026-08-22 — THIS ROW IS NO LONGER ABOUT RANKING, BECAUSE THE PRODUCT
         * NO LONGER RANKS. Everything below this block is the history of the
         * row as a RANKING row and is kept so the retracted claim cannot drift
         * back. What follows is why the row changed and what it may not lose.
         *
         * WHY THE RANKING FRAMING HAD TO GO RATHER THAN BE SOFTENED. The row
         * read "Ranking on one tissue · FULL ORDER NOT ROBUST". After 23e51f6
         * there is no ordering in the CLI, the API, the Studio, the protocol
         * language or the reports. A ledger row registering that a withdrawn
         * feature's output is not robust is not a disclosure — it is the
         * strongest possible hint that the feature is still shipping with a
         * caveat on it, which is precisely the reading the removal exists to
         * prevent. Retiring a caveat along with the claim it qualifies is not
         * the same edit as shrinking a caveat, and this is the only place on
         * the page where that distinction was applied.
         *
         * WHAT WAS RETIRED, EXACTLY, AND WHAT WAS KEPT. Retired: the ordering
         * breakpoints (36.3 on glioblastoma, 5.28 on breast, both below the
         * best constant of 39), the invariance list, and the delivery-law
         * top-pick clause — all four are statements ABOUT AN ORDERING and all
         * four are preserved verbatim in guardrail 1 at the head of this file,
         * which is now the record rather than an editing rule. KEPT, because it
         * is a live limit on every number the engine produces including the
         * map: the engine assumes a tissue capacity of 1.0 for every compound
         * while measured rat brain-slice values run 3.2 to 4184 (n = 33), and
         * that range is rat brain slices applied to human tissue, which is
         * itself an assumption. That provenance clause cost a line to buy on
         * 2026-08-21 and is not being handed back.
         *
         * WHAT THE ROW GAINED, AT NO COST IN LINES. The exact capacity
         * invariance and its exact exception. A converged gradient cannot move
         * with the capacity term at all — R divides diffusivity, velocity,
         * source and sink alike, so C_ss = S/k is fixed, verified to 1e-9
         * (results/capacity_ranking_robustness.md sec. 2b) — but the product
         * runs a SIX-HOUR window, which has not converged, and the falloff
         * shape does move inside it. That contrast is the disclosure: two
         * adjacent figures of different robustness cannot be fixed by
         * per-figure caveats, which is why the row states the invariance and
         * its failure in one breath rather than in two rows. It is stated
         * WITHOUT the measured movement (0.3772 absolute at R = 39, penetration
         * length collapsing 279.7 -> 76.7 um) because those figures live only
         * in app/server/reach.py and app/web/js/ui/reach.js and no results/
         * file carries them yet — see the citation guardrail at the head of
         * this file. The qualitative statement loses nothing.
         *
         * THE ROW IS 45 CHARACTERS SHORTER THAN THE ONE IT REPLACES. That is
         * not licence to spend the margin: at 390x844 beat 05's worst-case
         * bottom edge was measured at 832.9 of 844 before this pass, i.e. 11.1px
         * of headroom, which is under one ledger detail line. The four
         * disclosures this pass owed the page did not fit here and are on beat
         * 06 for that reason, each naming what it limits.
         *
         * ---- HISTORY BELOW THIS LINE ----
         *
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
          label: "Tissue capacity",
          state: "ASSUMED 1.0",
          detail:
            "Every compound is given a capacity of 1.0; measured rat brain-slice values run 3.2 to 4184 (n = 33), assumed to carry to human tissue. It cannot move a converged field at all — checked to one part in a billion — but the verdict in 03 is read off a finite exposure window, and there it does move.",
        },
        {
          label: "Close analogues",
          state: "0.498–0.839, NOT 0.718",
          detail:
            "0.718 was scored on the paper's donor counts, in a convention the model does not use. Redone on the engine's, 45 of 71 rows pin exactly and 26 have no published structure — until those are obtained, R² sits anywhere in that interval. It contains 0.718: nothing refuted, nothing confirmed, and a band that wide settles nothing. It still recovers only ~24% of the measured spread, from one cell line, under-separating a congeneric series ~4-fold. Its predecessor scored 0.022, below the whole interval; we published that refutation, then its replacement, and now this.",
        },
        /*
         * THE NULL COUNT WAS WRONG HERE UNTIL 2026-09-20, IN THE DIRECTION
         * THAT HELPS US. The row said "A null preserving the macrophages' own
         * clustering returns ~1.7" — one null, rounded. There are TWO primary
         * nulls in results/immune_sink_impact.md sec. 5.3, 1.6581 and 1.8522,
         * and the omitted one is the SECOND, which yields the SMALLER excess
         * (1.081x against 1.208x). So the singular and the rounding between
         * them concealed the low end of this row's own verdict column: nothing
         * on the page divided into 2.0 to give 1.08. Both now appear here and
         * in beat 04, at their measured precision, and the row came out 15
         * characters shorter for it.
         */
        {
          label: "Spatial spread",
          state: "1.08–1.21×, NOT 2.0×",
          detail:
            "A well-mixed model's 1.0 is not what bounds the raw 2.0×. Two spatially-aware nulls return 1.66 and 1.85, and tissue-edge distance alone — no immune input — reproduces 82% of it.",
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
    /*
     * "GET BACK A RANKING" WAS LEFT STANDING ON 2026-08-21 AND IS GONE AS OF
     * 2026-08-22, BECAUSE THE REASON IT WAS LEFT STANDING NO LONGER HOLDS. The
     * note that used to sit here argued the phrase "describes what the pipeline
     * emits, which is unchanged". As of 23e51f6 the pipeline does not emit it:
     * ranking is out of the CLI, the API, the Studio, the protocol language and
     * the reports. The sentence had become a promise of a feature that is not
     * there, which is a worse defect than the overclaim it was cleared of. See
     * the absolute prohibition at the head of this file.
     *
     * WHAT REPLACES IT IS THE PRODUCT'S ACTUAL OUTPUT, IN ITS ACTUAL ORDER: the
     * per-cell verdict first, and the exposure map second and explicitly AS THE
     * EVIDENCE FOR IT rather than as the thing being sold. That ordering is not
     * a stylistic choice — it is what app/web/js/ui/reach.js puts on the panel
     * and why, and a page that inverted it would sell the map and bury the
     * finding.
     *
     * "The same list, on your tissue, with your numbers in it" survives in
     * substance and is still the sentence that makes this beat honest: it is a
     * promise that the limits TRAVEL with the output. As of this pass it is
     * more than a promise here, because the limits are now literally below it.
     */
    body: "Vivum runs where your tissue already is: the solve is a single process on your own machine, and it sends no tissue or compound data anywhere. De-identified tissue only. Bring a section and a compound, and get back a per-cell verdict — delivery or potency — the map that is the evidence for it, and the limits below attached to both, on your tissue, with your numbers in them.",
    align: "center",
    /*
     * WHY A SECOND LEDGER EXISTS, AND WHY IT IS HERE RATHER THAN ON BEAT 05.
     * NOT because this is a good home for a limit. A limit belongs next to the
     * thing it limits, and beat 05 is the limits beat — this is the hand-off.
     * It is here because beat 05 is physically full and the alternative was to
     * ship the page softer than the product, which is the one thing that is not
     * allowed.
     *
     * THE MEASUREMENT, so that a future pass does not re-litigate it. Against a
     * production `next start` build at 390x844, 240 scroll positions, worst
     * case over every opacity: beat 05's bottom edge is 832.9px of 844. That is
     * 11.1px of headroom and a ledger detail line is 13.83px, so beat 05 could
     * not take ONE new row, let alone four. Nothing on beat 05 was compressed
     * to make room and nothing was bought out of an existing caveat. The one
     * row that changed there changed because its subject was withdrawn, and it
     * came out 45 characters SHORTER, which was left as margin.
     *
     * EVERY ROW BELOW NAMES WHAT IT LIMITS, in its own first clause. That is
     * the condition on parking a limit away from its claim: a reader must be
     * able to connect it back without scrolling up to find out what it is
     * about. Rows that could not satisfy that condition would not belong here.
     *
     * THE ORDER IS BY HOW BADLY A READER IS MISLED WITHOUT THE ROW, not by
     * topic. The delivery law is first because it is the only limit on this
     * page that can INVERT the verdict rather than bound it — a different
     * source term and the answer changes, not just its error bar. The second
     * and third rows are the two caveats beat 03's own sentence created and
     * could not afford to carry; they are here in preference to the barrier row
     * for that reason, because a caveat on a sentence this page states outranks
     * a general property of the map.
     *
     * WHAT DID NOT FIT, RECORDED AS A KNOWN OMISSION RATHER THAN AN OVERSIGHT.
     * The falloff SHAPE's measured movement under a capacity term is not a row
     * anywhere; it is folded into beat 05's "Tissue capacity" row as the
     * qualitative statement that the shape moves inside a six-hour window. It
     * is fully disclosed in the app directly beside the figure it bounds, which
     * is the right place for it. If this beat is ever given more room, the next
     * thing it should buy is the same "what the map cannot see" row split into
     * the two labelled quantities the commit message calls for — how much drug
     * crosses the wall, and where it goes once in — rather than one map and a
     * caveat.
     *
     * NO NUMBER BELOW IS QUOTED FROM A RUN. Every figure traces to a results/
     * file: 12.8x / 10.5x to results/breast_kp_validation.md; 7 of 12 and
     * 57.2 % to results/verdict_inversion_census.md sec. 0 and sec. 1; the
     * band-threshold reversal to results/mechanism_attribution_audit_2026-09-20
     * .md sec. 3.3 (0.2918082 against a shipped 0.1, and p50 = 0.2918052);
     * 2.2-3.0x, 2.45x and the direction of the error to
     * results/kyle2004_penetration.md sec. 6's CORRECTED block of 2026-08-22
     * (whose own words on the direction are "the model makes the drug look
     * like it penetrates better than it does"); and 34 of 38 to
     * results/synthesised_positions.md sec. 5. The barrier row and the shape
     * statement are deliberately NUMBERLESS — see the citation guardrail at the
     * head of this file.
     *
     * A TRAP THIS ROW SET AND SPRANG, RECORDED BECAUSE IT WILL RECUR: the 14x
     * was cited to sec. 0 item 3 of that file, and sec. 0 item 3 still says
     * 14x. The correction lives in sec. 6, in a quoted block, and the file's
     * own conclusions were never re-pointed at it. Citing a results/ file by
     * section is not the same as reading the file, and a number can be stale in
     * one paragraph of its own source while being correct in another.
     */
    ledger: {
      caption: "What the answer above rests on",
      rows: [
        /*
         * UPGRADED 2026-09-20 FROM "CAN INVERT IT" TO A MEASURED FACT, AND THE
         * DIRECTION OF THE EDIT IS THE POINT: the row got WORSE for us.
         *
         * "Can invert" was written when the evidence was one compound on one
         * section (results/permeability_split.md). results/verdict_inversion_
         * census.md (2026-08-25) ran the whole shipped panel and the answer is
         * not a possibility: the dominant limiter inverts for 7 of the 12
         * scored compounds, and the one the page's own lead is read from —
         * beat 03's section, beat 03's verdict — is one of them, with 57.2 %
         * of its tumour cells changing label. The census's own words: "The
         * single most-quoted run in this repository has a headline verdict
         * that is not robust to a term the model does not have."
         *
         * A hedge that a later measurement has settled is no longer a hedge.
         * It reads as caution and functions as concealment, because a reader
         * who is told a thing CAN happen assumes nobody has checked whether it
         * does. Somebody checked. The row says what they found.
         *
         * WHAT IS NOT CLAIMED, and the census is explicit about it (sec. 4):
         * the corrected arm is not the truth either. Making it honest needs a
         * per-compound wall permeability that swings entirely on an unmeasured
         * vessel radius. So the row says the verdict is not DETERMINED by the
         * data — not that the other answer is right.
         *
         * 12.8× / 10.5× stays, and now leads the row as the property of the law
         * rather than trailing it. It is the external anchor that says no single
         * constant repairs the law; what follows is what that costs the reader.
         *
         * A SECOND, INDEPENDENT ROUTE TO THE SAME REVERSAL WAS ADDED THE SAME
         * DAY, AND ITS SOURCE SAYS IT IS DISCLOSED NOWHERE.
         * results/mechanism_attribution_audit_2026-09-20.md sec. 3.3 sweeps
         * `unreached_ratio`, the constant that divides "never arrived" from
         * "arrived too weak". It ships at 0.1 and the file calls it, in its own
         * words, "an admittedly conventional constant". Beat 03's verdict
         * inverts DOSE -> DELIVERY at 0.2918082 — which is not an arbitrary
         * number either: it is the median surviving cell's own E/IC50 (p50 =
         * 0.2918052), i.e. the verdict says DOSE because slightly more than half
         * the survivors sit above a line we chose, and moving that line to the
         * middle of the distribution it is cutting moves the answer.
         *
         * THAT IS A WEAKER CLAIM THAN THE CENSUS'S AND IT IS WORTH MORE, WHICH
         * IS WHY IT EARNED THE LINE. The census needs a physics term the model
         * does not have; this needs nothing at all. Two independent routes to
         * one reversal is a different statement from one route — it is the
         * difference between "a missing term could overturn this" and "this
         * verdict is not determined by the tissue". The row's label and verdict
         * column changed to carry that: it is no longer a row about the
         * delivery law, it is a row about what the verdict rests on, and the
         * delivery law is one of the two things it rests on.
         *
         * NOT WRITTEN, because it is the product's disclosure rather than this
         * page's: the audit's actual complaint is that the Studio panel prints
         * the threshold's VALUE and never says how close the run is to it. That
         * belongs beside the number, in app/web/js/ui/. This page owes the
         * reader the fact, not the fix.
         */
        {
          label: "The 03 verdict",
          state: "INVERTS TWO WAYS",
          detail:
            "It comes from a one-way vascular leak with no path from tissue back to blood, which misses human microdialysis 12.8× high on one drug and 10.5× low on another, so no constant fixes both. Supply the missing term and the dominant limiter flips for 7 of 12 compounds, 03's included, relabelling 57.2% of its cells. Moving one conventional band threshold to the median survivor's own exposure flips it again.",
        },
        /*
         * THE NUMBER CHANGED ON 2026-09-20 AND THE ROW GOT SMALLER, WHICH IS
         * THE HARDER CORRECTION TO MAKE AND THE ONE MOST LIKELY TO BE SKIPPED.
         *
         * It read "14× TOO SHALLOW". That figure was faithfully copied from
         * results/kyle2004_penetration.md, and the same file corrected it on
         * 2026-08-22 in a block the copy never picked up: the 14× is a RATIO OF
         * TWO DIMENSIONLESS RATIOS and Kyle's near:far readout is EXPONENTIAL
         * in the transport length, so 14× is a 2.45× error in λ (2.00–2.98×
         * with the paper's own SDs) amplified by dlnR/dlnλ ≈ 2.9. On this
         * product's own readout the amplification does not transfer — the map's
         * half-distance responds to λ with exponent 0.92–1.22, near-linear
         * against 2.92 — so the figure on the MAP is 2.2–3.0×, not 14×. This
         * row is attached to a sentence about the map, so the map's figure is
         * the one it owes.
         *
         * WHY IT WAS NOT LEFT ALONE, since 14× is the more self-critical
         * number and nobody would have complained. Because it is not true, and
         * a page that overstates its own error has given up the only thing that
         * makes its understatements worth reading. It also inflates a defect
         * the same file shows is largely INHERITED rather than independent —
         * 3.6696 × 3.7374 = 13.71, two thirds of it the residual on the
         * compound α is fit to — and the two anchors that bound it are
         * disjoint (Primeau 57.71–72.13 µm, Kyle 30.73–55.68 µm), so "fit α
         * better" is provably not the repair.
         *
         * NOT SAID HERE, and it is the sharper caveat: kyle2004_penetration.md
         * also records that the anchor may not measure a LENGTH at all — zero-
         * loss retarded diffusion reproduces the measured value at λ = ∞ — in
         * which case inverting it to λ is invalid. It does not fit in this row
         * and it would weaken a limit rather than a claim, so it is recorded
         * here and left off the page rather than quietly dropped.
         */
        {
          label: "“Arrives almost everywhere”",
          state: "2.2–3.0× TOO SHALLOW",
          detail:
            "The map's one out-of-sample perivascular gradient is that much too shallow, from a transport length measured 2.45× too long. The error flatters the drug: it under-reports how much tumour goes untreated.",
        },
        {
          label: "“Under-dosed”",
          state: "NOT “RAISE THE DOSE”",
          detail:
            "Most of those cells are flagged resistant too, so more drug alone would not have killed them. And no delivery verdict has been reached non-degenerately: the only one came from a compound whose barrier term is zero — an input, not a discovery.",
        },
      ],
    },
  },
] as const;

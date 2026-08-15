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
 * Two guardrails narrowed by measurement, both of which read as weaker than
 * the obvious marketing line and must not drift back to it:
 *
 * 1. RANKING IS NOT UNCONDITIONALLY INVARIANT. It is invariant to
 *    common-mode perturbations (cell-state route, oxygen route, alpha,
 *    grid) and NOT to drug-specific transport parameters — correcting the
 *    transvascular delivery law changes the top-1 GBM pick. Say the narrow
 *    version; the broad one was tested and is false.
 *    See results/transvascular_exchange_revalidation.md.
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
 *    consumes logP and hydrogen-bond donor count. Blind and unfitted on those
 *    same 71 compounds it reaches R^2 = 0.718, 99 % within 3x, and orders 24 of
 *    31 exact-mass pairs correctly against the old model's 0 of 31. So
 *    "a compound's structure drives its transport" is now SAYABLE.
 *    THE LIMIT THAT REPLACES THE OLD ONE, AND IT IS NOT OPTIONAL: the model
 *    reproduces 8.3x of a measured 29.7x spread -- about a QUARTER -- so it
 *    UNDER-SEPARATES similar-mass analogues rather than resolving them, and the
 *    whole test is ONE cell line. Do NOT write "validated", "accurate", "gets
 *    transport right", or any absolute diffusivity figure. Say that structure
 *    reaches transport, say it was scored blind against 71 measured compounds,
 *    and say it recovers about a quarter of the measured spread.
 *    See results/two_route_promotion.md, results/diffusivity_falsification.md.
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
    body: "Convection, diffusion and reaction, solved on that geometry. Interstitial pressure builds where vessels leak and drainage fails and pushes back against the drug coming in; oxygen and hypoxia solve on the same grid, and response follows Hill kinetics from the concentration each cell sees. What that buys is the shape of the field, not its magnitude — and the difference is measured, not asserted.",
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
     */
    body: "A compound compiles from its own structure: two routes through the cell layer, paracellular and transcellular, fed by lipophilicity and hydrogen-bond donors as well as mass, with potency from measured GDSC2 screens. Structure reaches transport now — scored blind against 71 measured compounds — but it recovers only about a quarter of their measured spread, so close analogues are under-separated rather than resolved. Each run then decomposes survival cell by cell: a molecule that never arrives and one that arrives too weak look identical in a survival curve, and need opposite fixes.",
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
     */
    body: "Ranking compounds on the same tissue is the output to trust — far more robust than any absolute number, and still not unconditionally invariant. Here is where it breaks.",
    align: "right",
    ledger: {
      caption: "Carried by every run · none of this is buried",
      rows: [
        {
          label: "Ranking on one tissue",
          state: "CONDITIONAL",
          detail:
            "Invariant to cell-state route, oxygen route, alpha and grid — not to drug-specific transport. Correcting the delivery law changes the top pick on glioblastoma.",
        },
        {
          label: "Close analogues",
          state: "UNDER-SEPARATED",
          detail:
            "Two-route transport (logP + hydrogen-bond donors) scores R² = 0.718 blind on 71 measured compounds, but recovers only ~24% of their spread, from one cell line — under-separating a congeneric series ~4-fold. Its predecessor scored 0.022; we published that refutation, then its replacement.",
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
    body: "Codon runs where your tissue already is: the solve is a single process on your own machine, and it sends no tissue or compound data anywhere. De-identified tissue only. Bring a section and a candidate series, and get back a ranking, a per-cell mechanism attribution of what killed each compound, and the list above attached to both — the same list, on your tissue, with your numbers in it.",
    align: "center",
  },
] as const;

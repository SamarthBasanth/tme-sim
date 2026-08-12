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
 * Two guardrails narrowed by measurement, both of which read as weaker than
 * the obvious marketing line and must not drift back to it:
 *
 * 1. RANKING IS NOT UNCONDITIONALLY INVARIANT. It is invariant to
 *    common-mode perturbations (cell-state route, oxygen route, alpha,
 *    grid) and NOT to drug-specific transport parameters — correcting the
 *    transvascular delivery law changes the top-1 GBM pick. Say the narrow
 *    version; the broad one was tested and is false.
 *    See results/transvascular_exchange_revalidation.md.
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

export interface SectionCopy {
  kicker: string;
  eyebrow: string;
  headline: string;
  body: string;
  align: "left" | "right" | "center";
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
    headline: "Then run physics on it.",
    body: "Convection, diffusion and reaction, solved on that geometry. Interstitial fluid pressure builds where vessels leak and drainage fails, and pushes back against the drug coming in. Oxygen and hypoxia solve on the same grid. Response follows Hill kinetics from the concentration each cell actually sees.",
    align: "left",
  },
  {
    kicker: "03",
    eyebrow: "Delivery or potency",
    headline: "Which one killed the compound?",
    body: "A compound compiles from its own structure: descriptors to Stokes-Einstein diffusivity to hindered transport through crowded ECM, with potency from measured GDSC2 screens. Each run then decomposes survival cell by cell — not how many lived, but why each one did. A molecule that never arrives and a molecule that arrives too weak look identical in a survival curve and need opposite fixes.",
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
    headline: "Every run ships its own caveats.",
    body: "Ranking compounds on the same tissue is far more robust than any absolute number, and it is the output to trust. It is not unconditionally invariant, and we can tell you exactly where it breaks: correcting the transvascular delivery law re-weights the panel by permeability and changes the top pick on glioblastoma. That is written up, not buried. The spatial spread needs single-cell resolved tissue — on spot data the immune positions are synthesized rather than measured, so we do not report it there. Absolute kill fraction is not defensible at all, and every run says so in its own caveats panel. Vivum does not predict patient response, and its null results are published.",
    align: "right",
  },
  {
    kicker: "06",
    eyebrow: "Get started",
    headline: "Run it on your own tissue.",
    body: "Vivum runs inside your own infrastructure. Proprietary tissue and compound data never cross your firewall. Bring a section and a candidate series, and get back a ranking, a per-cell mechanism attribution, and the caveats that come with them.",
    align: "center",
  },
] as const;

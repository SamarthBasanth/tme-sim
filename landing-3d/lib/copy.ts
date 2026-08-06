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
 * is not a defensible output; relative ranking on the same tissue and
 * spatial delivery attribution are. Every number that appears below must
 * be traceable to a real run or a real dataset count.
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
    body: "The section behind this text is real. 40,887 cells from a human glioblastoma, each one where the assay found it, colored by measured cell type. Vivum rebuilds tumor architecture from spatial transcriptomics this way: cells, vasculature, stroma, ECM and oxygen, across 38 sections from 13 patients. Nothing about the shape was invented, including the tears. It reads flat because it is: one 8.1 by 5.1 mm section, with nothing measured along the third axis.",
    align: "center",
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
    body: "A compound compiles from its own structure: descriptors to Stokes-Einstein diffusivity to hindered transport through crowded ECM, with potency from measured GDSC2 screens. Each run then decomposes survival cell by cell. One returns delivery 0.172, dose 0.825, resistance 0.003. A molecule that never arrives and a molecule that arrives too weak need different fixes.",
    align: "right",
  },
  {
    kicker: "04",
    eyebrow: "Spatial attribution",
    headline: "The mean is not the finding. The spread is.",
    body: "In one section, tumor cells in the macrophage-densest quintile hold about 0.43 of drug exposure. In the sparsest, about 0.86. That is a 2.0x spread inside a single tumor, and it falls monotonically across every quintile between them. A well-mixed compartmental model returns exactly 1.0 there, by construction. Geometry produces quantities a compartment cannot express.",
    align: "left",
  },
  {
    kicker: "05",
    eyebrow: "Stated limits",
    headline: "Every run ships its own caveats.",
    body: "Ranking compounds on the same tissue is invariant and defensible. Absolute kill fraction is not, and every run says so: a caveats panel states what each number can and cannot be read as. Validation is against measured tissue penetration depths and tissue:plasma ratios, not clinical outcomes. Vivum does not predict patient response, and its null results are published.",
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

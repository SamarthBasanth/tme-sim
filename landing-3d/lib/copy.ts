/**
 * Narrative copy, shared between the cinematic scroll experience and the
 * prefers-reduced-motion fallback so the two never drift out of sync.
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
    eyebrow: "3D tumor visualization",
    headline: "See the Tumor Before You Treat It.",
    body: "A real-time, three-dimensional view of tumor anatomy — built to bring surgical and treatment planning closer to what's actually there.",
    align: "center",
  },
  {
    kicker: "02",
    eyebrow: "Anatomical context",
    headline: "Every mass, in its place.",
    body: "The camera pulls back to reveal the tumor inside its surrounding anatomy — never viewed in isolation from the structure around it.",
    align: "left",
  },
  {
    kicker: "03",
    eyebrow: "Segmentation",
    headline: "Core. Edema. Margin.",
    body: "The mass separates into its constituent regions, each carrying its own material signature — so what's dense, swollen, or boundary reads at a glance.",
    align: "right",
  },
  {
    kicker: "04",
    eyebrow: "AI analysis",
    headline: "A closer pass, region by region.",
    body: "A scanning plane sweeps the volume, surfacing regions worth a second look. The interface stays restrained on purpose — this is a reading aid, not a verdict.",
    align: "left",
  },
  {
    kicker: "05",
    eyebrow: "Treatment planning",
    headline: "Plan with millimeter-level context.",
    body: "A proposed trajectory is traced against the model in real time, with nearby critical structures brought into view before a plan is finalized.",
    align: "right",
  },
  {
    kicker: "06",
    eyebrow: "Get started",
    headline: "Bring this into your own workflow.",
    body: "The full model reassembles, resting alongside the case for one last look before you go.",
    align: "center",
  },
] as const;

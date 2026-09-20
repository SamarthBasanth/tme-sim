import type { Metadata } from "next";
import {
  Archivo,
  Atkinson_Hyperlegible_Next,
  Atkinson_Hyperlegible_Mono,
} from "next/font/google";
import "./globals.css";

/*
 * WHY THESE FACES.
 *
 * TWO PREVIOUS STACKS ARE NAMED HERE BECAUSE EACH FAILED A DIFFERENT WAY, AND
 * NEITHER FAILURE MODE MAY COME BACK.
 *
 * (1) A system stack (Helvetica Neue / SF Mono) with a webfont declared but
 *     never actually loaded. Every heading silently resolved to Helvetica Neue
 *     Semibold with negative tracking on a near-black ground. That is the
 *     single most generic dark-product look there is, and nobody noticed for
 *     the life of the page because a font that fails to load does not error --
 *     it just renders. This is why the three faces below were confirmed to
 *     FETCH (HTTP 200 from the Google Fonts CSS endpoint) before being chosen,
 *     and why the rendered page was screenshotted and read rather than assumed.
 *
 * (2) IBM Plex Sans / Condensed / Mono. Correct on subject lineage -- Plex was
 *     commissioned for technical and scientific documentation -- but it is the
 *     family Anthropic's Claude interface is set in, so this page read as a
 *     skin of a chat product rather than as an instrument. Do not return to it.
 *
 * The network objection that originally pinned this to system fonts does not
 * apply to `next/font/google`: it SELF-HOSTS, fetching the woff2 at build time
 * and serving it from this origin. No runtime request to Google, nothing for
 * the strict CSP to block, no privacy leak.
 *
 * ---------------------------------------------------------------------------
 * THE ARGUMENT FOR THIS PARTICULAR TRIAD, IN TERMS OF THE SUBJECT.
 *
 * This page's thesis is not that the simulator is good. It is that the
 * simulator's LIMITS are published, legible, and checkable -- beat 05 is a
 * spec sheet of the product's own failures. So the typographic problem is not
 * "look scientific". It is: a reader must be able to take a number away
 * without misreading it, and must actually finish a caveat. lib/copy.ts states
 * this outright -- "a limit a reader does not finish is not a disclosure".
 * That converts legibility from taste into a product requirement, and there is
 * a family drawn against exactly that requirement.
 *
 *   PROSE + MONO -- Atkinson Hyperlegible Next / Atkinson Hyperlegible Mono.
 *   Commissioned by the Braille Institute and drawn for low-vision readers,
 *   its explicit design goal is CHARACTER DISAMBIGUATION: 0 vs O, 1 vs l vs I,
 *   6 vs 8, 5 vs S are separated by structural difference -- a slashed zero, a
 *   tailed l, an open-tailed 1 -- rather than by a hairline that dies at 10px
 *   on a dark ground. Every quantity on this page is meant to be transcribed
 *   and checked against a results/ file. A face engineered so a digit cannot
 *   be confused for a letter is an instrument-grade choice, not a stylistic
 *   one, and it is the rare case where the type specimen and the product
 *   argument are the same sentence.
 *
 *   Taking BOTH text roles from one 2025-matched superfamily is also load
 *   bearing rather than tidy. The ledger sets a sans label and a mono verdict
 *   on ONE shared baseline (.ledger__row), and both instrument panels pair a
 *   sans <dt> against a mono <dd> the same way. A mono drawn against its own
 *   sans keeps x-height and baseline in agreement, so those pairs sit level
 *   without per-element nudging. A borrowed mono would not.
 *
 *   DISPLAY -- Archivo, variable, with the width axis loaded.
 *   Omnibus-Type's revival of late-19th-century American grotesques, drawn for
 *   HIGHLIGHTING text in newspapers and magazines: headline lettering whose job
 *   is to report a finding, not to sell one. That is the register these
 *   headlines are in -- "The shape is solved. The scale is not." is a result,
 *   not a promise.
 *
 *   The `wdth` axis is requested for a structural reason, not for choice. The
 *   hierarchy on this page is carried by WIDTH CONTRAST against the body face,
 *   not by size alone, and .scroll-pin clips at overflow:hidden -- so a display
 *   face that sets one line wider than its predecessor does not degrade, it
 *   silently truncates a disclosure. A continuous width axis lets the display
 *   be tuned to the column instead of hoping a fixed-width cut fits.
 *
 * WHAT WAS DELIBERATELY NOT CHOSEN: IBM Plex (above), and Inter / Space
 * Grotesk / Geist, which are the current house style of generated design.
 * Picking one of those would have traded one generic look for another.
 *
 * All three are variable fonts, so `weight: "variable"` ships one file per
 * family across the whole weight range rather than a static instance per step.
 *
 * ---------------------------------------------------------------------------
 * FALLBACK METRICS, AND WHY TWO OF THESE THREE CARRY THEM BY HAND.
 *
 * `next/font/google` normally emits a second, metric-matched @font-face -- a
 * local system face squeezed by `size-adjust` with `ascent-override` /
 * `descent-override` on top -- so the fallback shown before the woff2 lands
 * occupies the SAME BOX as the real face and the swap moves nothing. It reads
 * those numbers out of a precalculated table
 * (next/dist/server/capsize-font-metrics.json). That table has `archivo` and
 * the ORIGINAL `atkinsonHyperlegible`, and it does not have the two 2025
 * faces this page actually uses, so the build printed:
 *
 *   Failed to find font override values for font `Atkinson Hyperlegible Mono`
 *   Failed to find font override values for font `Atkinson Hyperlegible Next`
 *   -- Skipping generating a fallback font.
 *
 * With `display: "swap"` and no override metrics that is a real reflow on
 * first paint, and it lands on the one thing this page's typography exists to
 * guarantee: .ledger__row sets a sans label and a mono verdict on one shared
 * baseline, and both instrument panels pair a sans <dt> against a mono <dd>
 * the same way. Measured before the fix, on the fallback face, the whole
 * reduced-motion page reflowed by up to 52px vertically and the ledger's
 * mono verdicts changed width by 5.6px.
 *
 * So the overrides are supplied by hand, in app/globals.css, from metrics
 * MEASURED off the loaded faces rather than guessed -- see the block there for
 * the numbers, the method and the platform caveat. `adjustFontFallback: false`
 * is what stops Next from trying and failing; it is not a way of silencing the
 * warning. If those @font-face rules are ever deleted, turn this back on (it
 * still will not work for these families) and re-measure.
 */
/*
 * Archivo is the one of the three Next DOES have metrics for, and its
 * generated fallback was still wrong -- because the table stores one width and
 * this page never uses the default width. Headings are set at `font-stretch:
 * 80%`, so the fallback Next computed was ~18% too wide and the headings
 * visibly reflowed on swap. The replacement is four stretch-partitioned faces
 * in globals.css; this flag is what stops Next generating a fifth that would
 * compete with them.
 */
const archivo = Archivo({
  subsets: ["latin"],
  weight: "variable",
  axes: ["wdth"],
  variable: "--font-archivo",
  display: "swap",
  adjustFontFallback: false,
  fallback: ["Archivo Width Fallback"],
});

const atkinson = Atkinson_Hyperlegible_Next({
  subsets: ["latin"],
  weight: "variable",
  style: ["normal", "italic"],
  variable: "--font-atkinson",
  display: "swap",
  adjustFontFallback: false,
  fallback: ["Atkinson Hyperlegible Next Fallback"],
});

const atkinsonMono = Atkinson_Hyperlegible_Mono({
  subsets: ["latin"],
  weight: "variable",
  variable: "--font-atkinson-mono",
  display: "swap",
  adjustFontFallback: false,
  fallback: ["Atkinson Hyperlegible Mono Fallback"],
});

/*
 * "validated against real measured tissue benchmarks" was the last surviving
 * unqualified "validated" on any marketing surface, and it sat on the highest-
 * visibility one there is: this is the string search engines and link previews
 * render. It violated guardrail 1b in lib/copy.ts -- `Do NOT write "validated"`
 * -- in the same repository that forbids it, and the README records that
 * `vivum validate` reports more failures than passes.
 *
 * "Scored, with the failures published alongside the passes" is not a softened
 * version of the same boast. It is the actual differentiator: anyone can claim
 * a benchmark, and publishing the ones you lose is the part competitors do not
 * copy. Keep the failures clause attached to the scoring clause -- split apart,
 * the first half reads as the claim this replaced.
 *
 * REMOVED 2026-09-20: "Built for pharma, biotech, and research teams COMPARING
 * CANDIDATE THERAPIES." That is cross-compound comparison framing, and it was
 * sitting on the one surface a reader meets before the page loads -- a search
 * result, a link unfurl, a shared card. lib/copy.ts:146-149 admits no version
 * of it: "No ordering, ranked list, shortlist, 'top pick', 'which compound
 * wins', pairwise call or cross-compound comparison may appear on this page in
 * any form, hedged or not, positive or negative." The product removed the
 * FEATURE on 2026-08-22 (23e51f6) and this string outlived it by a month.
 * Flagged in results/claim_surface_audit_2026-08-30.md section 2.4 as
 * PROHIBITED, with the mechanical cause named: this file does not import
 * copy.ts, so nothing couples the guardrail to the string it governs. That is
 * still true and is recorded in the report rather than papered over -- the
 * coupling that would fix it is a build-time check, not a comment.
 *
 * What replaces it says what the product is FOR without ordering anything, and
 * both halves are already rendered on the page itself (beat 05's last two
 * ledger rows), so the snippet cannot drift from the body.
 */
export const metadata: Metadata = {
  title: "Vivum — See Where the Drug Actually Goes",
  description:
    "Vivum reconstructs real tumor tissue from spatial biology data and simulates drug transport physics through it — penetration, exposure, and delivery, scored against real measured tissue benchmarks with the failures published alongside the passes. For research use: not for diagnostic or treatment-planning use, and patient response is not predicted.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html
      lang="en"
      className={`h-full antialiased ${atkinson.variable} ${archivo.variable} ${atkinsonMono.variable}`}
    >
      <body className="min-h-full">{children}</body>
    </html>
  );
}

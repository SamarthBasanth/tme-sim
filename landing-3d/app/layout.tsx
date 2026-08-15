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
 */
const archivo = Archivo({
  subsets: ["latin"],
  weight: "variable",
  axes: ["wdth"],
  variable: "--font-archivo",
  display: "swap",
});

const atkinson = Atkinson_Hyperlegible_Next({
  subsets: ["latin"],
  weight: "variable",
  style: ["normal", "italic"],
  variable: "--font-atkinson",
  display: "swap",
});

const atkinsonMono = Atkinson_Hyperlegible_Mono({
  subsets: ["latin"],
  weight: "variable",
  variable: "--font-atkinson-mono",
  display: "swap",
});

/*
 * "validated against real measured tissue benchmarks" was the last surviving
 * unqualified "validated" on any marketing surface, and it sat on the highest-
 * visibility one there is: this is the string search engines and link previews
 * render. It violated guardrail 1b in lib/copy.ts -- `Do NOT write "validated"`
 * -- in the same repository that forbids it, and the README records that
 * `codon validate` reports more failures than passes.
 *
 * "Scored, with the failures published alongside the passes" is not a softened
 * version of the same boast. It is the actual differentiator: anyone can claim
 * a benchmark, and publishing the ones you lose is the part competitors do not
 * copy. Keep the failures clause attached to the scoring clause -- split apart,
 * the first half reads as the claim this replaced.
 */
export const metadata: Metadata = {
  title: "Codon — See Where the Drug Actually Goes",
  description:
    "Codon reconstructs real tumor tissue from spatial biology data and simulates drug transport physics through it — penetration, exposure, and delivery, scored against real measured tissue benchmarks with the failures published alongside the passes. Built for pharma, biotech, and research teams comparing candidate therapies.",
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

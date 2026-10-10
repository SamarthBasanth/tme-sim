import type { Metadata } from "next";
import localFont from "next/font/local";
import "./globals.css";

/*
 * WHY THESE FACES, AND WHY THE FILES ARE IN THE REPOSITORY.
 *
 * THREE PREVIOUS STACKS ARE NAMED HERE BECAUSE EACH FAILED A DIFFERENT WAY,
 * AND NO FAILURE MODE MAY COME BACK.
 *
 * (1) A system stack (Helvetica Neue / SF Mono) with a webfont declared but
 *     never actually loaded. Every heading silently resolved to Helvetica Neue
 *     Semibold with negative tracking on a near-black ground. That is the
 *     single most generic dark-product look there is, and nobody noticed for
 *     the life of the page because a font that fails to load does not error --
 *     it just renders. This is why every face below was confirmed to PARSE
 *     (axis ranges and vertical metrics read out of the binary, see
 *     app/fonts/README.md) and why the rendered page was screenshotted and
 *     read rather than assumed.
 *
 * (2) IBM Plex Sans / Condensed / Mono. Correct on subject lineage -- Plex was
 *     commissioned for technical and scientific documentation -- but it is the
 *     family Anthropic's Claude interface is set in, so this page read as a
 *     skin of a chat product rather than as an instrument. Do not return to it.
 *
 * (3) ARCHIVO AS THE DISPLAY FACE, which is the one this pass replaced, and
 *     the reason is worth stating precisely because the choice was argued and
 *     the argument was sound. Archivo is a grotesque drawn for highlighting
 *     text in newspapers; the register was right. What was wrong was the
 *     RESULT: set at `font-stretch: 80%` and weight 650, white, on #07080a, it
 *     renders as a heavy condensed grotesque headline -- which is the single
 *     most common headline treatment on a generated dark product page, and is
 *     a near relative of exactly the failure (1) is about. The page had three
 *     typographic roles and all three were grotesque or mono, so the entire
 *     hierarchy was carried by size, weight and width: three settings of one
 *     voice. Nothing distinguished a FINDING from a CONTROL except how big it
 *     was.
 *
 * ---------------------------------------------------------------------------
 * THE ARGUMENT FOR THIS TRIAD, IN TERMS OF THE SUBJECT.
 *
 * This page's thesis is not that the simulator is good. It is that the
 * simulator's LIMITS are published, legible, and checkable -- beat 05 is a
 * spec sheet of the product's own failures. So the typographic problem is not
 * "look scientific". It is: a reader must be able to take a number away
 * without misreading it, and must actually finish a caveat. lib/copy.ts states
 * this outright -- "a limit a reader does not finish is not a disclosure".
 * That converts legibility from taste into a product requirement.
 *
 *   PROSE + MONO -- Atkinson Hyperlegible Next / Atkinson Hyperlegible Mono.
 *   UNCHANGED BY THIS PASS, AND DELIBERATELY SO. Commissioned by the Braille
 *   Institute and drawn for low-vision readers, their explicit design goal is
 *   CHARACTER DISAMBIGUATION: 0 vs O, 1 vs l vs I, 6 vs 8, 5 vs S are
 *   separated by structural difference -- a slashed zero, a tailed l, an
 *   open-tailed 1 -- rather than by a hairline that dies at 10px on a dark
 *   ground. Every quantity on this page is meant to be transcribed and checked
 *   against a results/ file. A face engineered so a digit cannot be confused
 *   for a letter is an instrument-grade choice, not a stylistic one. That is a
 *   product argument and it did not stop being true because the display face
 *   was wrong, so it was not traded for a change of mood.
 *
 *   Taking BOTH text roles from one 2025-matched superfamily is also load
 *   bearing rather than tidy. The ledger sets a label and a mono verdict on
 *   ONE shared baseline (.ledger__row), and both instrument panels pair a sans
 *   <dt> against a mono <dd> the same way. A mono drawn against its own sans
 *   keeps x-height and baseline in agreement, so those pairs sit level without
 *   per-element nudging. A borrowed mono would not.
 *
 *   DISPLAY -- Newsreader, variable, with the OPTICAL SIZE axis.
 *   Production Type's screen serif for long-form news. Three reasons, in
 *   descending order of how much they are about this page rather than about
 *   taste:
 *
 *   (a) IT INTRODUCES A SECOND VOICE, which is the actual defect in (3). A
 *       serif heading over a sans paragraph over a mono verdict is a hierarchy
 *       of KIND, not of degree: a reader can tell at a glance which register
 *       they are in without having measured anything. The page's three
 *       information types -- the claim, the prose, the reading -- now have
 *       three different shapes rather than three sizes of one.
 *
 *   (b) THE REGISTER IS A REPORT, NOT A PITCH. "The shape is solved. The scale
 *       is not." is a finding. A text serif at display size is how a journal
 *       or a broadsheet sets a finding; a heavy condensed grotesque is how a
 *       product sets a promise. The copy has been pulled toward the first for
 *       two years and the type was pulling the other way.
 *
 *   (c) THE `opsz` AXIS IS A REAL REPLACEMENT FOR ARCHIVO'S `wdth`, NOT AN
 *       ABSENCE OF ONE. The width axis was requested for a structural reason:
 *       .scroll-pin clips at overflow:hidden, so a display face that cannot be
 *       tuned to the column does not degrade, it silently truncates a
 *       disclosure. `opsz` is a continuous axis serving the same need from the
 *       other end -- at 60 the drawing tightens, sharpens and sheds
 *       sidebearing; at 14.4 (the ledger label, the wordmark) it opens up and
 *       thickens its hairlines so a 15px serif on a dark ground does not
 *       vanish. globals.css pins it per role rather than leaving it on the
 *       browser default of `auto`, and the @font-face block there records the
 *       measurement that forced that: on this face `auto` swings the advance
 *       7.6% across the headings' own responsive clamp, which is about a line
 *       on a phone h1.
 *
 * WHAT WAS DELIBERATELY NOT CHOSEN: IBM Plex and Archivo (above), and Inter /
 * Space Grotesk / Geist / Poppins / Montserrat, which are the current house
 * style of generated design. Picking one of those would have traded one
 * generic look for another.
 *
 * All four files are variable, so one file per family covers the whole weight
 * range rather than a static instance per step.
 *
 * ---------------------------------------------------------------------------
 * WHY `next/font/local` AND NOT `next/font/google`.
 *
 * The old comment here defended `next/font/google` against the network
 * objection, correctly: it SELF-HOSTS, fetching the woff2 at build time and
 * serving it from this origin, so there is no runtime request to Google and
 * scripts/network-audit.mjs passed with it. That defence is still true and it
 * answers the wrong question. The requests it does not make are the runtime
 * ones. The build still had to reach fonts.gstatic.com, which means:
 *
 *   - the build did not work without egress, so "this page contacts nothing"
 *     was a property of the page and not of the thing that produces it;
 *   - the bytes were not pinned, and Google re-cuts families in place, so two
 *     builds a month apart could ship different outlines and a different line
 *     count with nothing in the diff -- on a page where three beats have under
 *     12px of vertical headroom;
 *   - the OFL licence, which requires the licence to travel with the Font
 *     Software, was never distributed with the binary.
 *
 * All four .woff2 files are now committed in app/fonts/ with their licences.
 * 221 KB in the repository buys an offline, byte-pinned, correctly-licensed
 * build. See app/fonts/README.md for provenance and how to re-obtain them.
 *
 * ---------------------------------------------------------------------------
 * FALLBACK METRICS, AND WHY ALL THREE CARRY THEM BY HAND.
 *
 * `next/font` normally emits a second, metric-matched @font-face -- a local
 * system face squeezed by `size-adjust` with `ascent-override` /
 * `descent-override` on top -- so the fallback shown before the woff2 lands
 * occupies the SAME BOX as the real face and the swap moves nothing. For
 * `next/font/local` it derives those numbers from Arial or Times New Roman
 * only, and it gets them wrong here for a reason that survives the switch to
 * local files: a size-adjust is a CONSTANT, and two of these three faces are
 * set at a variation this page chooses rather than at the family default.
 *
 * With `display: "swap"` and no override metrics that is a real reflow on
 * first paint, and it lands on the one thing this page's typography exists to
 * guarantee: .ledger__row sets a label and a mono verdict on one shared
 * baseline, and both instrument panels pair a <dt> against a mono <dd> the
 * same way. Measured before the fix, on the fallback face, the whole
 * reduced-motion page reflowed by up to 52px vertically and the ledger's mono
 * verdicts changed width by 5.6px.
 *
 * So the overrides are supplied by hand, in app/globals.css, from metrics
 * MEASURED off the loaded faces rather than guessed -- see the block there for
 * the numbers, the method and the platform caveat. `adjustFontFallback: false`
 * is what stops Next emitting a competing one. If those @font-face rules are
 * ever deleted, turn this back on and re-measure.
 */

/*
 * `weight` is the full axis range, not a point, so one file serves every
 * weight the stylesheet asks for. The `opsz` axis needs nothing declared
 * here: there is no @font-face descriptor for optical size, and globals.css
 * drives the axis directly with `font-variation-settings` at each of the five
 * display roles.
 */
const newsreader = localFont({
  src: "./fonts/Newsreader-roman-latin.woff2",
  weight: "200 800",
  style: "normal",
  variable: "--font-newsreader",
  display: "swap",
  adjustFontFallback: false,
  fallback: ["Newsreader Fallback"],
});

const atkinson = localFont({
  src: [
    { path: "./fonts/AtkinsonHyperlegibleNext-roman-latin.woff2", style: "normal" },
    { path: "./fonts/AtkinsonHyperlegibleNext-italic-latin.woff2", style: "italic" },
  ],
  weight: "200 800",
  variable: "--font-atkinson",
  display: "swap",
  adjustFontFallback: false,
  fallback: ["Atkinson Hyperlegible Next Fallback"],
});

/*
 * The mono file's DEFAULT instance is ExtraLight (wght default 200), not
 * Regular -- read out of its fvar, see app/fonts/README.md. Declaring the
 * range here is therefore not cosmetic: without `weight: "200 800"` the
 * browser treats the file as a 200-weight static face, every mono verdict on
 * beat 05 renders hairline on a near-black ground, and the one role on this
 * page that carries measured quantities becomes the least legible thing on it.
 */
const atkinsonMono = localFont({
  src: "./fonts/AtkinsonHyperlegibleMono-roman-latin.woff2",
  weight: "200 800",
  style: "normal",
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
      className={`h-full antialiased ${atkinson.variable} ${newsreader.variable} ${atkinsonMono.variable}`}
    >
      <body className="min-h-full">{children}</body>
    </html>
  );
}

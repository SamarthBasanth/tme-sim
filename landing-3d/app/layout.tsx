import type { Metadata } from "next";
import {
  IBM_Plex_Sans,
  IBM_Plex_Sans_Condensed,
  IBM_Plex_Mono,
} from "next/font/google";
import "./globals.css";

/*
 * WHY THESE FACES, AND WHY THIS REPLACED THE SYSTEM STACK.
 *
 * The previous stack was Helvetica Neue / Neue Haas Grotesk with a SF Mono
 * fallback, chosen deliberately for two stated reasons: an instrument-and-
 * wayfinding lettering lineage rather than product-marketing, and no Google
 * Fonts network dependency. Both reasons are honoured here, and the second
 * one is what unblocked the change: `next/font/google` SELF-HOSTS. It fetches
 * at build time and serves the woff2 from this origin, so there is no runtime
 * request to Google and no privacy leak -- the exact objection that pinned the
 * old stack to system fonts does not apply to this API.
 *
 * What the system stack actually produced was the problem. With no webfont
 * loading, every heading resolved to Helvetica Neue Semibold with negative
 * tracking on a near-black ground, which is the single most generic
 * dark-product look there is; the page read as a template regardless of what
 * the copy said.
 *
 * IBM Plex is the right family for this subject on the merits, not as a
 * default. It was commissioned for technical and scientific documentation --
 * the register this product actually lives in -- and it has real character
 * where a neo-grotesque has none: the flat-sided bowls, the notched 'g', the
 * angled terminals that come from drafting lettering. It is emphatically not
 * the current house style of generated design (not Inter, not Space Grotesk,
 * not a geometric). And the project already asked for it: `--font-mono` in
 * globals.css named "IBM Plex Mono" as its second choice all along, so the
 * mono here is the stack's own stated intent finally being loaded rather than
 * silently falling through to SF Mono.
 *
 * Three roles, one superfamily, so they cohere without being flat:
 *   Condensed 600  display -- density and instrument-panel authority
 *   Sans 400/500   running text
 *   Mono 400/500   every measured quantity, label and caption
 * The Condensed/Regular width contrast is what carries hierarchy, rather than
 * size alone, which is why this is not just "one font at three sizes".
 *
 * Weights are enumerated because none of these three is a variable font;
 * asking for a range would silently ship a fallback.
 */
const plexSans = IBM_Plex_Sans({
  subsets: ["latin"],
  weight: ["400", "500", "600"],
  style: ["normal", "italic"],
  variable: "--font-plex-sans",
  display: "swap",
});

const plexCondensed = IBM_Plex_Sans_Condensed({
  subsets: ["latin"],
  weight: ["500", "600", "700"],
  variable: "--font-plex-condensed",
  display: "swap",
});

const plexMono = IBM_Plex_Mono({
  subsets: ["latin"],
  weight: ["400", "500"],
  variable: "--font-plex-mono",
  display: "swap",
});

export const metadata: Metadata = {
  title: "Vivum — See Where the Drug Actually Goes",
  description:
    "Vivum reconstructs real tumor tissue from spatial biology data and simulates drug transport physics through it — penetration, exposure, and delivery, validated against real measured tissue benchmarks. Built for pharma, biotech, and research teams comparing candidate therapies.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html
      lang="en"
      className={`h-full antialiased ${plexSans.variable} ${plexCondensed.variable} ${plexMono.variable}`}
    >
      <body className="min-h-full">{children}</body>
    </html>
  );
}

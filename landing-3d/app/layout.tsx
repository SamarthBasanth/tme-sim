import type { Metadata } from "next";
import "./globals.css";

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
    <html lang="en" className="h-full antialiased">
      <body className="min-h-full">{children}</body>
    </html>
  );
}

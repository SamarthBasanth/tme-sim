import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Tumor Atlas — See the Tumor Before You Treat It",
  description:
    "A cinematic, scroll-driven 3D tumor visualization concept: anatomical context, segmentation, AI-assisted analysis, and treatment planning in one continuous scene.",
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

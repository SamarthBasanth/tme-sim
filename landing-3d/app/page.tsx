import ScrollExperience from "@/components/scroll/ScrollExperience";

export default function Home() {
  return (
    <>
      {/*
        OUTSIDE <main>, not inside it.
        A <header> nested in <main> is not the `banner` landmark — the HTML-AAM
        mapping only gives `banner` to a <header> whose nearest ancestor
        sectioning element is <body>. Inside <main> it degrades to a generic
        group, so the site mark was neither a landmark a screen-reader user
        could jump to nor part of the main content in any meaningful sense; it
        was simply an unnamed box in the middle of the document. Hoisting it
        gives the page the two landmarks it should have had all along, banner
        and main, and changes nothing visually: `.site-header` is
        `position: fixed` and never participated in this flow.
      */}
      <header className="site-header">
        <span className="site-header__mark">CODON</span>
      </header>
      <main>
        <ScrollExperience />
      </main>
    </>
  );
}

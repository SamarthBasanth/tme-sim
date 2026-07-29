import ScrollExperience from "@/components/scroll/ScrollExperience";

export default function Home() {
  return (
    <main>
      <header className="site-header">
        <span className="site-header__mark">TUMOR ATLAS</span>
        <span className="site-header__tagline">3D visualization concept</span>
      </header>
      <ScrollExperience />
    </main>
  );
}

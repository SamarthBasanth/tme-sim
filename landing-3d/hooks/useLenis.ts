"use client";

import { useEffect } from "react";
import gsap from "gsap";
import { ScrollTrigger } from "gsap/ScrollTrigger";
import Lenis from "lenis";

/**
 * Wires Lenis (smooth scroll) into GSAP's own ticker rather than running a
 * second requestAnimationFrame loop — this is the integration pattern Lenis
 * documents for GSAP consumers, and it matters here because ScrollTrigger's
 * `scrub` reads scroll position on GSAP's ticker tick, not the browser's
 * native scroll event. Driving Lenis from the same ticker keeps the two
 * perfectly in sync, which is what makes the scroll-linked timeline feel
 * directly tied to the scrollbar instead of laggy.
 *
 * A no-op when `disabled` (prefers-reduced-motion, or the simplified mobile
 * fallback) — native scrolling is used instead. `onInit` hands the caller
 * the live Lenis instance (e.g. so a CTA button can `lenis.scrollTo(...)`).
 */
export function useLenis(disabled = false, onInit?: (lenis: Lenis) => void): void {
  useEffect(() => {
    if (disabled) return;

    const lenis = new Lenis({
      duration: 1.1,
      smoothWheel: true,
      touchMultiplier: 1.4,
    });

    onInit?.(lenis);
    lenis.on("scroll", ScrollTrigger.update);

    const onTick = (time: number) => {
      lenis.raf(time * 1000);
    };
    gsap.ticker.add(onTick);
    gsap.ticker.lagSmoothing(0);

    return () => {
      gsap.ticker.remove(onTick);
      lenis.destroy();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [disabled]);
}

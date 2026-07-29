"use client";

import { createContext, useContext, useRef, type MutableRefObject, type ReactNode } from "react";

/**
 * A mutable ref, not React state. Scroll progress changes up to once per
 * animation frame; if we pushed it through `useState` every `<TumorModel>`,
 * `<CameraRig>`, and narrative section would re-render on every scroll tick.
 * Instead GSAP's ScrollTrigger writes into `progressRef.current` directly,
 * and R3F components read it inside their own `useFrame` loop — the two
 * render loops (DOM/GSAP and the WebGL canvas) stay decoupled, and nothing
 * in the React tree re-renders because of scrolling.
 */
export type ProgressRef = MutableRefObject<number>;

const ScrollProgressContext = createContext<ProgressRef | null>(null);

export function ScrollProgressProvider({ children }: { children: ReactNode }) {
  const progressRef = useRef(0);
  return (
    <ScrollProgressContext.Provider value={progressRef}>
      {children}
    </ScrollProgressContext.Provider>
  );
}

export function useScrollProgressRef(): ProgressRef {
  const ctx = useContext(ScrollProgressContext);
  if (!ctx) {
    throw new Error("useScrollProgressRef must be used within a ScrollProgressProvider");
  }
  return ctx;
}

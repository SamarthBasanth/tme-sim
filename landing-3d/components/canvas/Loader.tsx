"use client";

import { useEffect, useState } from "react";
import { useProgress } from "@react-three/drei";

/**
 * Visible loading state for the 3D scene, per the spec's requirement. Reads
 * drei's global `useProgress` (backed by three's LoadingManager), so it
 * will correctly track a real `/models/tumor.glb` once one exists.
 *
 * With today's procedural placeholder, nothing ever passes through three's
 * LoadingManager, so `total` (items ever queued) stays 0 forever — `active`
 * and `progress` alone would never report "done" in that case. `total > 0`
 * is what distinguishes "a real asset is loading, wait for it" from
 * "nothing to load, dismiss immediately".
 *
 * Stays mounted briefly after completion to let the fade-out transition
 * finish instead of popping. `dismissed` is only ever set from inside the
 * timeout callback (never synchronously in the effect body), and only ever
 * transitions false -> true.
 */
export default function Loader() {
  const { progress, active, total } = useProgress();
  const [dismissed, setDismissed] = useState(false);
  const isDone = total > 0 ? !active && progress >= 100 : true;

  useEffect(() => {
    if (!isDone) return;
    const timeout = setTimeout(() => setDismissed(true), 420);
    return () => clearTimeout(timeout);
  }, [isDone]);

  if (dismissed) return null;

  return (
    <div className={`scene-loader ${isDone ? "scene-loader--done" : ""}`}>
      <div className="scene-loader__mark" aria-hidden="true" />
      <div className="scene-loader__bar">
        <div className="scene-loader__bar-fill" style={{ width: `${Math.min(100, progress)}%` }} />
      </div>
      <span className="scene-loader__label">
        Loading model — {Math.min(100, Math.round(progress))}%
      </span>
    </div>
  );
}

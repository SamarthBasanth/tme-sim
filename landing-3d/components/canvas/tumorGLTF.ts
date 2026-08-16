"use client";

import { useLoader } from "@react-three/fiber";
import { GLTFLoader, DRACOLoader, MeshoptDecoder } from "three-stdlib";

/**
 * glTF loading for this page, with every decoder resolved from THIS origin.
 *
 * WHY THIS EXISTS INSTEAD OF `useGLTF` FROM DREI
 * ==============================================
 * The page publishes, as a product claim, that it sends no tissue or compound
 * data anywhere and contacts nothing but its own server. drei's `useGLTF`
 * attaches a `DRACOLoader` to every glTF load and, by default, points it at
 *
 *     https://www.gstatic.com/draco/versioned/decoders/1.5.5/
 *
 * `DRACOLoader` is lazy. It fetches the decoder only when it is actually
 * handed a primitive compressed with `KHR_draco_mesh_compression`.
 * `public/models/tumor.glb` declares no `extensionsUsed` and no
 * `extensionsRequired`, so nothing was ever fetched and the claim was true.
 *
 * True BY ACCIDENT. The claim's truth rested on a property of a binary asset
 * that no check was reading. Re-exporting a 660 KB mesh with Draco compression
 * is an ordinary size optimisation; the scene would have kept rendering, no
 * error would have surfaced, and a published privacy claim would have become
 * false in silence. This project's recurring failure is exactly that shape --
 * the engine stays sound while a claim surface drifts out from under it.
 *
 * So the decoders are vendored into `public/draco/` (see the README there:
 * Apache-2.0, 763 341 bytes on disk, zero bytes in the JS bundle and zero
 * bytes on the wire unless a compressed asset is actually loaded) and the
 * loader is pinned at them below.
 *
 * `useGLTF.setDecoderPath("/draco/")` would have repointed drei's loader, and
 * would have been enough at runtime. It is not used, for one reason: drei
 * initialises its `decoderPath` module variable to the gstatic literal, so
 * that string stays in the shipped bundle whatever we assign to it afterwards.
 * A static audit of build output would then have to carry a standing exemption
 * for a live CDN hostname -- and an exemption is a hole that a genuine
 * regression can hide in. Loading through `useLoader` directly means drei's
 * `Gltf` module is never imported, the literal is never emitted, and
 * `scripts/audit-origins.mjs` can hold `www.gstatic.com` as FORBIDDEN with no
 * exemption at all. The guard is stronger because the string is gone.
 *
 * This is the same underlying cache drei uses (`useLoader` from
 * @react-three/fiber) and the same `DefaultLoadingManager`, so `useProgress`
 * in `Loader.tsx` still tracks these loads.
 */

/**
 * Same-origin path serving the vendored Draco decoder. `DRACOLoader` appends
 * `draco_wasm_wrapper.js` + `draco_decoder.wasm` (or `draco_decoder.js` where
 * `WebAssembly` is absent) to this prefix, so all three files live there.
 *
 * A leading slash, not a relative path: the loader resolves it against the
 * document origin, which is precisely the guarantee being made.
 */
export const DRACO_DECODER_PATH = "/draco/";

/**
 * One shared decoder for the page. `DRACOLoader` spins up a worker pool and
 * compiles the WASM module once; a per-load instance would recompile it and
 * re-request the binary.
 */
let dracoLoader: DRACOLoader | null = null;

function configureLoader(loader: GLTFLoader): void {
  if (!dracoLoader) {
    dracoLoader = new DRACOLoader();
    dracoLoader.setDecoderPath(DRACO_DECODER_PATH);
  }
  loader.setDRACOLoader(dracoLoader);

  /*
   * `EXT_meshopt_compression` is the other compression extension a re-export
   * might plausibly turn on, and it is the same trap in a different shape --
   * except that three-stdlib's `MeshoptDecoder` carries its WASM inline as
   * base64 and fetches nothing. Wiring it here keeps parity with drei's
   * default and means a meshopt re-export also just works, from this origin.
   */
  loader.setMeshoptDecoder(
    typeof MeshoptDecoder === "function" ? MeshoptDecoder() : MeshoptDecoder,
  );

  /*
   * No KTX2 / Basis transcoder is attached, deliberately. drei's default for
   * that one points at a jsDelivr CDN, and this page has no textured asset:
   * the tissue section is a vertex-coloured point cloud. Attaching a
   * transcoder "just in case" would re-import exactly the kind of CDN default
   * this module exists to remove. `scripts/audit-origins.mjs` checks the GLB's
   * declared extensions against what is wired up here, so a textured
   * re-export fails the build with a message rather than reaching for a CDN.
   */
}

/**
 * Suspense-based glTF load. Returns the parsed document; callers here use
 * `.scene`.
 */
export function useTumorGLTF(url: string) {
  return useLoader(GLTFLoader, url, configureLoader);
}

/** Warm the cache so the swap from the procedural placeholder has no pop-in. */
export function preloadTumorGLTF(url: string): void {
  useLoader.preload(GLTFLoader, url, configureLoader);
}

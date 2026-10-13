# Vendored Draco decoder — served from this origin, never from a CDN

## Why these files are here

The page claims that it contacts nothing but its own server. Until this
directory existed, that claim was true **by accident**.

`@react-three/drei`'s `useGLTF` attaches a `DRACOLoader` to every glTF load and
points it, by default, at

    https://www.gstatic.com/draco/versioned/decoders/1.5.5/

`DRACOLoader` is lazy: it fetches the decoder only when it is actually handed a
mesh compressed with `KHR_draco_mesh_compression`. `public/models/tumor.glb`
declares no `extensionsUsed` and no `extensionsRequired`, so it is uncompressed,
so the decoder was never fetched, so no request ever left the origin.

That made the truth of a published claim depend on a property of an asset that
nobody was checking. Re-exporting a 660 KB mesh with Draco compression is a
routine size optimisation. The page would have kept working, the scene would
have kept rendering, and the claim would have quietly become false — a silent
failure, which is the exact failure mode this project keeps hitting.

So the decoder is served from here instead, and the loader is pinned at
`/draco/` in `components/canvas/tumorGLTF.ts`. A future Draco-compressed
re-export now simply works, from this origin, with no code change.

## What they are

Copied verbatim from `node_modules/three/examples/jsm/libs/draco/gltf/` as
shipped in **three@0.185.1**. The `gltf/` variant (not the sibling default
build) is the one targeted by the `KHR_draco_mesh_compression` extension, which
is the only Draco that can reach this page.

| file | bytes | fetched when |
| --- | --- | --- |
| `draco_wasm_wrapper.js` | 58 456 | a Draco-compressed glTF is loaded |
| `draco_decoder.wasm` | 192 420 | a Draco-compressed glTF is loaded |
| `draco_decoder.js` | 512 465 | only if `typeof WebAssembly !== "object"` |

Total on disk: 763 341 bytes (745 KiB).

**None of it is in the JavaScript bundle**, and none of it is requested by the
page as it ships today. These are files in `public/`; a browser fetches them
only if `DRACOLoader` asks for them, and `DRACOLoader` only asks when it meets a
compressed primitive. The cost of vendoring, for every visitor to the page as it
currently stands, is zero bytes.

The JS fallback is included even though `WebAssembly` is available in every
browser that can run this scene's WebGL2 renderer. Vendoring two of the three
files would have replaced one latent trap with another — "works everywhere
except the one case we dropped" — and 512 KB of unrequested disk is a poor
reason to leave that behind.

### Integrity

    8bb2952d2ba7d67e1414f8df819410cb0434a666be53f671fff75f68843d76f6  draco_wasm_wrapper.js
    a680d927bed9cb864ddbd63521868891af2bfbe755092761b4837487618df8ac  draco_decoder.wasm
    8625489da79a805f4f2a7d511c3e52d8b4085608a9d2a4d5f4f9de5db0aea04f  draco_decoder.js

`scripts/audit-origins.mjs` re-checks these digests on every build. A decoder
swapped for a different binary is a change to what this origin serves, and it
should not pass silently.

### Licence

Draco is © Google, [Apache License 2.0](https://github.com/google/draco/blob/master/LICENSE).
The files are redistributed unmodified. `LICENSE` in this directory is the
Apache 2.0 text as published by Google with Draco.

## Updating them

Re-copy from `node_modules/three/examples/jsm/libs/draco/gltf/` after a `three`
upgrade, then update the digests above and in `scripts/audit-origins.mjs`. The
decoder must stay version-compatible with the `DRACOLoader` in
`three-stdlib`, which is why it is taken from the installed `three` rather than
pinned independently.

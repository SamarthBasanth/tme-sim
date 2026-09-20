# Deploying `landing-3d` to Vercel

Written from a verification run, not from the Vercel docs. Every setting below
was arrived at by building and serving this app; the "How this was verified"
section at the end lists what was actually executed and what was not.

Scope: **this directory only.** `app/server` (the Python FastAPI Studio) is a
separate deployment target and does not belong on Vercel. Nothing in this
document applies to it.

---

## 1. Project settings

The repository root is `tme-sim/`, which is a Python project with a
`Dockerfile`, a `pyproject.toml` and a `tests/` tree. Vercel's framework
detection runs against the **Root Directory**, so if that is left at the repo
root Vercel will either fail to detect Next.js or try to build the wrong thing.
Setting it correctly is the single most important line in this file.

| Setting | Value | Notes |
| --- | --- | --- |
| **Root Directory** | `landing-3d` | Not the repo root. See above. |
| **Include files outside Root Directory** | **off** (the default) | Verified: nothing under `app/`, `components/`, `lib/`, `hooks/` or `scripts/` references a path outside `landing-3d`. |
| **Framework Preset** | Next.js | Auto-detected once Root Directory is right. |
| **Build Command** | leave as the default (`npm run build`) | Do **not** override it to `next build` — see §2. |
| **Output Directory** | leave blank / default | Next.js on Vercel does not use a static output dir; overriding this breaks the deployment. There is no `output: 'export'` in `next.config.ts`, so `out/` does not exist. |
| **Install Command** | leave as the default (`npm ci`) | `package-lock.json` is committed and was verified to install cleanly from scratch. |
| **Node.js Version** | 22.x | Pinned in-repo, see §3. |
| **Environment Variables** | one, and it matters — see §4 | |

## 2. Do not override the build command

`package.json` defines:

```
"build": "next build && node scripts/audit-origins.mjs"
```

The second half is not a lint step. `scripts/audit-origins.mjs` reads the
shipped build output and fails the build if any external hostname appears in it
that has not been explicitly declared with a written reason. It is what stops a
CDN font, an analytics snippet, a remote image loader, or a re-introduced
`drei/useGLTF` (which hardcodes a `gstatic.com` DRACO decoder URL) from shipping
unnoticed and quietly falsifying the page's own privacy claim.

If you set the Vercel Build Command to `next build`, that guard silently stops
running on every future deploy and nothing tells you. Leave it on the default so
Vercel runs `npm run build`.

**Running the guard on its own needs the env var, and its failure looks like a
regression.** `npm run audit:origins` derives the Studio host's declaration from
`NEXT_PUBLIC_STUDIO_URL`, so invoking it standalone in a shell that does not
export one fails — with a message about an undeclared host, which reads exactly
like the guard catching something. It is not. Either run the full
`npm run build` (which sets nothing but is normally invoked with the variable
present) or pass it explicitly:

```
NEXT_PUBLIC_STUDIO_URL=https://studio.example.org npm run audit:origins
```

Written down because it cost someone a confused minute, and a guard whose
*setup* failure is indistinguishable from a *real* finding is one people learn
to disregard.

## 3. Node version

`next@16.2.12` declares `engines: { node: ">=20.9.0" }`. That range would let
Vercel pick whatever its current default major happens to be, and a default that
moves under you is a classic Vercel-only breakage. So it is pinned explicitly:

- `package.json` → `"engines": { "node": "22.x" }` — **this is the one Vercel
  reads.**
- `.nvmrc` → `22` — for local developers using `nvm`/`fnm`. Vercel does not read
  this file; it is here so the two do not drift.

Verified by execution:

- Node **22.22.0** — `npm ci` + `npm run build` from a clean tree: pass.
- Node **24.13.0** — `npm ci` + `npm run build` from a clean tree: pass.
- Node **26.5.0** (the workstation this was verified on) — `npm run build`,
  `npm run lint`, `npm run audit:network`: pass.

So 24.x is a safe bump if Vercel ever retires 22.x; change the one line in
`engines` and the one line in `.nvmrc` together.

> **Expect a local warning.** If your workstation runs Node 24 or newer, `npm`
> will now print `EBADENGINE ... wanted: {"node":"22.x"}` on install. It is a
> warning, not an error — npm only enforces `engines` when `engine-strict=true`
> is set in an `.npmrc`, and there is no `.npmrc` here. Builds still work; the
> pin exists for Vercel, not for you.

## 3a. Verified against the real Studio, 2026-09-20

Everything above was written against placeholder hosts. This section records a
run against the deployment this page actually points at, because a runbook
verified only with `studio.example.org` has never proved the one thing that
matters — that the built page contains the production URL and not a leftover
loopback default.

```
NEXT_PUBLIC_STUDIO_URL=https://vivum.onrender.com npm run build   -> exit 0
NEXT_PUBLIC_STUDIO_URL=https://vivum.onrender.com npm run lint    -> exit 0
NEXT_PUBLIC_STUDIO_URL=... npm run audit:viewports  -> exit 0 (11 x 2)
NEXT_PUBLIC_STUDIO_URL=... npm run audit:network    -> exit 0
```

Checked in the build output itself, not inferred from the exit code:

| Check | Result |
|---|---|
| `https://vivum.onrender.com` in `.next/` | **5 occurrences** |
| `localhost:8000` in `.next/static/` | **0** |

The second row is the load-bearing one. `lib/env.ts` falls back to
`http://localhost:8000` when the variable is unset, and Next **inlines**
`process.env.NEXT_PUBLIC_STUDIO_URL` as a string at build time — so a build run
without the variable does not fail at runtime in an obvious way, it ships a
page whose primary call to action points at the visitor's own machine. Zero is
the only acceptable number here, and it is worth re-checking on any deploy
where the build log does not show the variable being read.

The live target itself was verified the same day: `/` serves the workspace,
`/studio` serves the Studio, `/index.html` is preserved, and `/api/*` answers
401 without a key.

## 4. Environment variables

**There is exactly one, and the default is wrong in production.**

| Variable | Read by | If unset |
| --- | --- | --- |
| `NEXT_PUBLIC_STUDIO_URL` | `lib/env.ts`, and independently by `scripts/audit-origins.mjs` at build time | **production build FAILS**; `npm run dev` falls back to `http://localhost:8000` |

That is the entire list. Grep confirms no other `process.env` reference exists
in `app/`, `components/`, `lib/` or `hooks/`. No API keys, no secrets, no
database URL, no analytics token. Nothing needs to be added to Vercel except
this one variable.

Two things about it will bite a first-time deployer:

1. **If you do not set it, the production build now FAILS — deliberately.**
   `lib/env.ts` throws, naming the consequence and both places to set the
   value.

   > **Corrected 2026-08-22.** This section previously said the opposite: *"the
   > deploy succeeds and the page's two primary 'Open the Studio' buttons point
   > at `http://localhost:8000`"*. That was true when this runbook was written
   > and was made false, in the same commit that added the runbook, by the
   > refusal in `lib/env.ts`. The old behaviour is exactly why the refusal
   > exists — the deploy was green in every direction (build passed, page
   > rendered, both audits passed) while the main call to action on a public
   > site pointed at the visitor's own machine, and nothing anywhere reported
   > it. `npm run dev` still falls back to localhost so local work needs no
   > setup.

   So the failure mode is now a loud build error instead of a silent bad
   deploy. **Set it before the first production deploy.**

2. **It is a `NEXT_PUBLIC_*` variable, so it is inlined into the JavaScript
   bundle at build time.** Adding or changing it in the Vercel dashboard has no
   effect on an already-built deployment — you must **redeploy** (and the
   redeploy must not reuse the build cache's output). Set it before the first
   production deploy, for the Production environment, and redeploy after any
   change.

Set it to the real Studio URL, e.g. `https://studio.your-domain.tld`. Verified:
building with `NEXT_PUBLIC_STUDIO_URL=https://studio.codon.example` makes the
origin audit re-declare that host as `NAVIGATION` and pass; setting it to a
value that is not a URL fails the build loudly rather than shipping a broken
`href`.

## 5. Do not add Vercel's optional integrations

The page's stated property — that it contacts nothing but its own origin — is
enforced by two scripts in this repo, and is a product claim, not a preference.

**Do not enable, and do not accept the dashboard's prompt to add:**

- **Vercel Web Analytics** (`@vercel/analytics`) — injects a script that beacons
  to `va.vercel-scripts.com`.
- **Vercel Speed Insights** (`@vercel/speed-insights`) — same shape, different
  host.
- Any third-party script, tag manager, chat widget, session recorder, remote
  font, or `next/image` remote loader.

Adding any of these will fail `npm run build` at the origin-audit step, which is
the intended behaviour — the build breaking is the guard working, not a bug to
route around. Do not "fix" it by adding the host to the declaration list.

The fonts are already safe: `next/font/google` downloads the woff2 files **at
build time** and serves them from this origin. Verified at runtime — the only
font requests the page makes are to `/_next/static/media/*.woff2`. The DRACO
decoder is vendored into `public/draco/` for the same reason.

Note for air-gapped or egress-restricted build environments: `next/font/google`
does need network access **during the build** to fetch those woff2 files.
Vercel's build container has it. That is a build-time fetch only; it produces
zero runtime requests.

## 6. Security headers

`next.config.ts` sets them for every path (`source: "/:path*"`), including
`/models/*.glb`, `/draco/*` and the static chunks. Vercel applies
`next.config.ts` headers natively — no `vercel.json` is needed, and none exists.

Headers set (all confirmed present on a live response from the built app):

```
Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-inline'
  'wasm-unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data:
  blob:; font-src 'self'; connect-src 'self'; media-src 'self'; worker-src
  'self' blob:; manifest-src 'self'; object-src 'none'; base-uri 'self';
  form-action 'self'; frame-ancestors 'none'
X-Frame-Options: DENY
X-Content-Type-Options: nosniff
Referrer-Policy: strict-origin-when-cross-origin
Permissions-Policy: camera=(), microphone=(), geolocation=(),
  browsing-topics=(), interest-cohort=()
Strict-Transport-Security: max-age=63072000
```

`poweredByHeader: false` is also set, removing `X-Powered-By: Next.js`.
Confirmed absent from responses.

The three non-obvious choices, each of which is annotated in `next.config.ts`
itself:

- **`'unsafe-inline'` in `script-src`.** The App Router streams its RSC payload
  as inline `<script>` blocks — five of them in the production HTML for `/`.
  Without this, React never hydrates and the scroll experience and canvas never
  start. Next's documented alternative is a per-request nonce from `proxy.ts`,
  which its own docs say **requires dynamic rendering**. This route is currently
  fully prerendered (`○ (Static)`, `x-nextjs-prerender: 1`), which is what lets
  a CDN serve it with no function invocation; nonces would trade that for a
  server render on every view. The page has no forms, no query parameters it
  echoes, no user input and no `dangerouslySetInnerHTML`, so there is no
  injection vector for a nonce to close. If any of that changes, revisit it.
  `'unsafe-eval'` is deliberately **not** present.
- **`'unsafe-inline'` in `style-src`.** The page writes computed `style`
  attributes for scroll-driven heights, progress-bar widths and specimen swatch
  colours — 22 in the pre-hydration HTML, 27 after. There is no formulation that
  keeps those without hashes Next does not emit.
- **`'wasm-unsafe-eval'` + `worker-src blob:`.** For the vendored DRACO decoder.
  Honest caveat: `public/models/tumor.glb` declares no
  `KHR_draco_mesh_compression`, so the decoder is never instantiated today and
  **this allowance was not exercised by the verification run.** It is set so
  that re-exporting the mesh with Draco compression later does not blank the
  scene via a header nobody remembered.

Deliberately **not** set: `upgrade-insecure-requests` (every subresource is a
same-origin relative path, so it would buy nothing, and it would break the local
`http://localhost:8000` Studio link), and `includeSubDomains` / `preload` on
HSTS (the Studio is expected on a sibling subdomain this deployment does not
control, and `preload` is effectively irreversible). Add those from the apex
domain only once every subdomain is known to be HTTPS-only.

## 7. First deploy, start to finish

1. Commit the working tree. `next.config.ts`, `package.json` and `.nvmrc` are
   the files this work touched.
2. Vercel → **Add New… → Project** → import `SamarthBasanth/tme-sim`.
3. **Set Root Directory to `landing-3d`.** Do this on the import screen; it is
   the step that is easy to miss and expensive to debug.
4. Leave Framework Preset (Next.js), Build Command, Output Directory and Install
   Command at their defaults.
5. Add `NEXT_PUBLIC_STUDIO_URL` = the real Studio URL, for **Production**
   (and Preview, if the Studio has a staging host). §4.
6. Deploy.
7. Decline every prompt to add Analytics or Speed Insights. §5.
8. After the deploy, verify from outside:
   ```
   curl -sI https://<your-domain>/ | grep -i -E 'content-security|x-frame|referrer|strict-transport|x-powered-by'
   ```
   Expect the five headers from §6 and **no** `X-Powered-By`. Then open the page
   with devtools on the Network tab, scroll to the bottom, and confirm every
   request is same-origin and the tumour scene renders.

---

## How this was verified

Everything below was executed on this machine against this code, not inferred.

| Check | Result |
| --- | --- |
| `npm run build` (Node 26.5.0, incremental) | pass — `✓ Compiled successfully`, TypeScript clean, `/` and `/_not-found` both `○ (Static)`; origin audit passed with 12 declared hosts, DRACO decoder unmodified, no drei `useGLTF` import |
| `npm run build` from a deleted `.next` | pass — cold build works, `next/font/google` fetched the woff2 files at build time |
| `npm ci` + `npm run build` in a clean copy, Node 24.13.0 | pass — full simulation of what Vercel does |
| `npm ci` + `npm run build` in a clean copy, Node 22.22.0 | pass |
| `npm run lint` | pass — zero output, zero findings |
| `npm run audit:network` | pass — headless Chromium loads the **built** page and scrolls it end to end at 1512×950 and 390×844, in both cinematic and reduced-motion modes; **zero requests to any origin but the page's own** in all four; the audit's own not-vacuous assertions confirmed a live WebGL context, `/models/tumor.glb` 200, the loader dismissed, and the `/draco/` files served from this origin |
| `npx next start` + real browser render | pass — page renders (screenshotted top and bottom), canvas live at 1 436 400 px, `gl.isContextLost() === false`, React hydrated, loader dismissed, **0 console errors, 0 page errors, 0 CSP violations** across desktop / mobile / reduced-motion |
| Security headers on a live response | pass — all six present, `X-Powered-By` absent |
| CSP does not break WebGL or inline styles | pass — this is why the render check above was run *after* the headers were added, in three viewports, with the scene loaded |
| `next.config.ts` Vercel compatibility | pass — no custom server, no `output: 'export'`, no rewrites, no redirects, no `serverExternalPackages`, no filesystem assumptions. It contains `poweredByHeader: false` and one `headers()` entry, both of which Vercel supports natively |
| Build with `NEXT_PUBLIC_STUDIO_URL` set to a real host | pass — origin audit re-declares that host as `NAVIGATION` |
| Build with `NEXT_PUBLIC_STUDIO_URL` set to a non-URL | fails the build loudly, as intended |
| Env var enumeration | one variable, `NEXT_PUBLIC_STUDIO_URL`. No secrets of any kind |

**Not verified, and why:**

- The `'wasm-unsafe-eval'` / `worker-src blob:` CSP allowances — no
  Draco-compressed asset ships, so the decoder never runs. §6.
- Behaviour on Vercel's actual build container. Everything here was verified
  locally, including a from-lockfile install on two Vercel-supported Node
  majors, but a local macOS/arm64 build is not the same machine as a Vercel
  linux/amd64 build container.
- Real-GPU rendering. The headless browser used SwiftShader (software WebGL).
  The scene initialises and draws; how it performs on a real GPU was not
  measured here.

# Vendored typefaces

Four `.woff2` files, committed to this repository, loaded through
`next/font/local` from `app/layout.tsx`. Nothing here is fetched at build time
and nothing here is fetched at runtime. The licence of every file is in this
directory beside it.

## Why these are committed rather than fetched

`next/font/google` self-hosts: it downloads the `.woff2` during `next build`
and serves it from this origin, so it was never a *runtime* privacy problem and
`scripts/network-audit.mjs` passed with it. What it is, is a **build-time
dependency on fonts.gstatic.com**. Three consequences, none of which the
runtime audit can see:

1. The build does not work offline, and did not work in any sandbox without
   egress. "The page makes no third-party requests" was true of the page and
   false of the thing that produces it.
2. The bytes were not pinned. Google re-cuts and re-versions families in place;
   two builds a month apart could ship different outlines, different metrics
   and therefore a different line count, with nothing in the diff.
3. The licence was never distributed with the binary. All three families are
   SIL OFL 1.1, which *requires* the licence to travel with the Font Software.
   Fetching at build time and serving the result satisfied nobody's reading of
   that.

Vendoring fixes all three at a cost of 221 KB in the repository.

## What is here

| File                                         | Family                     | Axes                          | Version | Role                 |
| -------------------------------------------- | -------------------------- | ----------------------------- | ------- | -------------------- |
| `Newsreader-roman-latin.woff2`               | Newsreader                 | `wght` 200–800, `opsz` 6–72   | v26     | `--font-display`     |
| `AtkinsonHyperlegibleNext-roman-latin.woff2` | Atkinson Hyperlegible Next | `wght` 200–800                | v7      | `--font-sans`        |
| `AtkinsonHyperlegibleNext-italic-latin.woff2`| Atkinson Hyperlegible Next | `wght` 200–800, italic        | v7      | `--font-sans` italic |
| `AtkinsonHyperlegibleMono-roman-latin.woff2` | Atkinson Hyperlegible Mono | `wght` 200–800                | v8      | `--font-mono`        |

All four are the **`latin` subset only** — unicode-range
`U+0000-00FF, U+0131, U+0152-0153, U+02BB-02BC, U+02C6, U+02DA, U+02DC, U+0304,
U+0308, U+0329, U+2000-206F, U+20AC, U+2122, U+2191, U+2193, U+2212, U+2215,
U+FEFF, U+FFFD`. That is the same subset `subsets: ["latin"]` was already
shipping, so nothing on the page lost a glyph. Audited rather than assumed:
the rendered sources (`lib/copy.ts`, `components/ui/*`, `components/scroll/*`,
`app/page.tsx`, comments stripped) contain exactly eight non-ASCII characters —
`±` U+00B1, `²` U+00B2, `·` U+00B7, `×` U+00D7, `–` U+2013, `—` U+2014,
`“` U+201C, `”` U+201D — and all eight are inside the range above. Anything
outside it falls through to the system stack in `--font-sans` / `--font-mono` /
`--font-display`, which is exactly what happened before.

No italic is vendored for the display or the mono face, because the stylesheet
asks for neither: the only two italic rules on the page
(`.attribution__footnote`, `.metric-overlay__footnote`) are both set in
`--font-sans`. Adding either file would be 18–130 KB that no selector reaches.

## Provenance, and how to re-obtain these exactly

Each file is the `latin` slice that the Google Fonts CSS v2 API serves to a
`woff2`-capable browser — i.e. byte-for-byte what `next/font/google` was
downloading at build time before this change, for the same families and the
same subset. Re-fetch with:

```
curl -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 \
  (KHTML, like Gecko) Chrome/120.0 Safari/537.36" \
  "https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,200..800"
```

and take the `src:` URL from the `@font-face` block whose `unicode-range`
begins `U+0000-00FF`. Same for
`Atkinson+Hyperlegible+Next:ital,wght@0,200..800;1,200..800` and
`Atkinson+Hyperlegible+Mono:wght@200..800`.

Verify what you got rather than trusting the filename — the axis ranges and the
vertical metrics in the table above are the check, and
`app/globals.css`'s hand-written fallback faces depend on the metrics being
these and not some later re-cut:

```
python3 -c "
from fontTools.ttLib import TTFont            # needs brotli
f = TTFont('Newsreader-roman-latin.woff2')
print(f['head'].unitsPerEm, [(a.axisTag,a.minValue,a.maxValue) for a in f['fvar'].axes])
print(f['hhea'].ascent, f['hhea'].descent, f['hhea'].lineGap)"
```

Measured, per em:

|                            | upm  | ascent | descent | lineGap | `line-height: normal` |
| -------------------------- | ---- | ------ | ------- | ------- | --------------------- |
| Newsreader                 | 2000 | 0.735  | 0.265   | 0       | 1.000 em              |
| Atkinson Hyperlegible Next | 1000 | 0.984  | 0.316   | 0       | 1.300 em              |
| Atkinson Hyperlegible Mono | 1000 | 0.984  | 0.316   | 0       | 1.300 em              |

Those three numbers in the last column are what the `@font-face` fallback
blocks in `app/globals.css` reproduce. If a re-fetch changes any of them, those
blocks are wrong and must be re-measured by the method documented there.

## Licences

All three families are under the **SIL Open Font License 1.1**, which permits
redistribution of the font files, embedded or not, provided the licence travels
with them and the fonts are not sold on their own.

- `OFL-Newsreader.txt` — Copyright 2020 The Newsreader Project Authors
  (Production Type). <https://github.com/productiontype/Newsreader>
- `OFL-AtkinsonHyperlegibleNext.txt` — Copyright 2020–2024 The Atkinson
  Hyperlegible Next Project Authors (Braille Institute of America).
  <https://github.com/googlefonts/atkinson-hyperlegible-next>
- `OFL-AtkinsonHyperlegibleMono.txt` — Copyright 2020–2024 The Atkinson
  Hyperlegible Mono Project Authors.
  <https://github.com/googlefonts/atkinson-hyperlegible-next-mono>

Reserved Font Names: none of the three declares one, so the files may also be
subsetted or re-cut under the same name. Nothing here does.

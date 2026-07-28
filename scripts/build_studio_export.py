"""Build results/studio_export.html — a single self-contained copy of the Studio.

The export is the SAME front end: every file under ``app/web/`` goes in
unmodified except ``js/api.js``, which is swapped for
``scripts/studio_export/api_baked.js`` (identical module surface, payloads read
from an inlined blob instead of fetched). ``js/mock.js`` is deliberately NOT
bundled — nothing in the export is synthetic.

Three things are inlined:

* the CSS, verbatim, plus ``scripts/studio_export/export.css`` for the banner;
* the ES modules, concatenated by the tiny bundler below (there is no npm in
  this project and adding one for a one-shot artifact would be worse than 90
  lines of regex);
* the run payloads: the RunResult JSON exactly as ``GET /api/runs/{id}``
  returns it, and the TMSG / TMSF binaries exactly as the ``/geometry`` and
  ``/field/{name}`` routes serve them, base64'd and deduplicated by content
  hash (the three runs share one reconstruction, so vascular / ECM / oxygen are
  one copy, not three).

Usage::

    ./.venv/bin/python scripts/build_studio_export.py r_00005 r_00006 r_00007
"""
from __future__ import annotations

import base64
import datetime
import hashlib
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

WEB = os.path.join(ROOT, "app", "web")
JS = os.path.join(WEB, "js")
CSS = os.path.join(WEB, "css")
EXPORT_SRC = os.path.join(ROOT, "scripts", "studio_export")
OUT = os.path.join(ROOT, "results", "studio_export.html")
#: Body-level fragment for the Artifact publisher, which supplies its own
#: <!doctype>/<html>/<head>/<body> skeleton. Same bytes of data, same bundle,
#: same banner — only the document wrapper differs.
OUT_ARTIFACT = os.path.join(ROOT, "results", "studio_export_artifact.html")

#: Stable across redeploys: the publisher reads it for the tab and gallery name.
ARTIFACT_TITLE = "tmesim Studio — Glioblastoma Microenvironment Simulator"

FIELD_NAMES = ("vascular_density", "ecm_density", "oxygen", "drug_auc")

#: ``api`` is served from the export tree, not from app/web.
MODULE_OVERRIDES = {
    "api": os.path.join(EXPORT_SRC, "api_baked.js"),
    # A therapy report is a backend JOB plus a backend-rendered HTML document.
    # Nothing about it can be baked, so the panel ships without its controls.
    "ui/reports": os.path.join(EXPORT_SRC, "reports_static.js"),
    # Projects/specimens are mutable server state and a server-side comparison.
    "ui/projects": os.path.join(EXPORT_SRC, "projects_static.js"),
    # Given a ui/ id so its '../util/dom.js' imports resolve like any other
    # panel module; it lives outside app/web because it exists only here.
    "ui/static_shell": os.path.join(EXPORT_SRC, "static_shell.js"),
}

#: Module entry points declared by index.html that are DELIBERATELY not
#: bundled, with the reason. Anything else index.html loads and this file does
#: not list is a build failure — a front end can grow a second entry point at
#: any time, and silently dropping one leaves its markup on screen with nothing
#: wired to it, which is exactly the dead control this export must not ship.
EXCLUDED_ENTRIES = {
    "ui/protocol": (
        "the protocol workspace is entirely server-backed — POST "
        "/api/protocol/{validate,upload,plan,start} and the sweep queue. "
        "Every control in it would 501, so the panel and its toolbar button "
        "are removed by static_shell.js instead of shipped dead."
    ),
}

#: Where the walk starts. Everything else is discovered from the imports, so a
#: front end that grows a module does not need this file edited — and a module
#: that is half-built and not yet reachable from main.js is not dragged in.
ENTRIES = ["main", "ui/static_shell"]


# --------------------------------------------------------------------------- #
# A ~90-line ES-module bundler.
#
# It handles exactly the syntax this codebase uses and refuses anything else,
# loudly. That is the whole trick: a general bundler is hard, a bundler for one
# known dialect is a page of regex that either works or raises.
# --------------------------------------------------------------------------- #
RE_IMPORT_NS = re.compile(
    r"^import\s+\*\s+as\s+([A-Za-z_$][\w$]*)\s+from\s+['\"]([^'\"]+)['\"];?\s*$")
RE_IMPORT_NAMED = re.compile(
    r"^import\s*\{([^}]*)\}\s*from\s+['\"]([^'\"]+)['\"];?\s*$")
RE_EXPORT_LIST = re.compile(r"^export\s*\{([^}]*)\};?\s*$")
RE_EXPORT_DECL = re.compile(
    r"^export\s+(async\s+function|function|class|const|let|var)\s+"
    r"([A-Za-z_$][\w$]*)")
RE_BAD = re.compile(r"^\s*(import|export)\b")


def module_path(mod_id: str) -> str:
    if mod_id in MODULE_OVERRIDES:
        return MODULE_OVERRIDES[mod_id]
    return os.path.join(JS, mod_id + ".js")


def resolve(from_id: str, spec: str) -> str:
    """'./toast.js' seen from 'ui/spec' -> 'ui/toast'."""
    base = os.path.dirname(from_id)
    target = os.path.normpath(os.path.join(base, spec))
    target = target.replace(os.sep, "/")
    if target.endswith(".js"):
        target = target[:-3]
    return target


def transform(mod_id: str, src: str) -> str:
    """Rewrite one ES module into a function body that returns its exports."""
    out, exports = [], []
    for lineno, line in enumerate(src.splitlines(), 1):
        m = RE_IMPORT_NS.match(line)
        if m:
            out.append(f"const {m.group(1)} = __req({resolve(mod_id, m.group(2))!r});")
            continue
        m = RE_IMPORT_NAMED.match(line)
        if m:
            names = []
            for part in m.group(1).split(","):
                part = part.strip()
                if not part:
                    continue
                if " as " in part:
                    a, b = (s.strip() for s in part.split(" as "))
                    names.append(f"{a}: {b}")
                else:
                    names.append(part)
            dep = resolve(mod_id, m.group(2))
            out.append(f"const {{ {', '.join(names)} }} = __req({dep!r});")
            continue
        m = RE_EXPORT_LIST.match(line)
        if m:
            for part in m.group(1).split(","):
                part = part.strip()
                if part:
                    exports.append(part.split(" as ")[-1].strip())
            continue
        m = RE_EXPORT_DECL.match(line)
        if m:
            exports.append(m.group(2))
            out.append(line[len("export "):])
            continue
        if RE_BAD.match(line) and not line.lstrip().startswith(("//", "*", "/*")):
            raise SystemExit(
                f"bundler: unsupported module syntax in {mod_id}.js line "
                f"{lineno}: {line.strip()!r}. Extend build_studio_export.py "
                f"rather than hand-editing the output."
            )
        out.append(line)

    body = "\n".join(out)
    ret = ", ".join(sorted(set(exports)))
    return (
        f"__defs[{mod_id!r}] = function () {{\n{body}\n"
        f"return {{ {ret} }};\n}};\n"
    )


def module_imports(mod_id: str, src: str):
    """Every module this one imports, as bundler ids."""
    deps = []
    for line in src.splitlines():
        m = RE_IMPORT_NS.match(line) or RE_IMPORT_NAMED.match(line)
        if m:
            deps.append(resolve(mod_id, m.group(2)))
    return deps


def walk():
    """Breadth-first over the import graph from ENTRIES. Returns [(id, src)]."""
    seen, order, queue = set(), [], list(ENTRIES)
    while queue:
        mod_id = queue.pop(0)
        if mod_id in seen:
            continue
        seen.add(mod_id)
        path = module_path(mod_id)
        if not os.path.exists(path):
            raise SystemExit(
                f"bundler: {mod_id!r} is imported but {path} does not exist")
        with open(path) as fh:
            src = fh.read()
        order.append((mod_id, src))
        queue.extend(module_imports(mod_id, src))
    return order


def check_api_parity(baked_src: str) -> None:
    """The baked adapter must cover every name app/web/js/api.js exports.

    Otherwise a front end that grows a route silently gets ``undefined`` at
    run time inside a 9 MB file nobody can debug. This turns that into a build
    failure with the missing names listed.
    """
    def names(src):
        out = set()
        for line in src.splitlines():
            m = RE_EXPORT_DECL.match(line)
            if m:
                out.add(m.group(2))
            m = RE_EXPORT_LIST.match(line)
            if m:
                out.update(p.split(" as ")[-1].strip()
                           for p in m.group(1).split(",") if p.strip())
        return out

    with open(os.path.join(JS, "api.js")) as fh:
        real = names(fh.read())
    missing = sorted(real - names(baked_src))
    if missing:
        raise SystemExit(
            "bundler: scripts/studio_export/api_baked.js does not export "
            f"{missing}, which app/web/js/api.js does. Add them (as honest "
            "refusals if they need a server) before shipping the export."
        )


def bundle() -> str:
    with open(MODULE_OVERRIDES["api"]) as fh:
        check_api_parity(fh.read())

    chunks = [
        "/* ---- tmesim Studio front end, bundled from app/web/js/ ---- */\n"
        "var __defs = {}, __cache = {};\n"
        "function __req(id) {\n"
        "  if (__cache[id]) return __cache[id];\n"
        "  var f = __defs[id];\n"
        "  if (!f) throw new Error('module not bundled: ' + id);\n"
        "  return (__cache[id] = f());\n"
        "}\n"
    ]
    mods = walk()
    for mod_id, src in mods:
        chunks.append(f"\n/* ===== {mod_id}.js ===== */\n")
        chunks.append(transform(mod_id, src))
    chunks.append("\n" + "".join(f"__req({e!r});\n" for e in ENTRIES))
    print(f"  modules       {len(mods)}: {', '.join(m for m, _ in mods)}")
    return "".join(chunks)


# --------------------------------------------------------------------------- #
# Payload baking
# --------------------------------------------------------------------------- #
def bake_runs(run_ids):
    """RunResults + deduplicated base64 binaries for the requested runs."""
    from app.server.cache import runs_dir

    runs, blobs, payloads, raw_bytes = [], {}, {}, 0
    for rid in run_ids:
        d = os.path.join(runs_dir(), rid)
        with open(os.path.join(d, "result.json")) as fh:
            doc = json.load(fh)
        if doc.get("status") != "done":
            raise SystemExit(f"{rid} is {doc.get('status')}, not done")
        runs.append(doc)

        pmap = {}
        for name, fname in ([("geometry", "geometry.bin"),
                             ("mechanism", "mechanism.bin"),
                             ("timelapse", "timelapse.bin")]
                            + [(f, f"field_{f}.bin") for f in FIELD_NAMES]):
            path = os.path.join(d, fname)
            if not os.path.exists(path):
                continue
            with open(path, "rb") as fh:
                data = fh.read()
            key = hashlib.sha256(data).hexdigest()[:16]
            if key not in blobs:
                blobs[key] = base64.b64encode(data).decode("ascii")
                raw_bytes += len(data)
            pmap[name] = key
        payloads[rid] = pmap

    runs.sort(key=lambda r: r.get("created_at") or 0, reverse=True)
    return runs, blobs, payloads, raw_bytes


def main(argv) -> int:
    from app.server.api import default_spec, drug_rows, validate_spec
    from app.server.cache import list_sections

    run_ids = argv[1:]
    if not run_ids:
        # Whatever scripts/export_studio_runs.py produced last.
        with open(os.path.join(ROOT, "results", "_export_run_ids.json")) as fh:
            run_ids = json.load(fh)
    runs, blobs, payloads, raw_bytes = bake_runs(run_ids)

    sections = list_sections()
    drugs = drug_rows()
    defaults = default_spec()

    # The real verdict for each precomputed spec, recorded now so the export can
    # show it instead of guessing.
    def canon(o):
        return json.dumps(o, sort_keys=True, separators=(",", ":"))

    validate = {canon(r["spec"]): validate_spec(r["spec"]) for r in runs}

    flagship = runs[-1]
    spec = flagship["spec"]
    rep = flagship["report"]
    section_id = spec["section"]
    section = next((s for s in sections if s["id"] == section_id), {})
    n_cells = int(rep["n_cells"])
    grid = rep["grid"]

    meta = {
        "section_id": section_id,
        "section_label": section.get("label") or section_id,
        "panel": section.get("panel"),
        "accession": section.get("accession"),
        "dx_um": spec["dx_um"],
        "sensitivity_mode": spec["sensitivity_mode"],
        "plasma_peak_uM": spec["pk"]["plasma_peak"],
        "total_time_h": round(spec["total_time_s"] / 3600.0, 2),
        "n_cells_total": n_cells,
        "n_cells_baked": n_cells,
        "subsample_stride": 1,
        "sampling_label": f"1:1 · {n_cells:,} of {n_cells:,} cells",
        "sampling_note": (
            f"NO SUBSAMPLING. All {n_cells:,} segmented cells of the section are "
            f"in this file, at full float32 precision — the geometry payload is "
            f"the byte-identical TMSG the server would serve."
        ),
        "field_rows": grid[0],
        "field_cols": grid[1],
        "payload_mib": round(raw_bytes / 1048576.0, 2),
        # The tab / gallery name. Baked rather than composed in the front end,
        # because the Artifact publisher reads the <title> ELEMENT while the
        # browser tab shows document.title — if static_shell.js rewrote one and
        # not the other they would disagree, and the gallery card would drift
        # from the tab on every redeploy.
        "document_title": None,   # filled per-variant below
        "built_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "repo": "tme-sim",
        "commands": [
            "./.venv/bin/python -m pip install -e '.[app]'   # fastapi + uvicorn",
            "./.venv/bin/python -m app.server                # http://127.0.0.1:8000",
        ],
    }

    baked = {
        "export_meta": meta,
        "sections": sections,
        "drugs": drugs,
        "defaults": defaults,
        "runs": runs,
        "validate": validate,
        "payloads": payloads,
        "blobs": blobs,
    }

    # Read the stylesheet list off index.html rather than keeping a copy of it
    # here, so a front end that adds a sheet cannot ship an unstyled export.
    with open(os.path.join(WEB, "index.html")) as fh:
        index = fh.read()
    sheets = re.findall(r'<link rel="stylesheet" href="css/([^"]+)">', index)
    if not sheets:
        raise SystemExit("no <link rel=stylesheet> found in app/web/index.html")
    css = []
    for name in sheets:
        with open(os.path.join(CSS, name)) as fh:
            css.append(f"/* ===== css/{name} ===== */\n" + fh.read())
    with open(os.path.join(EXPORT_SRC, "export.css")) as fh:
        css.append("/* ===== export.css ===== */\n" + fh.read())

    # Every module entry point index.html declares must be accounted for.
    declared = [m[len("js/"):-len(".js")] for m in
                re.findall(r'<script type="module" src="([^"]+)"></script>', index)]
    unhandled = [d for d in declared
                 if d not in ENTRIES and d not in EXCLUDED_ENTRIES]
    if unhandled:
        raise SystemExit(
            f"bundler: app/web/index.html loads module entry point(s) {unhandled} "
            "that this build neither bundles nor excludes. Add them to ENTRIES "
            "(if they work without a server) or to EXCLUDED_ENTRIES with a "
            "reason — do not ship their markup with nothing wired to it."
        )
    for d in declared:
        if d in EXCLUDED_ENTRIES:
            print(f"  excluded      {d}: {EXCLUDED_ENTRIES[d].split(chr(10))[0]}")

    body = index.split("<body>", 1)[1].rsplit("</body>", 1)[0]
    # Drop the module <script> tags; the bundle replaces them.
    body = re.sub(r'<script type="module"[^>]*></script>\s*', "", body)

    # The shipped favicon is an inline SVG, which must carry the SVG xmlns —
    # an http:// namespace identifier that is never fetched but does put the
    # substring "http://" inside an href. A PNG data URI of the same mark keeps
    # the export's "zero http(s) in src/href" property literally true.
    with open(os.path.join(EXPORT_SRC, "favicon.b64")) as fh:
        icon = fh.read().strip()
    standalone_title = (
        f"tmesim Studio — static export · {section_id} · dx {spec['dx_um']:g} um")

    def payload_for(title: str) -> str:
        baked["export_meta"]["document_title"] = title
        return json.dumps(baked, separators=(",", ":"), allow_nan=False)

    style_block = "<style>\n" + "\n".join(css) + "\n</style>\n"
    bundle_js = bundle()

    def data_block(title: str) -> str:
        # </script> cannot appear inside the JSON; escape defensively.
        return (
            '<script id="tmesim-baked" type="application/json">'
            + payload_for(title).replace("</", "<\\/")
            + "</script>\n"
            "<script>\n"
            'window.__TMESIM_BAKED = JSON.parse('
            'document.getElementById("tmesim-baked").textContent);\n'
            "</script>\n"
            "<script>\n" + bundle_js + "\n</script>\n"
        )

    # ---------------------------------------------------------------- #
    # 1. Standalone: a complete document you can double-click.
    # ---------------------------------------------------------------- #
    standalone = (
        "<!doctype html>\n"
        '<html lang="en" data-theme="dark">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, '
        'viewport-fit=cover">\n'
        '<meta name="color-scheme" content="dark light">\n'
        f"<title>{standalone_title}</title>\n"
        '<meta name="description" content="Single-file static export of tmesim '
        'Studio. Precomputed runs, no server.">\n'
        f'<link rel="icon" href="{icon}">\n'
        + style_block
        + "</head>\n<body>\n"
        + body
        + "\n" + data_block(standalone_title)
        + "</body>\n</html>\n"
    )

    # ---------------------------------------------------------------- #
    # 2. Artifact fragment: body-level content ONLY.
    #
    # No <!doctype>, <html>, <head> or <body> — the publisher wraps this in
    # its own skeleton, and a nested document element would be dropped by the
    # parser along with anything the author assumed it carried. Everything the
    # standalone build puts in <head> therefore has to survive at body level:
    #
    #   * <title>   — legal anywhere; the parser routes it through the "in
    #                 head" rules, so document.title is set either way. The
    #                 publisher reads it for the tab and the gallery card.
    #   * <style>   — legal in body; scope is the document, not the position.
    #   * favicon   — dropped. The publisher sets the icon itself.
    #   * charset / viewport / color-scheme — the publisher's own <head> owns
    #                 them; a second copy in <body> would be ignored anyway.
    #
    # data-theme cannot be stamped on <html> from here, and does not need to
    # be: css/tokens.css resolves an unset theme from prefers-color-scheme,
    # and main.js's initTheme() sets the attribute on documentElement at boot.
    # ---------------------------------------------------------------- #
    artifact = (
        f"<title>{ARTIFACT_TITLE}</title>\n"
        + style_block
        + body
        + "\n" + data_block(ARTIFACT_TITLE)
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for path, doc in ((OUT, standalone), (OUT_ARTIFACT, artifact)):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(doc)

    print(f"  runs          {[r['run_id'] for r in runs]}")
    print(f"  cells         {n_cells:,} (stride 1, no subsampling)")
    print(f"  field grid    {grid[0]} x {grid[1]} @ {spec['dx_um']:g} um")
    print(f"  blobs         {len(blobs)} unique, {raw_bytes:,} B raw")
    print(f"  stylesheets   {sheets}")
    for path, doc in ((OUT, standalone), (OUT_ARTIFACT, artifact)):
        size = os.path.getsize(path)
        # Static CSP surface of the emitted markup (the <script> bodies are JS,
        # not markup, so attribute regexes must not be run over them).
        markup = doc.split('<script id="tmesim-baked"')[0]
        inline_style = re.findall(r"\sstyle\s*=\s*[\"']", markup)
        inline_on = re.findall(r"\son[a-z]+\s*=\s*[\"']", markup)
        print(f"wrote {path}")
        print(f"  size          {size:,} B  ({size / 1048576.0:.2f} MiB)")
        print(f"  inline style= {len(inline_style)}   inline on*= {len(inline_on)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

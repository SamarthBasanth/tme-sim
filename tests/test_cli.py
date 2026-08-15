"""Tests for ``tmesim.cli`` — the ``codon`` console script / ``python -m
tmesim.cli`` entry point.

Two kinds of coverage:

1. **Argparse plumbing** — every subcommand is registered, ``--help`` works
   for each without crashing, and ``codon`` is what the tool calls itself now
   that it is a real installed command (``pyproject.toml``'s
   ``[project.scripts]``). The pre-rename ``vivum`` name stays registered as an
   alias so the commands quoted in ``results/*.md`` still run.
2. **``compare-drugs`` end-to-end** — reconstructs the real (open, cached)
   Xenium GBM section at a coarse grid and compares two real drugs on it,
   checking the printed table and the optional CSV output. Uses the same
   coarse ``dx`` as ``tests/conftest.py``'s ``real_micro`` fixture so this
   stays fast (a few seconds), and skips cleanly if the Xenium cache this
   repository ships with is not present.

``studio`` is tested with ``app.server.__main__.main`` monkeypatched rather
than actually binding a socket — this is a CLI-plumbing test, not a Studio
integration test (see ``tests/test_app_server.py`` for that).
"""
from __future__ import annotations

import csv
import json
import os
import pathlib

import pytest

import tmesim.cli as cli

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
XENIUM_GBM_CACHE = os.path.join(_REPO, "data_cache", "xenium_gbm")

# Coarse grid: matches conftest.py's real_micro fixture, so reconstruction
# finishes in a fraction of a second and the whole test file stays fast.
_FAST_DX = 80.0


def _csv_rows(path) -> list:
    """Data rows of an exported comparison CSV.

    The file now opens with a ``#`` provenance block (section, dx, window,
    alpha, delivery law, NOT FOR CLINICAL USE) that travels inside the artefact
    rather than beside it, so the reader skips comment lines — exactly what
    ``pandas.read_csv(path, comment='#')`` does. The same block is written as a
    ``.provenance.json`` sidecar for readers that cannot.
    """
    with open(path, newline="") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines))


def _skip_if_no_real_data() -> None:
    if not os.path.isdir(XENIUM_GBM_CACHE):
        pytest.skip(f"real Xenium GBM cache not present at {XENIUM_GBM_CACHE}")


# --------------------------------------------------------------------------- #
# Argparse plumbing
# --------------------------------------------------------------------------- #
class TestArgparsePlumbing:
    def test_all_subcommands_registered(self):
        parser = cli.build_parser()
        sub = next(a for a in parser._subparsers._group_actions
                   if hasattr(a, "choices"))
        assert set(sub.choices) == {
            "reconstruct", "simulate", "compare-drugs", "validate", "ingest",
            "molecules", "molecular-run", "acquire", "sections", "studio",
        }

    @pytest.mark.parametrize("args", [
        [],
        ["reconstruct", "--help"],
        ["simulate", "--help"],
        ["compare-drugs", "--help"],
        ["validate", "--help"],
        ["ingest", "--help"],
        ["molecules", "--help"],
        ["molecular-run", "--help"],
        ["acquire", "--help"],
        ["sections", "--help"],
        ["studio", "--help"],
    ])
    def test_help_exits_cleanly(self, args, capsys):
        with pytest.raises(SystemExit) as exc:
            cli.build_parser().parse_args(args)
        # argparse exits 0 for --help, 2 for a genuinely missing subcommand
        # (bare `codon` with no args, since the subparser is required).
        assert exc.value.code in (0, 2)

    def test_prog_name_is_codon(self):
        # cosmetic, but this is the name pyproject.toml registers as the
        # console script, and the QUICKSTART examples assume it.
        assert cli.build_parser().prog == "codon"

    def test_both_console_scripts_point_at_the_cli(self):
        """`codon` is the name; `vivum` is kept as a working alias.

        The rename happened after ~200 dated files in ``results/`` had already
        quoted `vivum ...` commands as their reproduce instructions. Dropping
        the old entry point would break every one of them, so both names are
        registered and both resolve to the same ``main``.
        """
        import re
        root = pathlib.Path(__file__).resolve().parents[1]
        text = (root / "pyproject.toml").read_text(encoding="utf-8")
        block = text.split("[project.scripts]", 1)[1].split("\n[", 1)[0]
        entries = dict(re.findall(r'^(\w+)\s*=\s*"([^"]+)"', block, re.M))
        assert entries.get("codon") == "tmesim.cli:main"
        assert entries.get("vivum") == "tmesim.cli:main"

    def test_compare_drugs_defaults(self):
        parser = cli.build_parser()
        args = parser.parse_args(["compare-drugs"])
        assert args.library == "gbm"
        assert args.drugs is None
        assert args.csv is None
        assert args.func is cli.cmd_compare_drugs

    def test_studio_defaults(self):
        parser = cli.build_parser()
        args = parser.parse_args(["studio"])
        assert args.host == "127.0.0.1"
        assert args.port == 8000
        assert args.reload is False
        assert args.func is cli.cmd_studio


# --------------------------------------------------------------------------- #
# compare-drugs, end-to-end on the real (cached) Xenium GBM section
# --------------------------------------------------------------------------- #
class TestCompareDrugsEndToEnd:
    def test_two_real_drugs_prints_a_sorted_table(self, capsys):
        _skip_if_no_real_data()
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "temozolomide,cisplatin",
        ])
        assert rc == 0
        out = capsys.readouterr().out
        assert "drug comparison" in out
        assert "temozolomide" in out
        assert "cisplatin" in out

    def test_writes_a_csv_with_one_row_per_drug(self, tmp_path, capsys):
        _skip_if_no_real_data()
        csv_path = tmp_path / "comparison.csv"
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "temozolomide,cisplatin", "--csv", str(csv_path),
        ])
        assert rc == 0
        assert csv_path.exists()
        rows = _csv_rows(csv_path)
        assert {r["drug"] for r in rows} == {"temozolomide", "cisplatin"}
        for r in rows:
            assert 0.0 <= float(r["kill_fraction"]) <= 1.0
            assert float(r["ic50_uM"]) > 0.0

    def test_the_exported_csv_carries_its_own_provenance(self, tmp_path):
        """The CSV is the only machine-readable export, so it is the artefact
        most likely to end up in a slide. It used to be eight columns and
        nothing else: no section, no dx, no window, no alpha, no delivery law,
        no 'not for clinical use'. Everything the run knows travels with it."""
        _skip_if_no_real_data()
        csv_path = tmp_path / "comparison.csv"
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "temozolomide", "--csv", str(csv_path),
        ])
        assert rc == 0

        # (a) in the file itself, where it cannot be separated from the numbers
        head = "\n".join(l for l in csv_path.read_text().splitlines()
                         if l.startswith("#"))
        assert "NOT FOR CLINICAL USE" in head
        for key in ("section_source", "dx_um", "exposure_window_s", "alpha",
                    "delivery_law", "sensitivity_mode", "seed"):
            assert f"{key}:" in head, f"{key} missing from the CSV header"
        # the law disclosure, not merely the word
        assert "S = PS*" in head
        assert "penetration_index is NOT in this file" in head

        # (b) machine-readable, for a reader that cannot skip comments
        side = json.loads((csv_path.parent / (csv_path.name +
                                              ".provenance.json")).read_text())
        assert side["exposure_window_s"] > 0
        assert side["alpha"] > 0
        assert side["dx_um"] == pytest.approx(_FAST_DX)
        assert side["delivery_law"].startswith("S = PS*")
        assert "NOT FOR CLINICAL USE" in side["not_for_clinical_use"]

        # (c) and the data is still one row per drug behind the comments
        rows = _csv_rows(csv_path)
        assert [r["drug"] for r in rows] == ["temozolomide"]

    def test_the_csv_carries_coverage_and_not_the_retired_index(self, tmp_path):
        """`penetration_index` is retired (results/delivery_metrics.md §4): it is
        a scale-invariant shape statistic that reports the section's vascular
        geometry, which every drug shares. `therapeutic_coverage` @1x IC50 was
        always computable and was never exported. They have swapped places."""
        _skip_if_no_real_data()
        csv_path = tmp_path / "comparison.csv"
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "temozolomide", "--csv", str(csv_path),
        ])
        assert rc == 0
        rows = _csv_rows(csv_path)
        assert "therapeutic_coverage_1x" in rows[0]
        assert "penetration_index" not in rows[0]
        assert 0.0 <= float(rows[0]["therapeutic_coverage_1x"]) <= 1.0

    def test_it_prints_the_window_and_the_delivery_law_it_ran_under(self, capsys):
        """The ranking is LAW-DEPENDENT and window-dependent, and this command
        — the core use case — printed neither. The disclosure comes from the one
        implementation that reads PDEConfig at call time, so it cannot go stale.
        """
        _skip_if_no_real_data()
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX), "--drugs", "temozolomide",
        ])
        assert rc == 0
        out = capsys.readouterr().out
        assert "EXPOSURE WINDOW" in out
        assert "1800 s" in out                 # this CLI's own default...
        assert "43200 s" in out                # ...named against the library's
        assert "Delivery law:" in out
        assert "S = PS*" in out
        assert "LAW" in out
        # and the retired metric is no longer a COLUMN of the printed table
        # (it survives only in the footnote that says it was retired, which is
        # history, not a number anyone can read off)
        header = next(l for l in out.splitlines() if "IC50(uM)" in l)
        assert "cov@IC50" in header
        assert "pen_idx" not in header

    def test_unknown_drug_is_skipped_not_fatal(self, capsys):
        # The molecular libraries are actively growing (see docs/QUICKSTART.md);
        # one bad/unwired name must not sink the whole comparison.
        _skip_if_no_real_data()
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "temozolomide,not-a-real-drug",
        ])
        assert rc == 0
        out = capsys.readouterr().out
        assert "SKIPPED" in out
        assert "temozolomide" in out

    def test_all_requested_drugs_unknown_is_a_clean_failure(self, capsys):
        _skip_if_no_real_data()
        rc = cli.main([
            "compare-drugs", "--dx", str(_FAST_DX),
            "--drugs", "not-a-real-drug",
        ])
        assert rc == 1

    def test_missing_breast_sample_default_fails_before_touching_disk(
        self, monkeypatch, capsys
    ):
        # No --sample given: compare-drugs resolves the default breast cache
        # path itself and checks it up front, rather than letting a deep
        # FileNotFoundError surface from inside the reconstruction pipeline.
        monkeypatch.setattr(cli, "_BREAST_SAMPLE_PATH", "data_cache/nope")
        rc = cli.main(["compare-drugs", "--library", "breast"])
        assert rc == 1
        out = capsys.readouterr().out
        assert "not present" in out
        assert "docs/QUICKSTART.md" in out

    def test_explicit_bad_sample_fails_with_a_clear_message(self, capsys):
        # An explicit --sample is trusted and passed straight through; if it is
        # wrong, the underlying loader's error should name the path it looked
        # in rather than an unrelated default.
        rc = cli.main([
            "compare-drugs", "--library", "breast",
            "--sample", "data_cache/definitely_not_a_real_cache_dir",
        ])
        assert rc == 1
        # main()'s top-level guard prints the failure + traceback to stderr.
        err = capsys.readouterr().err
        assert "data_cache/definitely_not_a_real_cache_dir" in err


# --------------------------------------------------------------------------- #
# studio: argv construction only (no real server bind)
# --------------------------------------------------------------------------- #
class TestStudioDelegatesToAppServer:
    def test_delegates_with_constructed_argv(self, monkeypatch):
        fastapi = pytest.importorskip("fastapi", reason="needs the 'app' extra")
        pytest.importorskip("uvicorn", reason="needs the 'app' extra")

        captured = {}

        def fake_main(argv):
            captured["argv"] = argv
            return 0

        import app.server.__main__ as server_main
        monkeypatch.setattr(server_main, "main", fake_main)

        rc = cli.main([
            "studio", "--host", "0.0.0.0", "--port", "9999",
            "--reload", "--cache-size", "2",
        ])
        assert rc == 0
        assert captured["argv"] == [
            "--host", "0.0.0.0", "--port", "9999", "--log-level", "info",
            "--reload", "--cache-size", "2",
        ]

    def test_missing_app_extra_gives_a_clear_message_not_a_traceback(
        self, monkeypatch, capsys
    ):
        import builtins

        real_import = builtins.__import__

        def blocked_import(name, *a, **kw):
            if name == "app.server.__main__" or name.startswith("app.server"):
                raise ImportError("simulated: fastapi not installed")
            return real_import(name, *a, **kw)

        monkeypatch.setattr(builtins, "__import__", blocked_import)
        rc = cli.main(["studio"])
        assert rc == 1
        out = capsys.readouterr().out
        assert "pip install -e '.[app]'" in out


# --------------------------------------------------------------------------- #
# acquire / sections — the first-run contract
#
# The failure these lock down: `codon acquire` used to fetch four TABULAR
# sources, report 4/4 ok, and leave every tissue section unavailable, with
# nothing in its output connecting those two facts. A user ran the command the
# README recommends and was returned to the start. Every assertion here is
# about SAYING SO; none of them touches the network.
# --------------------------------------------------------------------------- #
class TestAcquirePlanAndSections:
    def test_acquire_list_downloads_nothing_and_names_the_split(self, capsys):
        rc = cli.main(["acquire", "--list"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "nothing was downloaded" in out
        # the tabular sources are named, AND disclaimed
        for source, _what in cli._TABULAR_SOURCES:
            assert source in out
        assert "NONE of those four makes a tissue section runnable" in out
        # the section that does make one runnable is named, with its size
        assert "xenium_319" in out
        assert "16.6 MB" in out

    def test_acquire_list_sections_none_says_nothing_will_be_runnable(self, capsys):
        rc = cli.main(["acquire", "--list", "--sections", "none"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "NO spatial section" in out
        assert "will not make" in out

    def test_acquire_rejects_an_unknown_section_by_name(self, capsys):
        rc = cli.main(["acquire", "--sections", "not_a_section"])
        out = capsys.readouterr().out
        assert rc == 2
        assert "unknown section id(s): not_a_section" in out
        assert "xenium_319" in out  # tells you what IS known

    def test_every_unavailable_section_carries_an_action(self, capsys):
        """The point of `codon sections`: a remedy, not just a diagnosis."""
        from tmesim.data import section_catalog as SC

        rc = cli.main(["sections"])
        out = capsys.readouterr().out
        # rc is 1 only when NOTHING is runnable; both are legitimate here, since
        # this repo's data_cache may or may not be populated on a given machine.
        assert rc in (0, 1)
        for sid, rec in SC.SECTION_RECIPES.items():
            assert sid in out
            if f"[UNAVAILABLE] {sid}" not in out:
                continue
            assert "TO FIX IT:" in out
            if rec["mode"] == SC.AUTO:
                assert f"codon acquire --sections {sid}" in out
            else:
                assert "NOT auto-fetchable" in out

    def test_sections_json_is_machine_readable(self, capsys):
        rc = cli.main(["sections", "xenium_319", "--json"])
        doc = json.loads(capsys.readouterr().out)
        assert rc in (0, 1)
        assert doc["total"] == 1
        sec = doc["sections"][0]
        assert sec["id"] == "xenium_319"
        assert isinstance(sec["available"], bool)
        assert sec["verified"]  # every recipe records HOW it was verified

    def test_unknown_section_id_is_refused_not_ignored(self, capsys):
        rc = cli.main(["sections", "nope"])
        assert rc == 2
        assert "unknown section id(s): nope" in capsys.readouterr().out


class TestSectionCatalog:
    def test_ids_match_the_studio_registry_when_it_is_importable(self):
        """The catalog must not drift from ``app/server/cache.py``'s SECTIONS.

        A remedy for a section id the Studio does not have is a remedy nobody
        can use, and a Studio section with no catalog entry is exactly the dead
        end this module exists to remove.
        """
        cache = pytest.importorskip("app.server.cache")
        from tmesim.data import section_catalog as SC

        assert list(SC.SECTION_RECIPES) == list(cache.SECTIONS)

    def test_recipes_agree_with_the_registry_on_where_the_data_lives(self):
        cache = pytest.importorskip("app.server.cache")
        from tmesim.data import section_catalog as SC

        for sid, rec in SC.SECTION_RECIPES.items():
            assert rec["cache_dir"] == cache.SECTIONS[sid]["cache_dir"]
            assert rec["platform"] == cache.SECTIONS[sid]["source"]
            assert SC._SAMPLE_ARG[sid] == cache.SECTIONS[sid]["sample"]

    def test_auto_recipes_declare_a_byte_count_for_every_file(self):
        """Byte counts are the only integrity check on a streamed download."""
        from tmesim.data import section_catalog as SC

        for sid in SC.auto_section_ids():
            rec = SC.SECTION_RECIPES[sid]
            assert rec["files"]
            assert rec["bytes"] == sum(n for _u, _b, n in rec["files"])
            for url, basename, n in rec["files"]:
                assert url.startswith("https://")
                assert n > 0 and basename

    def test_manual_recipes_carry_commands_not_a_document_reference(self):
        from tmesim.data import section_catalog as SC

        for sid in SC.manual_section_ids():
            remedy = "\n".join(SC.remedy_for(sid))
            assert "curl" in remedy or "codon ingest" in remedy

    def test_resolve_section_arg(self):
        from tmesim.data import section_catalog as SC

        assert SC.resolve_section_arg(None) == ([SC.DEFAULT_SECTION], [])
        assert SC.resolve_section_arg("none") == ([], [])
        assert SC.resolve_section_arg("all") == (SC.auto_section_ids(), [])
        assert SC.resolve_section_arg("xenium_319,nope") == (["xenium_319"],
                                                             ["nope"])

    def test_fetch_section_refuses_a_manual_section_without_touching_the_network(
        self, tmp_path
    ):
        from tmesim.data import section_catalog as SC

        out = SC.fetch_section("visium_gbm1", cache_root=str(tmp_path))
        assert out["status"] == "skipped"
        assert "codon sections" in out["note"]
        assert not list(tmp_path.iterdir())

    def test_human_bytes_never_invents_a_number(self):
        from tmesim.data import section_catalog as SC

        assert SC.human_bytes(None) == "unknown size"
        assert SC.human_bytes(0) == "unknown size"
        assert SC.human_bytes(16_598_845) == "16.6 MB"

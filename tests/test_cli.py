"""Tests for ``tmesim.cli`` — the ``vivum`` console script / ``python -m
tmesim.cli`` entry point.

Two kinds of coverage:

1. **Argparse plumbing** — every subcommand is registered, ``--help`` works
   for each without crashing, and ``vivum`` is what the tool calls itself now
   that it is a real installed command (``pyproject.toml``'s
   ``[project.scripts]``). The superseded ``codon`` name stays registered as an
   alias so the commands quoted in ``results/*.md`` still run.
2. **``simulate`` end-to-end, and ``compare-drugs`` staying gone** —
   reconstructs the real (open, cached) Xenium GBM section at a coarse grid and
   maps ONE compound on it, checking the printed report. Uses the same coarse
   ``dx`` as ``tests/conftest.py``'s ``real_micro`` fixture so this stays fast
   (a few seconds), and skips cleanly if the Xenium cache this repository ships
   with is not present.

   ``compare-drugs`` was REMOVED on 2026-08-21 (see ``docs/QUICKSTART.md`` §1).
   It reconstructed one section and printed several compounds ranked against
   each other, and that cross-compound ordering was measured and does not hold
   (``results/pk_ranking_robustness.md``: exposure spans 113.1x and the
   published GBM order breaks at 1.53x; ``results/capacity_ranking_robustness.md``:
   the breast top-1 pick changes under this repo's own capacity rule). Every
   assertion that used to pin the command's behaviour is still here, converted
   to pin the new contract: the command does not exist, the parser refuses it
   and names ``simulate``, no ranked CSV is written by anything, and the
   disclosures the CSV used to carry still exist at their single source.

``studio`` is tested with ``app.server.__main__.main`` monkeypatched rather
than actually binding a socket — this is a CLI-plumbing test, not a Studio
integration test (see ``tests/test_app_server.py`` for that).
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import pathlib

import pytest

import tmesim.cli as cli

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
XENIUM_GBM_CACHE = os.path.join(_REPO, "data_cache", "xenium_gbm")
QUICKSTART = pathlib.Path(_REPO) / "docs" / "QUICKSTART.md"

# Coarse grid: matches conftest.py's real_micro fixture, so reconstruction
# finishes in a fraction of a second and the whole test file stays fast.
_FAST_DX = 80.0


def _skip_if_no_real_data() -> None:
    if not os.path.isdir(XENIUM_GBM_CACHE):
        pytest.skip(f"real Xenium GBM cache not present at {XENIUM_GBM_CACHE}")


@pytest.fixture(scope="module")
def one_map(tmp_path_factory):
    """ONE real ``vivum simulate`` run, shared by the tests below.

    This is the product now: one real tissue, one compound. The run is done
    once at module scope because reconstruct+solve is the expensive part, and
    several converted assertions read the same output. ``--out`` points at a
    tmp dir so the suite never writes into ``results/``.
    """
    _skip_if_no_real_data()
    out_dir = tmp_path_factory.mktemp("one_map")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = cli.main(["simulate", "--dx", str(_FAST_DX),
                       "--drug", "temozolomide", "--out", str(out_dir)])
    return rc, buf.getvalue(), out_dir


# --------------------------------------------------------------------------- #
# Argparse plumbing
# --------------------------------------------------------------------------- #
class TestArgparsePlumbing:
    def test_all_subcommands_registered(self):
        parser = cli.build_parser()
        sub = next(a for a in parser._subparsers._group_actions
                   if hasattr(a, "choices"))
        assert set(sub.choices) == {
            "reconstruct", "simulate", "validate", "ingest",
            "molecules", "molecular-run", "acquire", "sections", "studio",
        }
        # was: "compare-drugs" is one of the registered subcommands.
        # now:  it is not, and nothing re-registers it under another name.
        assert "compare-drugs" not in sub.choices
        assert not any("compare" in c for c in sub.choices)

    @pytest.mark.parametrize("args", [
        [],
        ["reconstruct", "--help"],
        ["simulate", "--help"],
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
        # (bare `vivum` with no args, since the subparser is required).
        assert exc.value.code in (0, 2)

    def test_prog_name_is_vivum(self):
        # cosmetic, but this is the name pyproject.toml registers as the
        # console script, and the QUICKSTART examples assume it.
        assert cli.build_parser().prog == "vivum"

    def test_both_console_scripts_point_at_the_cli(self):
        """`vivum` is the name; `codon` is kept as a working alias.

        The product name went Vivum -> Codon (2026-08-14) -> Vivum
        (2026-08-22). Both renames landed after dated files in ``results/`` had
        already quoted commands as their reproduce instructions, so one era of
        that record says `vivum ...` and the other says `codon ...`. Those files
        are not edited, so dropping EITHER entry point would break one era of
        them: both names are registered and both resolve to the same ``main``.
        """
        import re
        root = pathlib.Path(__file__).resolve().parents[1]
        text = (root / "pyproject.toml").read_text(encoding="utf-8")
        block = text.split("[project.scripts]", 1)[1].split("\n[", 1)[0]
        entries = dict(re.findall(r'^(\w+)\s*=\s*"([^"]+)"', block, re.M))
        assert entries.get("vivum") == "tmesim.cli:main"
        assert entries.get("codon") == "tmesim.cli:main"

    def test_compare_drugs_is_refused_and_the_help_points_at_simulate(self, capsys):
        """Converted from ``test_compare_drugs_defaults``.

        That test pinned the defaults of a command that ranked compounds
        against each other (``--library gbm``, the whole-panel ``--drugs``
        default, the ranked ``--csv`` export). The ordering was measured and
        does not hold (``results/pk_ranking_robustness.md``,
        ``results/capacity_ranking_robustness.md``), so the command is gone and
        each of those defaults becomes an assertion that the switch no longer
        exists anywhere in the parser. The replacement is ``simulate`` — one
        tissue, one compound — and the parser has to say so.
        """
        parser = cli.build_parser()

        # was: parse_args(["compare-drugs"]) succeeds.
        # now: argparse refuses it, exit code 2, not a traceback.
        with pytest.raises(SystemExit) as exc:
            parser.parse_args(["compare-drugs"])
        assert exc.value.code == 2
        err = capsys.readouterr().err
        assert "invalid choice: 'compare-drugs'" in err
        # ...and the error is ACTIONABLE: it names the command that replaced it.
        assert "simulate" in err

        # was: args.library == "gbm"       -> no --library on any subcommand
        # was: args.drugs is None          -> no --drugs on any subcommand
        # was: args.csv is None            -> no ranked --csv on any subcommand
        for dead_flag in ("--library", "--drugs", "--csv"):
            with pytest.raises(SystemExit) as exc:
                parser.parse_args(["simulate", dead_flag, "x"])
            assert exc.value.code == 2
            assert dead_flag in capsys.readouterr().err

        # was: args.func is cli.cmd_compare_drugs
        # now: the handler is gone from the module, and `simulate` is the
        # single-compound command the help points at.
        assert not hasattr(cli, "cmd_compare_drugs")
        assert parser.parse_args(["simulate"]).func is cli.cmd_simulate

    def test_studio_defaults(self):
        parser = cli.build_parser()
        args = parser.parse_args(["studio"])
        assert args.host == "127.0.0.1"
        assert args.port == 8000
        assert args.reload is False
        assert args.func is cli.cmd_studio


# --------------------------------------------------------------------------- #
# The single-compound MAP, end-to-end on the real (cached) Xenium GBM section,
# and the removed `compare-drugs` staying removed.
#
# Every test in this class is a CONVERSION of a test that pinned
# `compare-drugs` behaviour. The old contract was "these compounds, ranked".
# The new contract is "this tissue, this compound, and no ranking anywhere" —
# because the cross-compound ordering was measured and failed
# (results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md).
# Nothing that was asserted has been dropped; each assertion was rewritten to
# pin what is true now.
# --------------------------------------------------------------------------- #
class TestOneCompoundMapAndTheRemovedComparison:
    def test_one_real_drug_prints_a_map_and_no_cross_drug_table(self, one_map):
        """Converted from ``test_two_real_drugs_prints_a_sorted_table``.

        was: two drugs on one section print a table sorted by kill fraction,
             with both drug names in it.
        now: ONE drug on one section prints its own map; the second drug is
             simply not part of the run, and the sorted comparison table does
             not exist on any surface.
        """
        rc, out, _out_dir = one_map
        assert rc == 0
        # was: "drug comparison" in out   (the ranked table header)
        assert "drug comparison" not in out
        assert "sorted by predicted kill_fraction" not in out
        # was: "temozolomide" in out      (still true — it is THE compound)
        assert "temozolomide" in out
        assert "=== simulate results ===" in out
        # was: "cisplatin" in out         (a second compound in the same output)
        assert "cisplatin" not in out

    def test_no_ranked_csv_is_written_by_anything(self, one_map, tmp_path,
                                                  capsys):
        """Converted from ``test_writes_a_csv_with_one_row_per_drug``.

        was: a CSV with one row per drug, ordered, is the machine-readable
             export.
        now: that export IS the ranking, so it is gone — no ``--csv`` flag, no
             writer in the module, no file. The two value assertions it made
             (kill fraction in [0,1], a positive IC50) are kept: they are
             properties of a single run and are still readable off the map.
        """
        rc, out, out_dir = one_map
        assert rc == 0

        # was: csv_path.exists() after `--csv`
        with pytest.raises(SystemExit) as exc:
            cli.main(["simulate", "--csv", str(tmp_path / "comparison.csv")])
        assert exc.value.code == 2
        assert "--csv" in capsys.readouterr().err
        assert list(tmp_path.iterdir()) == []
        # was: the ranked writer + its provenance builder are module functions
        assert not hasattr(cli, "_comparison_provenance")
        assert not hasattr(cli, "_provenance_comment_lines")
        # and the run itself wrote no tabular export, only figures
        assert not any(p.suffix == ".csv" for p in out_dir.iterdir())

        # was: {r["drug"] for r in rows} == {"temozolomide", "cisplatin"}
        # now: exactly one compound is named as the drug of the run
        drug_lines = [l for l in out.splitlines()
                      if l.strip().startswith("drug (molecule)")]
        assert len(drug_lines) == 1
        assert "temozolomide" in drug_lines[0]

        # was: 0.0 <= float(r["kill_fraction"]) <= 1.0
        kill = float(next(l for l in out.splitlines()
                          if l.strip().startswith("kill_fraction")
                          ).split(":")[1])
        assert 0.0 <= kill <= 1.0
        # was: float(r["ic50_uM"]) > 0.0
        ic50 = float(next(l for l in out.splitlines()
                          if l.strip().startswith("ic50 / hill")
                          ).split(":")[1].split("uM")[0])
        assert ic50 > 0.0

    def test_the_provenance_the_csv_carried_did_not_evaporate_with_it(
        self, one_map
    ):
        """Converted from ``test_the_exported_csv_carries_its_own_provenance``.

        The CSV was the only machine-readable export and it carried the run's
        whole provenance inside itself. Deleting the export must not delete the
        disclosures — so every key that block asserted is re-asserted here
        against the surface that owns it now: the run's own printed report for
        the section conditions, and the single upstream implementation for the
        calibration constant and the delivery law.
        """
        rc, out, out_dir = one_map
        assert rc == 0

        # (a) the run states its own conditions, in its own output
        # was: "section_source:" in the CSV header
        assert "source                 xenium" in out
        # was: "dx_um:" in the CSV header
        assert f"dx_um                  {_FAST_DX:.4f}" in out
        # was: "exposure_window_s:" in the CSV header
        assert "total_time=1800s" in out
        # was: "sensitivity_mode:" in the CSV header
        assert "sensitivity" in out
        assert "upper_bound" in out

        # (b) the two values the CSV restated from elsewhere are still pinned
        #     at the ONE implementation each of them has, so neither can drift
        # was: "alpha:" in the header and side["alpha"] > 0
        from tmesim.validation.calibration import calibrate_alpha
        assert float(calibrate_alpha()) > 0.0
        # was: "S = PS*" in the header and side["delivery_law"].startswith(...)
        pytest.importorskip("app.server.delivery_law",
                            reason="needs the repo root importable")
        from app.server.delivery_law import delivery_law
        from tmesim.core import PDEConfig
        law = delivery_law(bool(getattr(PDEConfig(), "transvascular_exchange",
                                        False)))
        assert law["law"].startswith("S = PS*")
        # was: "seed:" in the header — still a global, still defaulting to 0
        assert cli.build_parser().parse_args(["simulate"]).seed == 0

        # (c) "NOT FOR CLINICAL USE" travelled inside the artefact. There is no
        #     artefact now, so the equivalent disclosure is pinned where a
        #     reader of the numbers meets them.
        qs = QUICKSTART.read_text(encoding="utf-8")
        assert "not a validated clinical" in " ".join(qs.split())

        # (d) was: "penetration_index is NOT in this file". That retirement was
        #     about CROSS-DRUG comparison — the index is scale-invariant, so it
        #     reported the section's vascular geometry, which every drug shares.
        #     Within ONE map that is exactly what it is for, and `simulate`
        #     prints it.
        assert "penetration_index" in out

        # (e) was: a `.provenance.json` sidecar next to the CSV.
        assert not any(p.name.endswith(".provenance.json")
                       for p in out_dir.iterdir())
        # was: [r["drug"] for r in rows] == ["temozolomide"]
        assert out.count("drug (molecule)") == 1

    def test_coverage_and_the_retired_index_keep_their_places(self, one_map):
        """Converted from ``test_the_csv_carries_coverage_and_not_the_retired_index``.

        `penetration_index` is retired AS A CROSS-DRUG NUMBER
        (results/delivery_metrics.md §4): it is a scale-invariant shape
        statistic that reports the section's vascular geometry, which every
        drug shares. That is a fatal objection to a ranking column and no
        objection at all to a single map, which is what this command now
        produces — so it is printed here, and `therapeutic_coverage` keeps its
        implementation for callers that want it. What is gone is the ranked
        export both of them used to be columns of.
        """
        rc, out, _out_dir = one_map
        assert rc == 0
        # was: "therapeutic_coverage_1x" in rows[0]
        from tmesim.validation.delivery_metrics import therapeutic_coverage
        assert callable(therapeutic_coverage)
        # was: "penetration_index" not in rows[0]  (it was a ranking column)
        # now: it IS the within-map statistic, printed with a value in [0, 1]
        assert "penetration_index" in out
        # was: 0.0 <= float(rows[0]["therapeutic_coverage_1x"]) <= 1.0
        pen = float(next(l for l in out.splitlines()
                         if l.strip().startswith("penetration_index")
                         ).split(":")[1])
        assert 0.0 <= pen <= 1.0
        # and the helper that built the ranked column is gone with the command
        assert not hasattr(cli, "_coverage_at_1x_ic50")

    def test_the_window_is_stated_and_the_law_is_still_single_sourced(
        self, one_map
    ):
        """Converted from ``test_it_prints_the_window_and_the_delivery_law_it_ran_under``.

        The old assertion existed because a RANKING is law-dependent and
        window-dependent and the command printed neither. There is no ranking
        left to be law-dependent, but a map is still window-dependent, so the
        window assertion is kept against `simulate`'s output; and the law
        disclosure is kept against the one implementation that reads PDEConfig
        at call time, so it still cannot go stale.
        """
        rc, out, _out_dir = one_map
        assert rc == 0
        # was: "EXPOSURE WINDOW" in out and "1800 s" in out
        assert "total_time=1800s" in out
        assert cli._FAST_TOTAL_TIME == 1800.0
        # was: "43200 s" in out — the library default, named against this one
        from tmesim.core import PDEConfig
        assert PDEConfig().total_time == 43200.0
        # was: "Delivery law:" / "S = PS*" / "LAW" in out
        pytest.importorskip("app.server.delivery_law",
                            reason="needs the repo root importable")
        from app.server.delivery_law import delivery_law
        for exchange in (True, False):
            assert delivery_law(exchange)["law"].startswith("S = PS*")
            assert delivery_law(exchange)["caveat"]
        # was: the printed table's header carries cov@IC50 and not pen_idx.
        # now: there is no ranked table header at all.
        assert not any("IC50(uM)" in l for l in out.splitlines())
        assert "cov@IC50" not in out
        assert "pen_idx" not in out

    def test_an_unknown_drug_is_now_fatal_rather_than_skipped(self, capsys):
        """Converted from ``test_unknown_drug_is_skipped_not_fatal``.

        was: one bad name in a multi-drug list is SKIPPED with a printed
             reason, so it cannot sink the whole comparison.
        now: there is no list. A run is one compound, so an unknown compound
             leaves nothing to compute and the run fails instead of quietly
             producing a partial result.
        """
        _skip_if_no_real_data()
        rc = cli.main(["simulate", "--dx", str(_FAST_DX),
                       "--drug", "not-a-real-drug"])
        # was: rc == 0 with a SKIPPED line
        assert rc == 1
        captured = capsys.readouterr()
        assert "SKIPPED" not in captured.out
        # was: "temozolomide" in out — the surviving drug of the pair. There is
        # no surviving drug now; the failure names the one that was asked for.
        assert "not-a-real-drug" in captured.out + captured.err

    def test_a_bad_drug_name_is_still_a_clean_failure(self, capsys):
        """Converted from ``test_all_requested_drugs_unknown_is_a_clean_failure``:
        the exit code contract is unchanged, only the command that carries it."""
        _skip_if_no_real_data()
        rc = cli.main(["simulate", "--dx", str(_FAST_DX),
                       "--drug", "not-a-real-drug"])
        assert rc == 1

    def test_the_breast_library_switch_is_gone_and_the_path_is_documented(
        self, capsys
    ):
        """Converted from
        ``test_missing_breast_sample_default_fails_before_touching_disk``.

        was: `--library breast` with no `--sample` resolves the breast cache
             path itself and refuses up front, naming the path and QUICKSTART,
             rather than letting a deep FileNotFoundError surface.
        now: `--library` was part of `compare-drugs` and went with it. Nothing
             resolves a non-default section implicitly any more, so the
             equivalent guarantee is that the user is TOLD the path to pass:
             the constant still names it and QUICKSTART still prints it.
        """
        # was: cli.main(["compare-drugs", "--library", "breast"]) -> rc 1
        with pytest.raises(SystemExit) as exc:
            cli.main(["simulate", "--library", "breast"])
        assert exc.value.code == 2
        assert "--library" in capsys.readouterr().err

        # was: "not present" in out — the message that named the missing cache
        assert cli._BREAST_SAMPLE_PATH == os.path.join(
            "data_cache", "xenium_breast_janesick", "rep1")
        # was: "docs/QUICKSTART.md" in out — the document that explains it
        qs = QUICKSTART.read_text(encoding="utf-8")
        assert "data_cache/xenium_breast_janesick/rep1" in qs
        assert "--sample" in qs

    def test_explicit_bad_sample_fails_with_a_clear_message(self, capsys):
        # An explicit --sample is trusted and passed straight through; if it is
        # wrong, the underlying loader's error should name the path it looked
        # in rather than an unrelated default. Converted only in which command
        # carries it: `simulate` instead of the removed `compare-drugs`.
        rc = cli.main([
            "simulate", "--source", "xenium",
            "--sample", "data_cache/definitely_not_a_real_cache_dir",
        ])
        assert rc == 1
        captured = capsys.readouterr()
        assert ("data_cache/definitely_not_a_real_cache_dir"
                in captured.err + captured.out)


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
# The failure these lock down: `vivum acquire` used to fetch four TABULAR
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
        """The point of `vivum sections`: a remedy, not just a diagnosis."""
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
                assert f"vivum acquire --sections {sid}" in out
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
            assert "curl" in remedy or "vivum ingest" in remedy

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
        assert "vivum sections" in out["note"]
        assert not list(tmp_path.iterdir())

    def test_human_bytes_never_invents_a_number(self):
        from tmesim.data import section_catalog as SC

        assert SC.human_bytes(None) == "unknown size"
        assert SC.human_bytes(0) == "unknown size"
        assert SC.human_bytes(16_598_845) == "16.6 MB"

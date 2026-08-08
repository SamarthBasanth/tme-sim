"""Tests for ``tmesim.cli`` — the ``vivum`` console script / ``python -m
tmesim.cli`` entry point.

Two kinds of coverage:

1. **Argparse plumbing** — every subcommand is registered, ``--help`` works
   for each without crashing, and ``vivum`` is what the tool calls itself now
   that it is a real installed command (``pyproject.toml``'s
   ``[project.scripts]``).
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
import os

import pytest

import tmesim.cli as cli

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
XENIUM_GBM_CACHE = os.path.join(_REPO, "data_cache", "xenium_gbm")

# Coarse grid: matches conftest.py's real_micro fixture, so reconstruction
# finishes in a fraction of a second and the whole test file stays fast.
_FAST_DX = 80.0


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
            "reconstruct", "simulate", "compare-drugs", "validate", "molecules",
            "molecular-run", "acquire", "studio",
        }

    @pytest.mark.parametrize("args", [
        [],
        ["reconstruct", "--help"],
        ["simulate", "--help"],
        ["compare-drugs", "--help"],
        ["validate", "--help"],
        ["molecules", "--help"],
        ["molecular-run", "--help"],
        ["acquire", "--help"],
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
        with open(csv_path, newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert {r["drug"] for r in rows} == {"temozolomide", "cisplatin"}
        for r in rows:
            assert 0.0 <= float(r["kill_fraction"]) <= 1.0
            assert float(r["ic50_uM"]) > 0.0

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

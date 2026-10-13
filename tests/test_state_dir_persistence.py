"""TMESIM_STATE_DIR: the setting that decides whether a hosted workspace keeps
anything, and the three files that have to agree for it to work.

WHY THIS FILE EXISTS. The variable was added on 2026-09-20 so that `runs/`,
`projects/` and `reports/` become children of one directory and a single
platform disk can cover all three. Measured the same day in a built container,
it did not work: a disk mounted at a path the image does not contain is created
ROOT-OWNED, the image runs as uid 10001, and the first use of the value is a
mkdir. The container exited 1 with a PermissionError out of os.makedirs while
printing its startup banner. The whole service, dead at boot, because of the
one setting whose purpose is to stop a redeploy destroying people's work.

Nothing in the suite could have caught it, because the failure needs three
things to be true at once that live in three files: the mount path in
render.yaml, the chowned directory in the Dockerfile, and the code that writes
there. Any one of them changing alone re-breaks it. That is what is pinned
here.
"""
from __future__ import annotations

import io
import os
import re
import stat
import subprocess
import time

import pytest

from app.server import cache

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(rel: str) -> str:
    with io.open(os.path.join(_ROOT, rel), encoding="utf-8") as fh:
        return fh.read()


class TestTheRefusalIsLegible:
    """A traceback is not a diagnosis. The sentence has to carry the fact."""

    def test_unset_is_not_a_problem(self, monkeypatch):
        monkeypatch.delenv("TMESIM_STATE_DIR", raising=False)
        assert cache.state_dir_problem() is None

    def test_a_writable_directory_is_not_a_problem(self, monkeypatch, tmp_path):
        monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path))
        assert cache.state_dir_problem() is None

    def test_the_probe_cleans_up_after_itself(self, monkeypatch, tmp_path):
        """A check that leaves litter in the operator's store is a bug of its
        own -- this directory is about to hold their projects."""
        monkeypatch.setenv("TMESIM_STATE_DIR", str(tmp_path))
        cache.state_dir_problem()
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.skipif(os.getuid() == 0, reason="root ignores the mode bits")
    def test_an_unwritable_directory_is_reported_not_raised(
            self, monkeypatch, tmp_path):
        """THE MEASURED OUTAGE, reproduced without a container.

        The production shape is a root-owned mount; the reproducible shape is
        a directory this process cannot write. Both reach os.makedirs with no
        permission, which is the fault that killed the boot.
        """
        d = tmp_path / "mount"
        d.mkdir()
        os.chmod(d, stat.S_IRUSR | stat.S_IXUSR)  # r-x------: no write
        monkeypatch.setenv("TMESIM_STATE_DIR", str(d))
        try:
            msg = cache.state_dir_problem()
        finally:
            os.chmod(d, stat.S_IRWXU)

        assert msg is not None, (
            "an unwritable state directory reported no problem; the server "
            "will go on to mkdir there and die in a traceback instead")
        assert str(d) in msg, msg
        # The one fact a traceback omits and the operator needs.
        assert f"uid {os.getuid()}" in msg, msg
        assert "owned by" in msg, msg
        # And what to actually do about it.
        assert "/var/lib/vivum" in msg, msg

    @pytest.mark.skipif(os.getuid() == 0, reason="root ignores the mode bits")
    def test_the_refusal_says_why_it_is_not_falling_back(
            self, monkeypatch, tmp_path):
        """Falling back to the ephemeral default would look like working
        persistence and destroy the work anyway. The message has to say that,
        or the next person to read it will 'fix' it by falling back."""
        d = tmp_path / "mount"
        d.mkdir()
        os.chmod(d, stat.S_IRUSR | stat.S_IXUSR)
        monkeypatch.setenv("TMESIM_STATE_DIR", str(d))
        try:
            msg = cache.state_dir_problem()
        finally:
            os.chmod(d, stat.S_IRWXU)
        assert msg is not None
        assert "redeploy" in msg.lower(), msg


class TestTheThreeFilesAgree:
    """The path is a contract between a YAML file, a Dockerfile and a module.
    None of them can see the other two at runtime."""

    def test_the_blueprint_mounts_where_the_image_chowns(self):
        blueprint = _read("render.yaml")
        dockerfile = _read("Dockerfile")

        m = re.search(r"^\s*mountPath:\s*(\S+)\s*$", blueprint, re.M)
        assert m, "render.yaml declares no disk mountPath"
        mount = m.group(1)

        chowned = re.findall(
            r"RUN mkdir -p (/\S+) && chown vivum:vivum (/\S+)", dockerfile)
        assert chowned, (
            "the Dockerfile pre-creates no owned mount point; a volume "
            "mounted on a path the image lacks is root-owned and uid 10001 "
            "cannot write to it")
        made, owned = chowned[0]
        assert made == owned, (made, owned)
        assert mount == made, (
            f"render.yaml mounts the disk at {mount} but the image only "
            f"pre-creates and chowns {made}. A volume at {mount} will be "
            "root-owned, the server will not be able to create its stores "
            "there, and it will refuse to start.")

    def test_the_blueprint_sets_the_variable_that_uses_the_disk(self):
        blueprint = _read("render.yaml")
        m = re.search(r"^\s*mountPath:\s*(\S+)\s*$", blueprint, re.M)
        mount = m.group(1)
        assert re.search(
            r"key:\s*TMESIM_STATE_DIR\s*\n\s*value:\s*" + re.escape(mount),
            blueprint), (
            "the blueprint provisions a disk and never points "
            "TMESIM_STATE_DIR at it, so the server writes to the container "
            "filesystem and the disk holds nothing while appearing to work")

    def test_the_error_message_names_the_path_the_blueprint_uses(self):
        """If the mount path ever moves, the advice in the refusal becomes
        wrong -- and the refusal is the only thing the operator will see."""
        blueprint = _read("render.yaml")
        mount = re.search(r"^\s*mountPath:\s*(\S+)\s*$", blueprint, re.M).group(1)
        src = _read("app/server/cache.py")
        body = src[src.index("def state_dir_problem"):src.index("def runs_dir")]
        assert mount in body, (
            f"state_dir_problem tells the operator to use a path other than "
            f"{mount}, which is where render.yaml actually mounts the disk")


class TestTheServerRefusesRatherThanCrashes:
    def test_boot_exits_with_a_config_code_and_a_sentence(self, tmp_path):
        """End to end, in a real process: an unwritable state dir must produce
        EX_CONFIG and one readable line, not a traceback."""
        if os.getuid() == 0:
            pytest.skip("root ignores the mode bits")
        d = tmp_path / "mount"
        d.mkdir()
        os.chmod(d, stat.S_IRUSR | stat.S_IXUSR)
        env = dict(os.environ)
        env["TMESIM_STATE_DIR"] = str(d)
        env["TMESIM_STUDIO_API_KEY"] = "k" * 24
        try:
            proc = subprocess.run(
                [os.path.join(_ROOT, ".venv", "bin", "python"),
                 "-m", "app.server", "--host", "127.0.0.1", "--port", "8123"],
                cwd=_ROOT, env=env, capture_output=True, text=True, timeout=90)
        finally:
            os.chmod(d, stat.S_IRWXU)

        assert proc.returncode == 78, (
            f"expected EX_CONFIG (78), got {proc.returncode}\n"
            f"stdout: {proc.stdout[-2000:]}\nstderr: {proc.stderr[-2000:]}")
        assert "REFUSING TO START" in proc.stderr, proc.stderr[-2000:]
        assert "Traceback" not in proc.stderr, (
            "the operator got a traceback instead of a diagnosis:\n"
            + proc.stderr[-2000:])


class TestTheServerAcceptsItsOwnDefaults:
    """`GET /api/defaults` must return a spec `POST /api/runs` will take.

    MEASURED ON THE LIVE DEPLOYMENT 2026-09-20. It did not. `/api/defaults`
    served `dx_um = 12.0`; the same server answered that spec with 400,
    because `TMESIM_MAX_GRID_CELLS=auto` resolves to 100,000 on a 512 MB host
    and 12 um on `xenium_319` is ~284,498 voxels. The first click a hosted
    visitor makes failed.

    Nothing coupled the two: `default_spec` was a literal, and the ceiling is
    read from the environment at import. This is the coupling.
    """

    @staticmethod
    def _reloaded(ceiling: str):
        import importlib
        os.environ["TMESIM_MAX_GRID_CELLS"] = ceiling
        import app.server.api as api
        return importlib.reload(api)

    @pytest.mark.parametrize("ceiling", ["100000", "250000", "1000000",
                                         "40000000"])
    def test_the_defaults_validate_under_every_ceiling(self, ceiling):
        prior = os.environ.get("TMESIM_MAX_GRID_CELLS")
        try:
            api = self._reloaded(ceiling)
            spec = api.default_spec()
            verdict = api.validate_spec(spec)
            assert verdict["ok"] is True, (
                f"with TMESIM_MAX_GRID_CELLS={ceiling} the server serves "
                f"defaults it then refuses (dx_um={spec['dx_um']}): "
                f"{verdict['errors']}")
        finally:
            if prior is None:
                os.environ.pop("TMESIM_MAX_GRID_CELLS", None)
            else:
                os.environ["TMESIM_MAX_GRID_CELLS"] = prior
            self._reloaded(prior or "40000000")

    def test_an_unconstrained_server_still_offers_12um(self):
        """No existing installation moves. The coarsening is for hosts that
        would otherwise refuse, and nowhere else."""
        prior = os.environ.get("TMESIM_MAX_GRID_CELLS")
        try:
            api = self._reloaded("40000000")
            assert api.default_spec()["dx_um"] == 12.0
        finally:
            if prior is None:
                os.environ.pop("TMESIM_MAX_GRID_CELLS", None)
            else:
                os.environ["TMESIM_MAX_GRID_CELLS"] = prior
            self._reloaded(prior or "40000000")


class TestTheDrugCatalogueIsCompiledOnce:
    """`GET /api/drugs` did not answer inside 120 s on the live deployment.

    MEASURED 2026-09-20, repeatedly, while `/healthz` answered the SAME
    container in 0.2 s and the same image answered in 0.012 s on a laptop.

    The memo made that look impossible, but the memo is consulted only at
    ENTRY: before the first result exists, every concurrent caller ran the
    whole compile independently. The startup warm thread is one such caller,
    and the Studio's boot does `Promise.all([getSections, getDrugs,
    getDefaults])`, so the whole workstation waits on it. On a fraction of a
    CPU the retries pile up faster than the work completes.
    """

    def test_concurrent_callers_compile_the_library_once(self):
        import threading as _t
        import app.server.api as api

        calls = {"n": 0}
        real = api._compile_drug_rows

        def counting():
            calls["n"] += 1
            time.sleep(0.25)   # make the window wide enough to collide in
            return real()

        prior = api._drug_rows_cache
        api._drug_rows_cache = None
        api._compile_drug_rows = counting
        try:
            out = []
            threads = [_t.Thread(target=lambda: out.append(api.drug_rows()))
                       for _ in range(8)]
            for th in threads:
                th.start()
            for th in threads:
                th.join()
        finally:
            api._compile_drug_rows = real
            api._drug_rows_cache = prior

        assert calls["n"] == 1, (
            f"8 concurrent callers ran the compile {calls['n']} times. Each "
            "one is seconds of work; on a throttled host they starve each "
            "other and no caller finishes.")
        assert len(out) == 8
        assert all(r == out[0] for r in out), "callers disagreed about the rows"

    def test_the_memo_still_short_circuits_without_touching_the_lock(self):
        """The warm path must not pay for the fix. A lock acquired on every
        read would put contention on the common case to protect the rare one.
        """
        import app.server.api as api
        api.drug_rows()                      # ensure warm
        assert api._drug_rows_cache is not None
        t0 = time.perf_counter()
        for _ in range(50):
            api.drug_rows()
        elapsed = time.perf_counter() - t0
        assert elapsed < 1.0, f"50 warm reads took {elapsed:.3f}s"


class TestTheBodyCapCannotBeLiedTo:
    """The 4 MB cap was enforced on Content-Length, a client-supplied header.

    Its own docstring said so honestly -- "may be absent or a lie" -- and a
    security pass on 2026-09-20 showed the lie is free: the same body that
    earns a 413 with a declared length sailed through under
    `Transfer-Encoding: chunked`, because no header declares a size. Only
    `protocol_api` used `_capped_body`; the other `request.json()` sites were
    uncapped in practice.
    """

    @staticmethod
    def _client():
        os.environ["TMESIM_STUDIO_API_KEY"] = "k" * 24
        from fastapi.testclient import TestClient
        import app.server.api as api
        return TestClient(api.create_app()), api.MAX_REQUEST_BODY_BYTES

    def test_an_undeclared_oversized_body_is_refused(self):
        import json as _json
        client, cap = self._client()
        headers = {"X-API-Key": "k" * 24, "Content-Type": "application/json"}
        big = _json.dumps({"name": "x", "notes": "A" * (cap + 1000)}).encode()

        def chunked():
            for i in range(0, len(big), 65536):
                yield big[i:i + 65536]

        r = client.post("/api/projects", headers=headers, content=chunked())
        assert r.status_code == 413, (
            "a body over the cap was accepted because it declared no length; "
            "the cap is enforced on a header the client controls")
        # The refusal must not echo a size -- the read STOPS at the cap, so any
        # number would be where we stopped, not how big the body was.
        body = r.text
        assert str(len(big)) not in body, body[:200]

    def test_a_declared_oversized_body_is_still_refused(self):
        import json as _json
        client, cap = self._client()
        headers = {"X-API-Key": "k" * 24, "Content-Type": "application/json"}
        big = _json.dumps({"name": "x", "notes": "A" * (cap + 1000)}).encode()
        assert client.post("/api/projects", headers=headers,
                           content=big).status_code == 413

    def test_ordinary_requests_are_untouched(self):
        """The cap must not become the outage. A body-reading middleware that
        breaks normal traffic is worse than the exhaustion it prevents."""
        client, _ = self._client()
        headers = {"X-API-Key": "k" * 24}
        assert client.post("/api/projects", headers=headers,
                           json={"name": "small one"}).status_code == 201
        assert client.get("/api/workspace", headers=headers).status_code == 200
        assert client.get("/healthz").status_code == 200

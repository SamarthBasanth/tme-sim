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

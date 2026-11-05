"""Confirms the package never ships COSMIC data by default.

plotting_A's review flagged that COSMIC's published terms restrict
redistributing their signature catalogues, so the sdist/wheel this project
builds and publishes must never contain any -- regardless of what happens
to be sitting in a developer's local checkout (e.g. test fixtures, or a
--cosmic-dir someone configured for themselves). This builds a real sdist
with the project's own `setup.py sdist` and inspects its actual contents,
rather than just checking source files, so a change to packaging config
that accidentally starts bundling data again is caught here.
"""

import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def built_sdist(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("sdist")
    result = subprocess.run(
        [sys.executable, "setup.py", "sdist", "-d", str(out_dir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"setup.py sdist failed:\n{result.stdout}\n{result.stderr}")
    archives = list(out_dir.glob("*.tar.gz"))
    assert len(archives) == 1, f"expected exactly one sdist, found {archives}"
    with tarfile.open(archives[0]) as tar:
        return tar.getnames()


def test_sdist_does_not_contain_cosmic_catalogue_files(built_sdist):
    cosmic_files = [name for name in built_sdist if "cosmic_signatures" in name or "COSMIC_v" in name]
    assert cosmic_files == [], f"sdist must not bundle COSMIC data, found: {cosmic_files}"


def test_cosmic_signatures_directory_does_not_exist_in_the_repo():
    # The directory itself was removed, not just excluded from packaging --
    # nothing should be able to reintroduce bundled data by accident.
    directory = REPO_ROOT / "sigProfilerPlotting" / "web" / "cosmic_signatures"
    assert not directory.exists()


def test_setup_py_and_manifest_do_not_reference_cosmic_signatures():
    setup_text = (REPO_ROOT / "setup.py").read_text(encoding="utf-8")
    manifest_text = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "cosmic_signatures" not in setup_text
    assert "cosmic_signatures" not in manifest_text

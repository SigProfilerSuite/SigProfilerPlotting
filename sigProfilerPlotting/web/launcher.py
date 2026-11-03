"""Idempotent venv management and process lifecycle for the local dashboard.

This module is deliberately split into two entry points because it is reached
two different ways:

* The double-click launcher scripts (``launch_dashboard.sh`` /
  ``.command`` / ``.bat``) manage a private, repo-local virtualenv and must
  install/upgrade the package into it -- but only when something relevant
  changed, never on every launch. That is :func:`managed_main`, reached via
  ``python -m sigProfilerPlotting.web.launcher`` once the venv exists, or by
  executing this file directly by path (``--ensure-deps-only``/``--stop``
  only use the standard library, so that works even *before* the package has
  been installed into the fresh venv -- see the bootstrap step those scripts
  run first).
* The ``sigprofilerplotting-web`` console script, installed by ``pip`` into
  whatever environment the user chose, where dependency resolution already
  happened at install time. That is :func:`main` -- it never touches a venv
  or runs ``pip`` itself.

Both share the same argument parser, the same pidfile-based ``--stop``, and
the same server start-up, so the two entry points only disagree about
whether a dependency-install step runs first.

All user-facing status lines use ``flush=True`` (or write straight to a
subprocess's own stdout). Launch scripts and CI jobs often capture this
process's output through a pipe/file rather than a TTY, where Python's
default stdout buffering is fully block-buffered rather than line-buffered --
without an explicit flush, "already installed"/"dashboard running at ..."
can sit unflushed for a long-running `serve_forever()` process and never
reach a log a caller is tailing or grepping.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import signal
import subprocess
import sys
import webbrowser
from pathlib import Path

_PACKAGING_FILES = ("setup.py", "pyproject.toml")


def _repo_root() -> Path:
    # sigProfilerPlotting/web/launcher.py -> sigProfilerPlotting/web -> sigProfilerPlotting -> repo root
    return Path(__file__).resolve().parents[2]


def compute_fingerprint(repo_root: Path) -> str:
    """Hash everything that determines what needs to be installed.

    Hashing the packaging files themselves (rather than trying to parse and
    hash just ``install_requires``) catches any packaging-relevant change --
    a new/updated dependency, a version bump, a changed entry point -- without
    fragile parsing or a second, easily-drifting copy of the dependency list.
    """
    hasher = hashlib.sha256()
    for name in _PACKAGING_FILES:
        candidate = repo_root / name
        if candidate.is_file():
            hasher.update(candidate.read_bytes())
    hasher.update(f"{sys.version_info.major}.{sys.version_info.minor}".encode("ascii"))
    return hasher.hexdigest()


def needs_install(repo_root: Path, fingerprint_path: Path, force: bool) -> bool:
    """Pure decision function: should an install run? (unit-testable without a real venv)"""
    if force:
        return True
    if not fingerprint_path.is_file():
        return True
    try:
        stored = fingerprint_path.read_text(encoding="utf-8").strip()
    except OSError:
        return True
    return stored != compute_fingerprint(repo_root)


def _fingerprint_path() -> Path:
    return Path(sys.prefix) / ".fingerprint"


def ensure_dependencies(repo_root: Path, force: bool = False) -> None:
    fingerprint_path = _fingerprint_path()
    if needs_install(repo_root, fingerprint_path, force):
        print("Installing dashboard dependencies...", flush=True)
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "--quiet", "-e", str(repo_root)],
            check=True,
        )
        fingerprint_path.write_text(compute_fingerprint(repo_root), encoding="utf-8")
    else:
        print("Dependencies already installed; skipping install.", flush=True)


def _state_dir() -> Path:
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", str(Path.home() / "AppData" / "Local")))
        directory = base / "SigProfilerPlotting"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache")))
        directory = base / "sigprofilerplotting"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _pidfile() -> Path:
    return _state_dir() / "dashboard.pid"


def _read_pidfile(pidfile: Path) -> tuple[int, int] | None:
    """Return (pid, port) from the pidfile, or None if it's absent/unparseable."""
    if not pidfile.is_file():
        return None
    try:
        pid_text, port_text = pidfile.read_text(encoding="utf-8").strip().split(":", 1)
        return int(pid_text), int(port_text)
    except (OSError, ValueError):
        return None


def _is_process_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, just not ours to signal with 0 -- still alive
    return True


def _running_dashboard_url(pidfile: Path) -> str | None:
    """URL of an already-running dashboard, or None if none is live.

    As a side effect, removes the pidfile if it names a process that is no
    longer running (so the caller can safely claim it next).
    """
    existing = _read_pidfile(pidfile)
    if existing is None:
        return None
    pid, port = existing
    if _is_process_alive(pid):
        return f"http://127.0.0.1:{port}/"
    pidfile.unlink(missing_ok=True)
    return None


def stop_server() -> None:
    pidfile = _pidfile()
    existing = _read_pidfile(pidfile)
    if existing is None:
        pidfile.unlink(missing_ok=True)
        print("No running dashboard server found.", flush=True)
        return
    pid, _port = existing
    try:
        os.kill(pid, signal.SIGTERM)
        print(f"Stopped dashboard server (pid {pid}).", flush=True)
    except ProcessLookupError:
        print("No running dashboard server found.", flush=True)
    except PermissionError:
        print(f"Could not stop process {pid}: permission denied.", flush=True)
    finally:
        # The stopper (not the stopped process) owns pidfile cleanup: on
        # Windows, os.kill(..., SIGTERM) maps to TerminateProcess, which kills
        # the target immediately without running its own cleanup code.
        pidfile.unlink(missing_ok=True)


def _run_server(no_launch: bool) -> None:
    pidfile = _pidfile()
    existing_url = _running_dashboard_url(pidfile)
    if existing_url is not None:
        # Judgment call: refuse rather than stop-and-restart, even under
        # --reinstall -- a double-click/double-run should never silently
        # kill a dashboard the user (or another tab) may be using.
        print(f"SigProfilerPlotting dashboard is already running at {existing_url}", flush=True)
        print("Run with --stop first if you want to restart it.", flush=True)
        return

    from .server import make_server  # deferred: requires the package to be importable

    server = make_server("127.0.0.1", 0)
    port = server.server_address[1]
    pidfile.write_text(f"{os.getpid()}:{port}", encoding="utf-8")

    def _handle_term(signum, frame):
        raise SystemExit(0)

    if hasattr(signal, "SIGTERM"):
        signal.signal(signal.SIGTERM, _handle_term)

    url = f"http://127.0.0.1:{port}/"
    print(f"SigProfilerPlotting dashboard running at {url}", flush=True)
    print("Press Ctrl+C to stop, or run with --stop from another terminal.", flush=True)
    if not no_launch:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except (KeyboardInterrupt, SystemExit):
        pass
    finally:
        server.server_close()
        pidfile.unlink(missing_ok=True)


def _build_parser(prog: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog, description="Launch the local SigProfilerPlotting interactive dashboard."
    )
    parser.add_argument(
        "--no-launch", action="store_true", help="Start the server without opening a browser tab."
    )
    parser.add_argument(
        "--reinstall",
        action="store_true",
        help="Force a fresh dependency install even if the fingerprint matches.",
    )
    parser.add_argument(
        "--stop", action="store_true", help="Stop a running dashboard server and exit."
    )
    # Internal: used only by the launcher shell scripts, before the package is
    # necessarily importable yet, to run just the install-fingerprint check.
    parser.add_argument("--ensure-deps-only", action="store_true", help=argparse.SUPPRESS)
    return parser


def managed_main(argv: list[str] | None = None) -> int:
    """Entry point for the venv-managed flow (the launcher shell scripts)."""
    args = _build_parser("sigProfilerPlotting dashboard launcher").parse_args(argv)
    if args.stop:
        stop_server()
        return 0
    if args.ensure_deps_only:
        ensure_dependencies(_repo_root(), force=args.reinstall)
        return 0
    ensure_dependencies(_repo_root(), force=args.reinstall)
    _run_server(no_launch=args.no_launch)
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point for the installed ``sigprofilerplotting-web`` console script.

    Runs directly in whatever environment ``pip`` installed into -- no venv
    creation and no dependency installation, since pip already resolved
    dependencies there.
    """
    args = _build_parser("sigprofilerplotting-web").parse_args(argv)
    if args.stop:
        stop_server()
        return 0
    if args.reinstall:
        print(
            "--reinstall has no effect for the installed console command; "
            "reinstall the package with pip directly, e.g. "
            "`pip install --force-reinstall sigProfilerPlotting`.",
            flush=True,
        )
    _run_server(no_launch=args.no_launch)
    return 0


if __name__ == "__main__":
    sys.exit(managed_main())

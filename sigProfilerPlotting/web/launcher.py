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
import json
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.request
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


def _lockfile() -> Path:
    return _state_dir() / "dashboard.lock"


def _read_pidfile(pidfile: Path) -> tuple[int, int, str] | None:
    """Return (pid, port, token) from the pidfile, or None if absent/unparseable."""
    if not pidfile.is_file():
        return None
    try:
        pid_text, port_text, token = pidfile.read_text(encoding="utf-8").strip().split(":", 2)
        return int(pid_text), int(port_text), token
    except (OSError, ValueError):
        return None


def _write_pidfile_atomic(pidfile: Path, pid: int, port: int, token: str) -> None:
    # Write-then-rename is atomic on POSIX and Windows alike, so a concurrent
    # reader never observes a half-written pidfile.
    tmp = pidfile.with_name(pidfile.name + f".tmp{pid}")
    tmp.write_text(f"{pid}:{port}:{token}", encoding="utf-8")
    os.replace(tmp, pidfile)


def _is_process_alive(pid: int) -> bool:
    """Whether a process with this pid is still running.

    Not os.kill(pid, 0) on Windows: there, signal 0 is CTRL_C_EVENT, so the
    call sends Ctrl+C to every process on this console (the caller included)
    instead of probing the pid.
    """
    if pid <= 0:
        return False
    if sys.platform == "win32":
        return _windows_process_alive(pid)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, just not ours to signal with 0 -- still alive
    return True


def _windows_process_alive(pid: int) -> bool:
    import ctypes
    from ctypes import wintypes

    synchronize, query_limited_information = 0x00100000, 0x1000
    error_access_denied, wait_timeout = 5, 0x102

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.WaitForSingleObject.argtypes = (wintypes.HANDLE, wintypes.DWORD)
    kernel32.WaitForSingleObject.restype = wintypes.DWORD
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    kernel32.CloseHandle.restype = wintypes.BOOL

    handle = kernel32.OpenProcess(synchronize | query_limited_information, False, pid)
    if not handle:
        # Access denied: the process exists but belongs to someone else.
        return ctypes.get_last_error() == error_access_denied
    try:
        # A process that has exited is signalled, even while a handle to it
        # (such as a Popen's) keeps its pid from being reused.
        return kernel32.WaitForSingleObject(handle, 0) == wait_timeout
    finally:
        kernel32.CloseHandle(handle)


def _fetch_identity(port: int, timeout: float = 1.0) -> dict | None:
    """GET /api/identity on a candidate dashboard, or None if unreachable/invalid."""
    try:
        with urllib.request.urlopen(
            f"http://127.0.0.1:{port}/api/identity", timeout=timeout
        ) as response:
            return json.loads(response.read().decode("utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None


def _verify_identity(pid: int, port: int, token: str) -> bool:
    """Confirm the process actually listening on ``port`` is *this* dashboard.

    A pidfile naming a live PID is not proof by itself -- the PID could have
    been reused by an unrelated process, or the port taken over by something
    else entirely. The dashboard's own loopback-only ``/api/identity``
    endpoint echoes back its real PID and a random per-launch token minted at
    startup; both must match what the pidfile recorded before we ever trust
    it enough to reuse or stop it.
    """
    identity = _fetch_identity(port)
    return bool(identity) and identity.get("pid") == pid and identity.get("token") == token


def _running_dashboard_url(pidfile: Path) -> str | None:
    """URL of an already-running dashboard, or None if none is live.

    As a side effect, removes the pidfile if it does not name a dashboard
    that is verifiably still running under that exact pid/port/token triple
    (so the caller can safely claim the pidfile next).
    """
    existing = _read_pidfile(pidfile)
    if existing is None:
        return None
    pid, port, token = existing
    if _verify_identity(pid, port, token):
        return f"http://127.0.0.1:{port}/"
    pidfile.unlink(missing_ok=True)
    return None


def _remove_pidfile_if_matches(pidfile: Path, pid: int, port: int, token: str) -> None:
    """Delete ``pidfile`` only if it still names this exact pid/port/token.

    A bare ``unlink()`` in a cleanup path is unsafe: between this instance
    writing its pidfile and reaching cleanup, a *different* instance (a
    replacement started after this one began shutting down) may have
    overwritten the pidfile with its own record. Comparing before deleting
    means an old instance's cleanup can never erase a newer instance's
    record -- it only ever removes the record it itself wrote.
    """
    if _read_pidfile(pidfile) == (pid, port, token):
        pidfile.unlink(missing_ok=True)


def _wait_for_identity_gone(port: int, timeout: float = 5.0) -> bool:
    """Wait for ``/api/identity`` on ``port`` to stop responding.

    A successful shutdown POST only means the request was accepted -- the
    server may take a moment to unwind ``serve_forever()`` and release the
    port. Returns True once it has actually gone, False if ``timeout``
    elapses first -- callers must treat a timeout as "still possibly
    live," never as "stopped," since those require different cleanup.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if _fetch_identity(port) is None:
            return True
        time.sleep(0.05)
    return _fetch_identity(port) is None


def _acquire_startup_lock(timeout: float = 5.0) -> Path | None:
    """Atomically claim the right to check-and-start a dashboard.

    ``os.open`` with ``O_CREAT | O_EXCL`` is an atomic test-and-set on both
    POSIX and Windows, so two launchers racing to start at the same instant
    converge on exactly one of them doing the check-then-write sequence in
    :func:`_run_server`; the loser waits (or gives up) instead of each
    independently deciding no dashboard exists yet and both binding a server.
    """
    lock_path = _lockfile()
    deadline = time.monotonic() + timeout
    while True:
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(os.getpid()).encode("ascii"))
            os.close(fd)
            return lock_path
        except FileExistsError:
            try:
                holder_pid = int(lock_path.read_text(encoding="utf-8").strip())
            except (OSError, ValueError):
                holder_pid = None
            if holder_pid is None or not _is_process_alive(holder_pid):
                lock_path.unlink(missing_ok=True)
                continue
            if time.monotonic() >= deadline:
                return None
            time.sleep(0.05)


def stop_server() -> bool:
    """Stop the running dashboard, if any. Returns False only on a real failure.

    A missing/stale/unverifiable pidfile is a clean no-op (True) -- there is
    nothing to stop. Only a verified-live dashboard that fails to actually
    shut down is a failure (False); in that case the pidfile is deliberately
    *kept*, since it may still be the only record of a genuinely running
    server -- erasing it here would let a later launch start a second
    instance on top of the first.
    """
    pidfile = _pidfile()
    existing = _read_pidfile(pidfile)
    if existing is None:
        pidfile.unlink(missing_ok=True)
        print("No running dashboard server found.", flush=True)
        return True
    pid, port, token = existing
    if not _verify_identity(pid, port, token):
        # The pidfile does not name a verifiable dashboard -- it is stale
        # (process exited) or the pid/port has been reused by something
        # else. Never signal a PID solely because it appears in a file.
        pidfile.unlink(missing_ok=True)
        print("No running dashboard server found.", flush=True)
        return True
    try:
        request = urllib.request.Request(
            f"http://127.0.0.1:{port}/api/shutdown",
            data=b"{}",
            method="POST",
            headers={
                "Origin": f"http://127.0.0.1:{port}",
                "Content-Type": "application/json",
                "X-Dashboard-Token": token,
            },
        )
        urllib.request.urlopen(request, timeout=5.0)
    except (OSError, urllib.error.HTTPError) as exc:
        print(
            f"Could not stop the dashboard server: {exc} "
            "(it may still be running -- the pidfile has been kept; try again).",
            flush=True,
        )
        return False

    if not _wait_for_identity_gone(port):
        # The shutdown POST was accepted, but the server never actually
        # went away within the timeout -- it may be stalled mid-shutdown or
        # simply still alive. Treat this exactly like a failed shutdown:
        # keep the pidfile, since it may be the only record of a server
        # that is, in fact, still running.
        print(
            f"Dashboard server (pid {pid}) did not confirm shutdown in time "
            "-- it may still be running; the pidfile has been kept.",
            flush=True,
        )
        return False

    # Coordinate with the startup lock so a launch racing this cleanup can
    # never observe a half-stopped server as "nothing running" before this
    # pidfile removal (or its own replacement write) has settled. If the
    # lock itself cannot be acquired, skip the removal entirely rather than
    # compare-and-delete unlocked -- that would race the lock holder's own
    # read of this same pidfile. The next launch or `--stop` will verify
    # and clean it up.
    lock_path = _acquire_startup_lock()
    if lock_path is None:
        print(
            f"Stopped dashboard server (pid {pid}); its pidfile will be "
            "cleaned up by the next launch.",
            flush=True,
        )
        return True
    try:
        _remove_pidfile_if_matches(pidfile, pid, port, token)
    finally:
        lock_path.unlink(missing_ok=True)
    print(f"Stopped dashboard server (pid {pid}).", flush=True)
    return True


def _run_server(no_launch: bool, cosmic_dir: str | None = None) -> None:
    if cosmic_dir:
        from . import cosmic

        cosmic.configure_catalogue_dir(cosmic_dir)

    pidfile = _pidfile()
    lock_path = _acquire_startup_lock()
    if lock_path is None:
        print(
            "Another dashboard launch is already in progress; try again in a moment.",
            flush=True,
        )
        return
    try:
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
        token = server.instance_token
        _write_pidfile_atomic(pidfile, os.getpid(), port, token)
    finally:
        lock_path.unlink(missing_ok=True)

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
        # Compare-and-delete, under the startup lock: a replacement instance
        # may have already written its own pid/port/token into this pidfile
        # (e.g. it started right after this server began shutting down but
        # before this `finally` ran) -- this must never delete that record.
        # If the lock itself cannot be acquired, skip the removal entirely
        # rather than compare-and-delete unlocked, which would race the
        # lock holder's own read of this same pidfile; the next launch or
        # `--stop` will verify and clean it up.
        cleanup_lock = _acquire_startup_lock()
        if cleanup_lock is not None:
            try:
                _remove_pidfile_if_matches(pidfile, os.getpid(), port, token)
            finally:
                cleanup_lock.unlink(missing_ok=True)


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
    parser.add_argument(
        "--cosmic-dir",
        metavar="PATH",
        default=None,
        help=(
            "Folder holding your own authorized COSMIC SBS reference files "
            "(named like COSMIC_v3.4_SBS_GRCh37.txt) to enable the dashboard's "
            "COSMIC comparison panel. SigProfilerPlotting does not bundle or "
            "redistribute this data -- see docs/Plotting-Dashboard.md for where "
            "to obtain it. Defaults to the SIGPROFILERPLOTTING_COSMIC_DIR "
            "environment variable if this flag is omitted."
        ),
    )
    # Internal: used only by the launcher shell scripts, before the package is
    # necessarily importable yet, to run just the install-fingerprint check.
    parser.add_argument("--ensure-deps-only", action="store_true", help=argparse.SUPPRESS)
    return parser


def managed_main(argv: list[str] | None = None) -> int:
    """Entry point for the venv-managed flow (the launcher shell scripts)."""
    args = _build_parser("sigProfilerPlotting dashboard launcher").parse_args(argv)
    if args.stop:
        return 0 if stop_server() else 1
    if args.ensure_deps_only:
        ensure_dependencies(_repo_root(), force=args.reinstall)
        return 0
    ensure_dependencies(_repo_root(), force=args.reinstall)
    _run_server(no_launch=args.no_launch, cosmic_dir=args.cosmic_dir)
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point for the installed ``sigprofilerplotting-web`` console script.

    Runs directly in whatever environment ``pip`` installed into -- no venv
    creation and no dependency installation, since pip already resolved
    dependencies there.
    """
    args = _build_parser("sigprofilerplotting-web").parse_args(argv)
    if args.stop:
        return 0 if stop_server() else 1
    if args.reinstall:
        print(
            "--reinstall has no effect for the installed console command; "
            "reinstall the package with pip directly, e.g. "
            "`pip install --force-reinstall sigProfilerPlotting`.",
            flush=True,
        )
    _run_server(no_launch=args.no_launch, cosmic_dir=args.cosmic_dir)
    return 0


if __name__ == "__main__":
    sys.exit(managed_main())

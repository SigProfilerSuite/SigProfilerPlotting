import contextlib
import io
import os
import subprocess
import sys
import threading
import time

import pytest

from sigProfilerPlotting.web import launcher


def _write_packaging_files(root):
    (root / "setup.py").write_text("setup(name='x')\n", encoding="utf-8")
    (root / "pyproject.toml").write_text("[build-system]\n", encoding="utf-8")


def test_fingerprint_is_stable_for_unchanged_inputs(tmp_path):
    _write_packaging_files(tmp_path)
    assert launcher.compute_fingerprint(tmp_path) == launcher.compute_fingerprint(tmp_path)


def test_fingerprint_changes_when_setup_py_changes(tmp_path):
    _write_packaging_files(tmp_path)
    before = launcher.compute_fingerprint(tmp_path)
    (tmp_path / "setup.py").write_text("setup(name='x', version='2')\n", encoding="utf-8")
    after = launcher.compute_fingerprint(tmp_path)
    assert before != after


def test_needs_install_when_no_fingerprint_file_exists(tmp_path):
    _write_packaging_files(tmp_path)
    fingerprint_path = tmp_path / ".fingerprint"
    assert launcher.needs_install(tmp_path, fingerprint_path, force=False) is True


def test_needs_install_is_false_once_fingerprint_matches(tmp_path):
    _write_packaging_files(tmp_path)
    fingerprint_path = tmp_path / ".fingerprint"
    fingerprint_path.write_text(launcher.compute_fingerprint(tmp_path), encoding="utf-8")
    assert launcher.needs_install(tmp_path, fingerprint_path, force=False) is False


def test_needs_install_true_when_forced_even_if_matching(tmp_path):
    _write_packaging_files(tmp_path)
    fingerprint_path = tmp_path / ".fingerprint"
    fingerprint_path.write_text(launcher.compute_fingerprint(tmp_path), encoding="utf-8")
    assert launcher.needs_install(tmp_path, fingerprint_path, force=True) is True


def test_needs_install_true_when_stored_fingerprint_is_stale(tmp_path):
    _write_packaging_files(tmp_path)
    fingerprint_path = tmp_path / ".fingerprint"
    fingerprint_path.write_text("not-a-real-fingerprint", encoding="utf-8")
    assert launcher.needs_install(tmp_path, fingerprint_path, force=False) is True


def test_ensure_dependencies_skips_pip_when_fingerprint_matches(tmp_path, monkeypatch):
    _write_packaging_files(tmp_path)
    fake_prefix = tmp_path / "venv"
    fake_prefix.mkdir()
    monkeypatch.setattr(launcher.sys, "prefix", str(fake_prefix))
    (fake_prefix / ".fingerprint").write_text(
        launcher.compute_fingerprint(tmp_path), encoding="utf-8"
    )

    calls = []
    monkeypatch.setattr(launcher.subprocess, "run", lambda *a, **k: calls.append((a, k)))
    launcher.ensure_dependencies(tmp_path, force=False)
    assert calls == []


def test_ensure_dependencies_installs_when_fingerprint_is_missing(tmp_path, monkeypatch):
    _write_packaging_files(tmp_path)
    fake_prefix = tmp_path / "venv"
    fake_prefix.mkdir()
    monkeypatch.setattr(launcher.sys, "prefix", str(fake_prefix))

    calls = []

    def fake_run(cmd, check):
        calls.append(cmd)

    monkeypatch.setattr(launcher.subprocess, "run", fake_run)
    launcher.ensure_dependencies(tmp_path, force=False)
    assert len(calls) == 1
    assert (fake_prefix / ".fingerprint").is_file()
    assert (fake_prefix / ".fingerprint").read_text() == launcher.compute_fingerprint(tmp_path)


def test_help_flag_exits_cleanly():
    parser = launcher._build_parser("test")
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--help"])
    assert exc.value.code == 0


def test_flags_parse_as_expected():
    parser = launcher._build_parser("test")
    args = parser.parse_args(["--no-launch", "--reinstall", "--stop"])
    assert args.no_launch is True
    assert args.reinstall is True
    assert args.stop is True
    assert args.ensure_deps_only is False


def test_stop_with_no_pidfile_is_a_clean_noop(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    launcher.stop_server()
    assert "No running dashboard server found." in capsys.readouterr().out


def test_stop_with_a_stale_pidfile_cleans_up_without_error(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    pidfile = tmp_path / "dashboard.pid"
    # A pid astronomically unlikely to be a real running process, and a port
    # nothing is listening on -- /api/identity is unreachable either way.
    pidfile.write_text("999999999:51000:sometoken", encoding="utf-8")
    launcher.stop_server()
    assert not pidfile.exists()
    assert "No running dashboard server found." in capsys.readouterr().out


def test_stop_with_a_corrupt_pidfile_cleans_up_without_error(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("not-a-pid", encoding="utf-8")
    launcher.stop_server()
    assert not pidfile.exists()
    assert "No running dashboard server found." in capsys.readouterr().out


def test_stop_terminates_a_real_dashboard_and_removes_pidfile(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    from sigProfilerPlotting.web.server import make_server

    server = make_server("127.0.0.1", 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    port = server.server_address[1]
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text(f"{os.getpid()}:{port}:{server.instance_token}", encoding="utf-8")
    try:
        launcher.stop_server()
        assert not pidfile.exists()
        assert f"Stopped dashboard server (pid {os.getpid()})" in capsys.readouterr().out
        thread.join(timeout=5)
        assert not thread.is_alive()
    finally:
        if thread.is_alive():
            server.shutdown()
        server.server_close()


def test_stop_does_not_kill_an_unrelated_process_named_by_a_stale_pidfile(
    tmp_path, monkeypatch, capsys
):
    # The exact scenario plotting_A's review flagged: a pidfile naming a
    # real, live, but completely unrelated process (not a dashboard at all --
    # nothing is even listening on the recorded port). stop_server() must
    # never signal a PID solely because a file names it; it must verify via
    # /api/identity first and, finding nothing there, leave the process alone.
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(10)"])
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text(f"{process.pid}:51000:not-the-real-token", encoding="utf-8")
    try:
        launcher.stop_server()
        assert not pidfile.exists()
        assert "No running dashboard server found." in capsys.readouterr().out
        assert process.poll() is None, "an unrelated live process must never be signaled"
    finally:
        process.kill()
        process.wait(timeout=5)


def test_managed_main_stop_short_circuits_before_ensure_dependencies(monkeypatch):
    calls = []
    monkeypatch.setattr(launcher, "stop_server", lambda: calls.append("stop"))
    monkeypatch.setattr(
        launcher, "ensure_dependencies", lambda *a, **k: calls.append("ensure")
    )
    assert launcher.managed_main(["--stop"]) == 0
    assert calls == ["stop"]


# -- Regression: launching twice without --stop in between must not orphan
# the first server (found by the user running `./launch_dashboard.sh
# --no-launch` twice in a row and reproducing it directly: the pidfile was
# silently overwritten by the second launch, leaving the first server
# unreachable via --stop). --------------------------------------------------


def test_is_process_alive_true_for_the_current_process():
    assert launcher._is_process_alive(os.getpid()) is True


def test_is_process_alive_false_for_an_unused_pid():
    assert launcher._is_process_alive(999999999) is False


def test_running_dashboard_url_is_none_without_a_pidfile(tmp_path):
    assert launcher._running_dashboard_url(tmp_path / "dashboard.pid") is None


def test_running_dashboard_url_reports_a_verified_live_dashboard(tmp_path, monkeypatch):
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000:thetoken", encoding="utf-8")
    monkeypatch.setattr(
        launcher, "_verify_identity", lambda pid, port, token: (pid, port, token) == (4242, 51000, "thetoken")
    )
    assert launcher._running_dashboard_url(pidfile) == "http://127.0.0.1:51000/"
    assert pidfile.exists()  # a verified-live dashboard's pidfile is left alone


def test_running_dashboard_url_cleans_up_when_identity_cannot_be_verified(tmp_path, monkeypatch):
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000:thetoken", encoding="utf-8")
    monkeypatch.setattr(launcher, "_verify_identity", lambda pid, port, token: False)
    assert launcher._running_dashboard_url(pidfile) is None
    assert not pidfile.exists()


def test_verify_identity_rejects_a_token_mismatch(tmp_path, monkeypatch):
    from sigProfilerPlotting.web.server import make_server

    server = make_server("127.0.0.1", 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        assert launcher._verify_identity(os.getpid(), port, server.instance_token) is True
        assert launcher._verify_identity(os.getpid(), port, "wrong-token") is False
        assert launcher._verify_identity(999999999, port, server.instance_token) is False
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_run_server_refuses_a_second_instance_without_touching_the_pidfile(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000:thetoken", encoding="utf-8")  # a "live" first instance
    monkeypatch.setattr(launcher, "_pidfile", lambda: pidfile)
    monkeypatch.setattr(launcher, "_running_dashboard_url", lambda p: "http://127.0.0.1:51000/")

    def fail_if_called(*a, **k):
        raise AssertionError("make_server must not run a second instance")

    import sigProfilerPlotting.web.server as server_module

    monkeypatch.setattr(server_module, "make_server", fail_if_called)

    launcher._run_server(no_launch=True)

    out = capsys.readouterr().out
    assert "already running at http://127.0.0.1:51000/" in out
    # The original instance's record must survive untouched -- this is
    # exactly what the reported bug overwrote.
    assert pidfile.read_text(encoding="utf-8") == "4242:51000:thetoken"


def test_acquire_startup_lock_is_exclusive_until_released(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    held = launcher._acquire_startup_lock(timeout=0.2)
    assert held is not None
    assert launcher._acquire_startup_lock(timeout=0.2) is None  # already held
    held.unlink()
    reacquired = launcher._acquire_startup_lock(timeout=0.2)
    assert reacquired is not None
    reacquired.unlink()


def test_acquire_startup_lock_recovers_from_a_stale_lock(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    lock_path = tmp_path / "dashboard.lock"
    lock_path.write_text("999999999", encoding="utf-8")  # holder is not alive
    acquired = launcher._acquire_startup_lock(timeout=1.0)
    assert acquired is not None
    acquired.unlink()


def test_concurrent_launches_leave_exactly_one_managed_server(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    # signal.signal() only works on the process's main thread; _run_server
    # normally only ever runs there (one launcher per process). Running it on
    # a background thread here, just to simulate two racing launches, would
    # otherwise raise ValueError on this call -- stub it out instead.
    monkeypatch.setattr(launcher.signal, "signal", lambda *a, **k: None)
    pidfile = tmp_path / "dashboard.pid"

    first_output = io.StringIO()

    def run_first():
        with contextlib.redirect_stdout(first_output):
            launcher._run_server(no_launch=True)

    thread1 = threading.Thread(target=run_first, daemon=True)
    thread1.start()
    try:
        deadline = time.monotonic() + 5
        while not pidfile.is_file() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert pidfile.is_file(), "the first launch never reached a bound server"

        second_output = io.StringIO()
        with contextlib.redirect_stdout(second_output):
            launcher._run_server(no_launch=True)
        assert "already running at" in second_output.getvalue()
    finally:
        launcher.stop_server()
        thread1.join(timeout=5)
    assert not thread1.is_alive()
    assert "already running at" not in first_output.getvalue()

import os
import subprocess
import sys

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
    # A pid astronomically unlikely to be a real running process.
    pidfile.write_text("999999999:51000", encoding="utf-8")
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


def test_stop_terminates_a_real_process_and_removes_pidfile(tmp_path, monkeypatch, capsys):
    # subprocess.Popen (not os.fork directly) to avoid the known
    # fork()-in-a-multithreaded-process hazard inside a pytest run.
    monkeypatch.setattr(launcher, "_state_dir", lambda: tmp_path)
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(10)"])
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text(f"{process.pid}:51000", encoding="utf-8")
    try:
        launcher.stop_server()
        assert not pidfile.exists()
        assert f"Stopped dashboard server (pid {process.pid})" in capsys.readouterr().out
        assert process.wait(timeout=5) is not None
    finally:
        if process.poll() is None:
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


def test_running_dashboard_url_reports_a_live_process(tmp_path, monkeypatch):
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000", encoding="utf-8")
    monkeypatch.setattr(launcher, "_is_process_alive", lambda pid: pid == 4242)
    assert launcher._running_dashboard_url(pidfile) == "http://127.0.0.1:51000/"
    assert pidfile.exists()  # a live process's pidfile is left alone


def test_running_dashboard_url_cleans_up_a_dead_pid(tmp_path, monkeypatch):
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000", encoding="utf-8")
    monkeypatch.setattr(launcher, "_is_process_alive", lambda pid: False)
    assert launcher._running_dashboard_url(pidfile) is None
    assert not pidfile.exists()


def test_run_server_refuses_a_second_instance_without_touching_the_pidfile(
    tmp_path, monkeypatch, capsys
):
    pidfile = tmp_path / "dashboard.pid"
    pidfile.write_text("4242:51000", encoding="utf-8")  # a "live" first instance
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
    assert pidfile.read_text(encoding="utf-8") == "4242:51000"

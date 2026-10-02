
try:
    from ._version import version as _scm_version
except Exception:
    _scm_version = None

try:
    from importlib.metadata import PackageNotFoundError, version as _dist_version
except Exception:  # pragma: no cover
    PackageNotFoundError = Exception  # type: ignore

    def _dist_version(_name: str) -> str:  # type: ignore
        raise PackageNotFoundError()


def _resolve_version() -> str:
    if _scm_version is not None:
        return _scm_version

    try:
        return _dist_version("sigProfilerPlotting")
    except PackageNotFoundError:
        return "0+unknown"


short_version = _resolve_version()
version = short_version
update = "Managed by setuptools_scm from git tags"

"""The fixed menu of matrices bundled with the package itself, for the
dashboard's "Try an example" buttons -- something for a first-time reader to
explore before they have their own matrix ready.

Every file here is real SigProfilerMatrixGenerator output already shipped
under ``sigProfilerPlotting/examples/input`` (the same files
``plot_example.py`` and the README's quick start use), named explicitly
below rather than discovered by scanning that directory -- the set of
loadable names is fixed at import time and never influenced by a request,
so there is no path-traversal surface to defend the way an arbitrary
"browse the server's filesystem" feature would need.
"""

from __future__ import annotations

from pathlib import Path

from .matrix import MatrixError, load_matrix

_INPUT_DIR = Path(__file__).resolve().parent.parent / "examples" / "input"

# Each covers a different mutation class, and is confirmed by
# tests/test_web_examples.py to load cleanly through the dashboard's own
# strict parser -- not every file under examples/input does (some predate
# it and use an older column layout); only list one here once that is
# checked.
_EXAMPLES: tuple[dict[str, str], ...] = (
    {
        "name": "sbs96",
        "label": "Single-base substitutions",
        "filename": "breast_cancer_samples_example.SBS96.all",
    },
    {
        "name": "sbs288",
        "label": "Strand-tagged substitutions",
        "filename": "test288.all",
    },
    {
        "name": "dbs78",
        "label": "Doublet substitutions",
        "filename": "dbs78_example.all",
    },
    {
        "name": "id83",
        "label": "Small insertions and deletions",
        "filename": "breast_cancer_samples_example.ID94.all",
    },
    {
        "name": "cnv48",
        "label": "Copy-number alterations",
        "filename": "breast_cancer_samples_example.CNV48.all",
    },
    {
        "name": "sv32",
        "label": "Structural variants",
        "filename": "breast_cancer_samples_example.RS32.all",
    },
)
_BY_NAME = {entry["name"]: entry for entry in _EXAMPLES}


def _path_for(entry: dict[str, str]) -> Path:
    return _INPUT_DIR / entry["filename"]


def available() -> list[dict[str, str]]:
    """The bundled examples present in this install, each with the mutation
    context detected from its own rows. A pip install from before these
    were added to package_data (or any other install missing one of these
    files) simply offers fewer buttons rather than a broken one."""
    found = []
    for entry in _EXAMPLES:
        path = _path_for(entry)
        if not path.is_file():
            continue
        try:
            matrix = load_matrix(path)
        except MatrixError:
            continue
        found.append(
            {"name": entry["name"], "label": entry["label"], "context": matrix.context}
        )
    return found


def example_path(name: str) -> Path | None:
    """The bundled file for a known, present example, or None if `name` is
    unrecognized or its file is missing from this install."""
    entry = _BY_NAME.get(name)
    if entry is None:
        return None
    path = _path_for(entry)
    return path if path.is_file() else None

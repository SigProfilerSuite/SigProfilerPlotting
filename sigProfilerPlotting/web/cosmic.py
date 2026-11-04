"""COSMIC SBS96 reference-signature loading and cosine-similarity matching.

Bundled signature catalogues live in web/cosmic_signatures/ as plain
tab-delimited text (one "Type" column of 96 trinucleotide-context labels,
then one column per named signature) -- the exact format COSMIC itself
publishes, and the same 96-context label format the dashboard's own SBS96
renderer already uses, so no reshaping is needed to compare the two.

Read with plain os.path, not importlib.resources: this directory, like
reference_formats, has no __init__.py, which importlib.resources.files()
resolves to a module spec with origin=None on Python 3.9 --
pathlib.Path(None) then raises outright.
"""

from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

from .matrix import MatrixError

_ROOT = Path(os.path.dirname(os.path.abspath(__file__))) / "cosmic_signatures"

# Matches the single base immediately adjacent to the mutated base on each
# side of a "[C>A]"-style substitution, regardless of what precedes or
# follows it in the full label -- e.g. SBS96's own "A[C>A]A", SBS288's
# "T:A[C>A]A" (a strand tag prefix), SBS1536's "AA[C>A]AA" (wider
# pentanucleotide flanks), and SBS4608's "T:AA[C>A]AA" (both) all search to
# the same "A[C>A]A" trinucleotide context.
_INNER_TRINUCLEOTIDE_RE = re.compile(r"([ACGT])\[([ACGT]>[ACGT])\]([ACGT])")


def collapse_to_sbs96(raw_by_label: dict[str, float]) -> dict[str, float]:
    """Sum an arbitrary SBS-context {label: value} mapping down to the 96
    trinucleotide contexts COSMIC's own catalogue is expressed in. A no-op
    (aside from reshaping) for a matrix that is already SBS96.
    """
    collapsed: dict[str, float] = {}
    for label, value in raw_by_label.items():
        match = _INNER_TRINUCLEOTIDE_RE.search(label)
        if match is None:
            continue
        left, substitution, right = match.groups()
        key = f"{left}[{substitution}]{right}"
        collapsed[key] = collapsed.get(key, 0.0) + value
    return collapsed

# The genome build each bundled file covers -- the dashboard offers these as
# the "version" of COSMIC signatures to compare against. All are COSMIC
# v3.6 SBS; only the reference genome (and so the per-context mutation
# opportunity, hence the fitted signature weights) differs between them.
_FILENAME_BY_BUILD = {
    "GRCh37": "COSMIC_v3.6_SBS_GRCh37.txt",
    "GRCh38": "COSMIC_v3.6_SBS_GRCh38.txt",
    "mm9": "COSMIC_v3.6_SBS_mm9.txt",
    "mm10": "COSMIC_v3.6_SBS_mm10.txt",
    "mm39": "COSMIC_v3.6_SBS_mm39.txt",
    "rn6": "COSMIC_v3.6_SBS_rn6.txt",
    "rn7": "COSMIC_v3.6_SBS_rn7.txt",
}

_cache: dict[str, dict[str, dict[str, float]]] = {}


def available_builds() -> list[str]:
    """Genome builds actually bundled with this install, in a fixed,
    human-sensible order (human builds first) rather than filesystem order."""
    return [build for build in _FILENAME_BY_BUILD if (_ROOT / _FILENAME_BY_BUILD[build]).is_file()]


def _load_build(build: str) -> dict[str, dict[str, float]]:
    """{signature_name: {context_label: value}} for one genome build,
    loaded once and cached -- these files don't change at runtime."""
    if build in _cache:
        return _cache[build]
    filename = _FILENAME_BY_BUILD.get(build)
    if filename is None:
        raise MatrixError(
            f"Unknown COSMIC genome build {build!r}. Choose one of: "
            + ", ".join(available_builds())
        )
    path = _ROOT / filename
    if not path.is_file():
        raise MatrixError(f"The COSMIC reference file for {build!r} is not bundled with this install.")

    with open(path, encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise MatrixError(f"The bundled COSMIC reference file for {build!r} is empty.")
    header = rows[0]
    signature_names = header[1:]
    signatures: dict[str, dict[str, float]] = {name: {} for name in signature_names}
    for row in rows[1:]:
        if not row:
            continue
        label = row[0].strip()
        for name, cell in zip(signature_names, row[1:]):
            signatures[name][label] = float(cell)
    _cache[build] = signatures
    return signatures


def _cosine_similarity(a: list[float], b: list[float]) -> float:
    dot = sum(x * y for x, y in zip(a, b))
    norm_a = math.sqrt(sum(x * x for x in a))
    norm_b = math.sqrt(sum(y * y for y in b))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return dot / (norm_a * norm_b)


def top_matches(query_by_label: dict[str, float], build: str, n: int = 5) -> list[dict]:
    """The `n` COSMIC signatures from `build` most cosine-similar to
    `query_by_label` ({96-context label: value}), ranked descending.

    Cosine similarity is scale-invariant (cos(k*a, b) == cos(a, b) for any
    k > 0), so it makes no difference whether the query is raw counts,
    percentages, or proportions -- the caller's chosen display unit need
    not be considered here at all.
    """
    signatures = _load_build(build)
    labels = list(query_by_label.keys())
    query_vector = [query_by_label[label] for label in labels]
    scored = []
    for name, values_by_label in signatures.items():
        signature_vector = [values_by_label.get(label, 0.0) for label in labels]
        score = _cosine_similarity(query_vector, signature_vector)
        scored.append((score, name, values_by_label))
    scored.sort(key=lambda entry: entry[0], reverse=True)
    return [
        {"name": name, "score": score, "values": values_by_label}
        for score, name, values_by_label in scored[:n]
    ]

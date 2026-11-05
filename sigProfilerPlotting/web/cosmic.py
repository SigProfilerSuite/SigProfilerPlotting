"""COSMIC SBS reference-signature loading and cosine-similarity matching.

SigProfilerPlotting does not bundle or redistribute COSMIC's signature
catalogues -- COSMIC's published terms restrict that
(https://www.cosmickb.org/terms/,
https://cancer.sanger.ac.uk/cancergenome/assets/COSMIC_non-commercial_license_May2021.pdf).
Instead, a user who has their own authorized copy of the COSMIC SBS
reference files (downloaded directly from COSMIC) points the dashboard at
the local folder holding them -- see ``configure_catalogue_dir`` -- and
nothing is ever written into this package or its distributed artifacts.

Expected file naming inside that folder, matching COSMIC's own download
naming convention: ``COSMIC_v<version>_SBS_<build>.txt``, tab-delimited,
with a first column of 96 trinucleotide-context labels (COSMIC's own
"Type" column) and one column per named signature -- the same 96-context
label format the dashboard's own SBS96 renderer already uses, so no
reshaping is needed to compare the two.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

from .matrix import MatrixError

_FILENAME_RE = re.compile(r"^COSMIC_v(?P<version>[0-9][0-9A-Za-z.]*)_SBS_(?P<build>\w+)\.txt$")

# Canonical casing for recognized genome builds, and the fixed, human-sensible
# display order (human builds first) rather than filesystem order.
_KNOWN_BUILDS = ("GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6", "rn7")
_BUILD_BY_LOWER = {name.lower(): name for name in _KNOWN_BUILDS}

_catalogue_dir: Path | None = None
_env_dir = os.environ.get("SIGPROFILERPLOTTING_COSMIC_DIR")
if _env_dir:
    _catalogue_dir = Path(_env_dir).expanduser()


def configure_catalogue_dir(path) -> None:
    """Point the dashboard at a local folder of user-supplied COSMIC SBS
    reference files (or clear it with ``None``). Also settable at launch
    time via the ``SIGPROFILERPLOTTING_COSMIC_DIR`` environment variable or
    the launcher's ``--cosmic-dir`` flag.
    """
    global _catalogue_dir
    _catalogue_dir = Path(path).expanduser() if path else None
    _cache.clear()


def catalogue_dir() -> Path | None:
    return _catalogue_dir


@dataclass(frozen=True)
class CatalogueInfo:
    build: str
    version: str
    source_path: str
    checksum: str
    signature_count: int


def _discover() -> dict[str, tuple[Path, str]]:
    """{build: (path, version)} for every recognized file in the configured
    directory. Not cached -- this is a cheap directory listing, and re-
    running it means a file added or removed from the folder after startup
    is picked up on the next request rather than needing a server restart.
    """
    directory = _catalogue_dir
    if directory is None or not directory.is_dir():
        return {}
    found: dict[str, tuple[Path, str]] = {}
    for entry in sorted(directory.iterdir()):
        match = _FILENAME_RE.match(entry.name)
        if not match or not entry.is_file():
            continue
        build = _BUILD_BY_LOWER.get(match.group("build").lower())
        if build is None:
            continue
        found[build] = (entry, match.group("version"))
    return found

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

_cache: dict[str, dict[str, dict[str, float]]] = {}
_metadata_cache: dict[str, CatalogueInfo] = {}


def available_builds() -> list[str]:
    """Genome builds whose reference file is present in the configured
    catalogue directory right now (empty if none is configured), in a
    fixed, human-sensible order (human builds first) rather than
    filesystem order."""
    found = _discover()
    return [build for build in _KNOWN_BUILDS if build in found]


def catalogue_metadata() -> list[CatalogueInfo]:
    """Build/version/source-path/checksum/signature-count for every
    recognized file currently in the configured directory -- so the UI can
    show the user exactly which file and version they are comparing
    against, not just an opaque build name.
    """
    return [_load_metadata(build) for build in available_builds()]


def _checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_metadata(build: str) -> CatalogueInfo:
    cached = _metadata_cache.get(build)
    path, version = _discover()[build]
    if cached is not None and cached.source_path == str(path):
        return cached
    signatures = _load_build(build)
    info = CatalogueInfo(
        build=build,
        version=version,
        source_path=str(path),
        checksum=_checksum(path),
        signature_count=len(signatures),
    )
    _metadata_cache[build] = info
    return info


def _load_build(build: str) -> dict[str, dict[str, float]]:
    """{signature_name: {context_label: value}} for one genome build,
    loaded once and cached by its source path -- a file replaced at the
    same path (different content, same name) still invalidates correctly
    since the cache key folds in its checksum."""
    found = _discover()
    if build not in found:
        raise MatrixError(
            f"Unknown or unavailable COSMIC genome build {build!r}. Choose one of: "
            + ", ".join(available_builds())
        )
    path, _version = found[build]
    cache_key = f"{path}:{_checksum(path)}"
    if cache_key in _cache:
        return _cache[cache_key]

    with open(path, encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or len(rows[0]) < 2:
        raise MatrixError(
            f"The COSMIC reference file for {build!r} ({path}) must have a header row "
            "with a context-label column followed by one column per signature."
        )
    header = rows[0]
    signature_names = [name.strip() for name in header[1:]]
    if any(not name for name in signature_names):
        raise MatrixError(f"The COSMIC reference file for {build!r} ({path}) has an unnamed signature column.")
    signatures: dict[str, dict[str, float]] = {name: {} for name in signature_names}
    for line_number, row in enumerate(rows[1:], 2):
        if not row or all(not str(cell).strip() for cell in row):
            continue
        if len(row) != len(header):
            raise MatrixError(
                f"The COSMIC reference file for {build!r} ({path}) has {len(row)} fields "
                f"on line {line_number}; the header has {len(header)}."
            )
        label = row[0].strip()
        for name, cell in zip(signature_names, row[1:]):
            try:
                value = float(cell)
            except (TypeError, ValueError) as exc:
                raise MatrixError(
                    f"The COSMIC reference file for {build!r} ({path}) has a non-numeric "
                    f"value {cell!r} on line {line_number}, column {name!r}."
                ) from exc
            signatures[name][label] = value
    _cache[cache_key] = signatures
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

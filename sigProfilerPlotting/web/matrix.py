"""Strict matrix loading and automatic SigProfiler context detection.

Uploaded matrices are validated once by streaming through the file with
csv.reader directly off an open file handle (never reading the whole file
into one decoded string, nor materializing a list of every row, nor a
tuple of every cell's float value) and then read again, lazily, one sample
at a time, on each later render request. SigProfilerPlotting's own largest
matrix shape, SBS4608, can realistically have on the order of 50,000
sample columns; holding a Python float for every cell of every sample at
once (hundreds of millions of floats, several times over while parsing)
is not something an ordinary laptop should be asked to do just to upload
the file and look at one sample's plot.
"""

from __future__ import annotations

import csv
import io
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from ..sigProfilerPlotting import INDEX_VALS

# Reuse the package's own canonical header names (sigProfilerPlotting.py's
# INDEX_VALS, e.g. the "index" name pandas writes by default) rather than
# maintaining a second, drift-prone list here -- plus "context", which the
# dashboard accepts as an extra, more descriptive alias. "" is also accepted:
# `DataFrame.to_csv(sep="\t")` (no index name set) writes a blank first
# header cell, which is a very common way this exact kind of matrix gets
# produced -- rejecting it just because the corner cell is blank, when every
# data row still has a real mutation-context label, is needless friction.
MUTATION_TYPE_HEADER_ALIASES = {name.lower() for name in INDEX_VALS} | {"context", ""}


class MatrixError(ValueError):
    """An actionable error in a matrix supplied to the dashboard."""


@dataclass(frozen=True)
class MatrixData:
    context: str
    labels: tuple[str, ...]
    columns: tuple[str, ...]
    # Exactly one of these two is set, never both:
    # - `values`: the whole matrix already in memory (small test fixtures
    #   built directly in Python, and the pandas-DataFrame input path,
    #   which is already fully in memory as numpy arrays before it ever
    #   reaches here -- there is no streaming benefit left to capture).
    # - `_source_path`/`_delimiter`: every uploaded file actually ends up
    #   here (see load_matrix/_load_from_path) -- sample_values() re-reads
    #   just the one requested column from disk instead.
    values: tuple[tuple[float, ...], ...] | None = None
    _source_path: Path | None = None
    _delimiter: str = "\t"

    def sample_values(self, sample: str) -> list[float]:
        try:
            index = self.columns.index(sample)
        except ValueError as exc:
            raise MatrixError(
                f"Unknown sample/signature {sample!r}. Choose one of: "
                + ", ".join(self.columns)
            ) from exc
        if self.values is not None:
            return [row[index] for row in self.values]
        return self._read_column_from_disk(index)

    def _read_column_from_disk(self, index: int) -> list[float]:
        # Already fully validated once, at upload time (see
        # _load_from_path) -- this re-read trusts the file and only
        # converts the one column actually requested, so peak memory here
        # is proportional to one row's width, not the whole matrix.
        by_label: dict[str, float] = {}
        with open(self._source_path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=self._delimiter)
            next(reader, None)  # header
            for row in reader:
                if not row or all(not str(cell).strip() for cell in row):
                    continue
                by_label[row[0].strip()] = float(row[index + 1])
        return [by_label[label] for label in self.labels]

    def raw4608(self, sample: str) -> list[dict[str, object]]:
        return [
            {"mutationType": label, "contribution": value}
            for label, value in zip(self.labels, self.sample_values(sample))
        ]


def _references() -> dict[str, tuple[str, ...]]:
    # sigProfilerPlotting/reference_formats has no __init__.py, so it's an
    # implicit namespace package -- importlib.resources.files() on Python
    # 3.9 resolves that to a spec with origin=None, which pathlib.Path()
    # then rejects outright ("expected str, bytes or os.PathLike object,
    # not NoneType"). Use the same plain os.path-relative-to-__file__
    # approach sigProfilerPlotting.py itself already uses for this exact
    # directory (SPP_REFERENCE), which has no such namespace-package pitfall.
    root = Path(os.path.dirname(os.path.abspath(__file__))).parent / "reference_formats"
    result = {}
    for item in sorted(root.iterdir()):
        if item.name.endswith(".txt"):
            rows = tuple(
                line.strip()
                for line in item.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
            result[item.name[:-4].upper()] = rows
    return result


REFERENCE_CONTEXTS = _references()


def detect_context(labels: Iterable[str]) -> str:
    """Return the one reference context whose complete label set is present."""
    rows = tuple(str(label).strip() for label in labels)
    if not rows:
        raise MatrixError("The matrix contains no mutation-context rows.")
    seen: set[str] = set()
    duplicate = next((label for label in rows if label in seen or seen.add(label)), None)
    if duplicate is not None:
        raise MatrixError(f"Duplicate mutation context {duplicate!r}; each row must be unique.")

    row_set = set(rows)
    matches = [
        name for name, expected in REFERENCE_CONTEXTS.items()
        if len(rows) == len(expected) and row_set == set(expected)
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise MatrixError(
            "The rows match more than one context ("
            + ", ".join(sorted(matches))
            + "); specify a matrix with an unambiguous SigProfiler label set."
        )

    ranked = sorted(
        REFERENCE_CONTEXTS.items(),
        key=lambda pair: (
            len(set(pair[1]) - row_set) + len(row_set - set(pair[1])),
            abs(len(pair[1]) - len(rows)),
        ),
    )
    details = []
    for name, expected in ranked[:3]:
        expected_set = set(expected)
        missing = expected_set - row_set
        extra = row_set - expected_set
        text = f"{name}: {len(missing)} missing, {len(extra)} unrecognized"
        examples = list(sorted(missing))[:2] + list(sorted(extra))[:2]
        if examples:
            text += " (for example " + ", ".join(repr(v) for v in examples) + ")"
        details.append(text)
    raise MatrixError(
        f"Could not detect a supported context from {len(rows)} rows. Closest matches: "
        + "; ".join(details)
        + ". Use the complete, non-duplicated row set from a SigProfiler matrix."
    )


def _decode(source) -> tuple[list[list[str]], str]:
    """Fully materialize a bytes/file-object source as a list of rows.

    Only reached for sources that are already fully in memory by the time
    the caller has them (``bytes``, or a file-like object whose ``.read()``
    the caller already chose to call) -- there is no streaming benefit to
    chase here, so this keeps the simple, original, eager behavior. A real
    uploaded file always arrives as a path instead (see
    ``_materialize_source``/``_load_from_path``), which is what actually
    needs -- and gets -- the memory-bounded treatment.
    """
    name = getattr(source, "name", "uploaded matrix")
    if isinstance(source, bytes):
        raw = source
    elif hasattr(source, "read"):
        raw = source.read()
        if isinstance(raw, str):
            raw = raw.encode("utf-8")
    else:
        raise TypeError("matrix source must be a path, bytes, file object, or DataFrame")
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MatrixError(f"{name} is not UTF-8 text.") from exc
    if not text.strip():
        raise MatrixError(f"{name} is empty.")
    first = next((line for line in text.splitlines() if line.strip()), "")
    delimiter = "\t" if "\t" in first else "," if "," in first else None
    if delimiter is None:
        raise MatrixError(
            "The matrix must be tab-delimited (preferred) or comma-delimited with a header row."
        )
    try:
        return list(csv.reader(io.StringIO(text), delimiter=delimiter)), str(name)
    except csv.Error as exc:
        raise MatrixError(f"Could not parse {name}: {exc}") from exc


def _from_dataframe(frame) -> tuple[list[list[object]], str]:
    columns = [str(c).strip() for c in frame.columns]
    if columns and columns[0].lower() in MUTATION_TYPE_HEADER_ALIASES:
        labels = frame.iloc[:, 0].tolist()
        sample_columns = columns[1:]
        data = frame.iloc[:, 1:].values.tolist()
    else:
        index = frame.index
        if getattr(index, "__class__", type(index)).__name__ == "RangeIndex":
            raise MatrixError(
                "DataFrame input needs mutation contexts in its first named column or index."
            )
        labels = list(index)
        sample_columns = columns
        data = frame.values.tolist()
    rows = [["MutationType", *sample_columns]]
    rows.extend([[label, *values] for label, values in zip(labels, data)])
    return rows, "DataFrame"


def _validate_header(header: list[str]) -> list[str]:
    """Shared header checks for both the eager and streaming loaders.
    Returns the sample/signature column names."""
    width = len(header)
    if width < 2:
        raise MatrixError("The header must contain a context column and at least one sample/signature.")
    if header[0].lower() not in MUTATION_TYPE_HEADER_ALIASES:
        raise MatrixError(
            "The first column must identify mutation contexts (MutationType, Mutation Types, "
            "classification, index, context, or left blank)."
        )
    columns = header[1:]
    if any(not column for column in columns):
        raise MatrixError("Every sample/signature column needs a non-empty name.")
    if len(set(columns)) != len(columns):
        duplicate = next(c for c in columns if columns.count(c) > 1)
        raise MatrixError(f"Duplicate sample/signature column {duplicate!r}; names must be unique.")
    return columns


def _canonical_order(context: str, labels: list[str]) -> tuple[str, ...]:
    order = {label: i for i, label in enumerate(REFERENCE_CONTEXTS[context])}
    return tuple(sorted(labels, key=lambda label: order[label]))


def _load_from_rows(rows: list[list], name: str) -> MatrixData:
    """Build a fully in-memory MatrixData from an already-materialized list
    of rows (the DataFrame and bytes/file-object input paths -- both are
    already fully in memory before reaching this function, so there is no
    streaming benefit to chase here)."""
    if not rows:
        raise MatrixError(f"{name} is empty.")
    header = [str(value).strip() for value in rows[0]]
    width = len(header)
    columns = _validate_header(header)

    labels: list[str] = []
    values: list[tuple[float, ...]] = []
    for line_number, row in enumerate(rows[1:], 2):
        if not row or all(str(cell).strip() == "" for cell in row):
            continue
        if len(row) != width:
            raise MatrixError(
                f"Line {line_number} has {len(row)} fields; the header has {width}. "
                "Check for missing values or stray delimiters."
            )
        label = str(row[0]).strip()
        if not label:
            raise MatrixError(f"Line {line_number} has an empty mutation context.")
        parsed = []
        for column, cell in zip(columns, row[1:]):
            text = str(cell).strip()
            try:
                number = float(text)
            except (TypeError, ValueError) as exc:
                raise MatrixError(
                    f"Non-numeric value {text!r} at line {line_number}, column {column!r}."
                ) from exc
            if not math.isfinite(number):
                raise MatrixError(
                    f"Value at line {line_number}, column {column!r} must be finite; got {text!r}."
                )
            if number < 0:
                raise MatrixError(
                    f"Value at line {line_number}, column {column!r} must be non-negative; got {text!r}."
                )
            parsed.append(number)
        labels.append(label)
        values.append(tuple(parsed))
    context = detect_context(labels)
    order = {label: i for i, label in enumerate(REFERENCE_CONTEXTS[context])}
    sorted_rows = sorted(zip(labels, values), key=lambda pair: order[pair[0]])
    return MatrixData(
        context=context,
        labels=tuple(pair[0] for pair in sorted_rows),
        columns=tuple(columns),
        values=tuple(pair[1] for pair in sorted_rows),
    )


def _load_from_path(path: Path, name: str) -> MatrixData:
    """Validate every cell by streaming through `path` with csv.reader
    directly off the open file handle -- never reading the whole file into
    one decoded string, a list of every row, or a tuple of every cell's
    float value. The returned MatrixData reads lazily from this same path
    afterward (MatrixData._read_column_from_disk), one sample at a time,
    which is the only part of this that matters for real SBS4608-scale
    uploads: a render request only ever pays for the one sample it asked
    for, never for every other sample sitting unused in the matrix.
    """
    try:
        with open(path, encoding="utf-8-sig", newline="") as probe:
            first_line = probe.readline()
    except UnicodeDecodeError as exc:
        raise MatrixError(f"{name} is not UTF-8 text.") from exc
    if not first_line.strip():
        raise MatrixError(f"{name} is empty.")
    delimiter = "\t" if "\t" in first_line else "," if "," in first_line else None
    if delimiter is None:
        raise MatrixError(
            "The matrix must be tab-delimited (preferred) or comma-delimited with a header row."
        )

    labels: list[str] = []
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            try:
                header = next(reader)
            except StopIteration:
                raise MatrixError(f"{name} is empty.")
            header = [str(value).strip() for value in header]
            width = len(header)
            columns = _validate_header(header)

            for line_number, row in enumerate(reader, 2):
                if not row or all(str(cell).strip() == "" for cell in row):
                    continue
                if len(row) != width:
                    raise MatrixError(
                        f"Line {line_number} has {len(row)} fields; the header has {width}. "
                        "Check for missing values or stray delimiters."
                    )
                label = str(row[0]).strip()
                if not label:
                    raise MatrixError(f"Line {line_number} has an empty mutation context.")
                for column, cell in zip(columns, row[1:]):
                    text = str(cell).strip()
                    try:
                        number = float(text)
                    except (TypeError, ValueError) as exc:
                        raise MatrixError(
                            f"Non-numeric value {text!r} at line {line_number}, column {column!r}."
                        ) from exc
                    if not math.isfinite(number):
                        raise MatrixError(
                            f"Value at line {line_number}, column {column!r} must be finite; "
                            f"got {text!r}."
                        )
                    if number < 0:
                        raise MatrixError(
                            f"Value at line {line_number}, column {column!r} must be non-negative; "
                            f"got {text!r}."
                        )
                labels.append(label)
    except UnicodeDecodeError as exc:
        raise MatrixError(f"{name} is not UTF-8 text.") from exc

    context = detect_context(labels)
    return MatrixData(
        context=context,
        labels=_canonical_order(context, labels),
        columns=tuple(columns),
        values=None,
        _source_path=path,
        _delimiter=delimiter,
    )


def _materialize_source(source) -> tuple[Path, str]:
    """A real on-disk path for any (str, Path) source -- the only kind
    load_matrix now reads with the memory-bounded streaming loader. Raises
    if the source is anything else so the caller falls back to the eager,
    already-in-memory path (bytes, file-object, DataFrame)."""
    path = Path(source).expanduser()
    if not path.is_file():
        raise MatrixError(f"Matrix file does not exist: {path}")
    return path, str(path)


def load_matrix(source) -> MatrixData:
    """Load a file/bytes/file-object or pandas DataFrame and validate all cells.

    A path/str is streamed from disk with bounded memory, regardless of
    size -- this is what the dashboard's own upload handling always passes
    (see server.py), by writing the incoming request straight to a temp
    file rather than holding it in memory. bytes/file-object/DataFrame
    sources are already fully in memory by the time they reach here, so
    they use the simpler, original, eager code path unchanged.
    """
    if isinstance(source, (str, Path)):
        path, name = _materialize_source(source)
        return _load_from_path(path, name)
    if hasattr(source, "columns") and hasattr(source, "iloc"):
        rows, name = _from_dataframe(source)
    else:
        rows, name = _decode(source)
    return _load_from_rows(rows, name)

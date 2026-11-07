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
from typing import Iterable, Sequence

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
    # Logical sample columns can be a subset of the physical file columns
    # when structurally incomplete R-padded profiles were omitted.  Entries
    # are zero-based within the value columns (i.e. excluding the context).
    _source_column_indices: tuple[int, ...] = ()
    warnings: tuple[str, ...] = ()

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
        source_index = (
            self._source_column_indices[index]
            if self._source_column_indices
            else index
        )
        return self._read_column_from_disk(source_index)

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
                label = _canonicalize_label(row[0])
                by_label[label] = float(row[index + 1])
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

# Koh et al.'s original ID89 tables used the closed repeat interval 6--9.
# mSigSpectra now deliberately canonicalizes that one spelling to a right-open
# interval so the category remains valid when repeats above 9 are retained.
# Keep this exact, input-boundary alias narrow; no fuzzy label matching.
_LEGACY_LABEL_ALIASES = {"Del(C):R(6,9)": "Del(C):R(6,)"}


def _canonicalize_label(value: object) -> str:
    label = str(value).strip()
    return _LEGACY_LABEL_ALIASES.get(label, label)


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
    nonempty_lines = (line for line in text.splitlines() if line.strip())
    first = next(nonempty_lines, "")
    second = next(nonempty_lines, "")
    delimiter = _detect_delimiter(first, second)
    if delimiter is None:
        raise MatrixError(
            "The matrix must be tab-delimited (preferred) or comma-delimited with a header row."
        )
    try:
        return list(csv.reader(io.StringIO(text), delimiter=delimiter)), str(name)
    except csv.Error as exc:
        raise MatrixError(f"Could not parse {name}: {exc}") from exc


def _detect_delimiter(header_line: str, first_data_line: str = "") -> str | None:
    """Detect a delimiter, including one-sample R row-name exports.

    Such exports have no delimiter at all in their one-cell header, so the
    first data line is the only place where the delimiter can be observed.
    Header-first detection remains authoritative whenever it is available.
    """
    for line in (header_line, first_data_line):
        if "\t" in line:
            return "\t"
        if "," in line:
            return ","
    return None


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


def _prepare_header(header: list[str], first_data_row: Sequence[object] | None) -> list[str]:
    """Normalize R matrices whose top-left context cell was omitted.

    R's row-name exports can write only the sample names in the header while
    still writing the row name as the first field of every data row.  Accept
    that layout only when the first non-empty data row proves the header is
    short by exactly one field.  A recognized context header is never
    reinterpreted as a sample name, so genuinely malformed rows remain errors.
    """
    if (
        first_data_row is not None
        and header
        and header[0].lower() not in MUTATION_TYPE_HEADER_ALIASES
        and len(first_data_row) == len(header) + 1
    ):
        return ["MutationType", *header]
    return header


def _format_column_examples(columns: Sequence[str]) -> str:
    examples = ", ".join(repr(column) for column in columns[:5])
    if len(columns) > 5:
        examples += f", and {len(columns) - 5} more"
    return examples


def _resolve_missing_columns(
    context: str,
    labels: Sequence[str],
    columns: Sequence[str],
    missing_by_column: Sequence[set[str]],
) -> tuple[tuple[int, ...], tuple[str, ...]]:
    """Return usable physical columns and warnings for R ``NA`` padding.

    The Liu compilation combines genuine 476-channel signatures with
    lower-resolution 83-value profiles.  R pads those shorter columns with
    393 ``NA`` cells (the first 83 physical rows contain the original values),
    but those values are *not* ID476 measurements and cannot be scientifically
    zero-filled into the remaining channels.  Omit that one exact structural
    shape, plus wholly empty columns, and reject every other missing-value
    pattern.
    """
    keep: list[int] = []
    lower_resolution: list[str] = []
    empty: list[str] = []
    label_set = set(labels)
    id476_prefix = set(REFERENCE_CONTEXTS.get("ID476", ())[:83])

    for index, (column, missing) in enumerate(zip(columns, missing_by_column)):
        if not missing:
            keep.append(index)
            continue
        if context != "ID476":
            raise MatrixError(
                f"Missing value 'NA' in column {column!r}. NA padding is only "
                "supported for the documented ID476 R-export layouts."
            )
        if missing == label_set:
            empty.append(column)
            continue
        present = label_set - missing
        if present == id476_prefix and len(missing) == 393:
            lower_resolution.append(column)
            continue
        raise MatrixError(
            f"Column {column!r} has {len(present)} of 476 ID476 values. "
            "Only complete columns, wholly empty columns, or the documented "
            "83-value/393-NA R padding layout are supported; other missing "
            "values may indicate a corrupted matrix."
        )

    if not keep:
        raise MatrixError(
            "The matrix has no complete profiles to plot after omitting "
            "incomplete R-padded columns."
        )

    warnings: list[str] = []
    if lower_resolution:
        warnings.append(
            f"Omitted {len(lower_resolution)} lower-resolution "
            f"{'column' if len(lower_resolution) == 1 else 'columns'} with only "
            "83 of 476 values; they are not valid ID476 profiles and were not "
            f"zero-filled: {_format_column_examples(lower_resolution)}."
        )
    if empty:
        warnings.append(
            f"Omitted {len(empty)} entirely missing "
            f"{'column' if len(empty) == 1 else 'columns'}: "
            f"{_format_column_examples(empty)}."
        )
    return tuple(keep), tuple(warnings)


def _parse_cell(cell: object, line_number: int, column: str) -> float | None:
    text = str(cell).strip()
    if text == "NA":
        return None
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
    return number


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
    first_data_row = next(
        (
            row
            for row in rows[1:]
            if row and any(str(cell).strip() for cell in row)
        ),
        None,
    )
    header = _prepare_header(
        [str(value).strip() for value in rows[0]], first_data_row
    )
    width = len(header)
    columns = _validate_header(header)

    labels: list[str] = []
    values: list[tuple[float | None, ...]] = []
    missing_by_column = [set() for _ in columns]
    for line_number, row in enumerate(rows[1:], 2):
        if not row or all(str(cell).strip() == "" for cell in row):
            continue
        if len(row) != width:
            raise MatrixError(
                f"Line {line_number} has {len(row)} fields; the header has {width}. "
                "Check for missing values or stray delimiters."
            )
        label = _canonicalize_label(row[0])
        if not label:
            raise MatrixError(f"Line {line_number} has an empty mutation context.")
        parsed: list[float | None] = []
        for index, (column, cell) in enumerate(zip(columns, row[1:])):
            number = _parse_cell(cell, line_number, column)
            if number is None:
                missing_by_column[index].add(label)
            parsed.append(number)
        labels.append(label)
        values.append(tuple(parsed))
    context = detect_context(labels)
    keep, warnings = _resolve_missing_columns(
        context, labels, columns, missing_by_column
    )
    order = {label: i for i, label in enumerate(REFERENCE_CONTEXTS[context])}
    sorted_rows = sorted(zip(labels, values), key=lambda pair: order[pair[0]])
    return MatrixData(
        context=context,
        labels=tuple(pair[0] for pair in sorted_rows),
        columns=tuple(columns[index] for index in keep),
        values=tuple(
            tuple(pair[1][index] for index in keep)  # type: ignore[misc]
            for pair in sorted_rows
        ),
        warnings=warnings,
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
            first_data_line = next((line for line in probe if line.strip()), "")
    except UnicodeDecodeError as exc:
        raise MatrixError(f"{name} is not UTF-8 text.") from exc
    if not first_line.strip():
        raise MatrixError(f"{name} is empty.")
    delimiter = _detect_delimiter(first_line, first_data_line)
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
            first_data: tuple[int, list[str]] | None = None
            for line_number, row in enumerate(reader, 2):
                if row and any(str(cell).strip() for cell in row):
                    first_data = (line_number, row)
                    break
            header = _prepare_header(
                header, first_data[1] if first_data is not None else None
            )
            width = len(header)
            columns = _validate_header(header)
            missing_by_column = [set() for _ in columns]

            pending = [first_data] if first_data is not None else []
            row_iterator = iter(reader)
            for line_number, row in pending:
                if len(row) != width:
                    raise MatrixError(
                        f"Line {line_number} has {len(row)} fields; the header has {width}. "
                        "Check for missing values or stray delimiters."
                    )
                label = _canonicalize_label(row[0])
                if not label:
                    raise MatrixError(f"Line {line_number} has an empty mutation context.")
                for index, (column, cell) in enumerate(zip(columns, row[1:])):
                    if _parse_cell(cell, line_number, column) is None:
                        missing_by_column[index].add(label)
                labels.append(label)

            start_line = (first_data[0] + 1) if first_data is not None else 2
            for line_number, row in enumerate(row_iterator, start_line):
                if not row or all(str(cell).strip() == "" for cell in row):
                    continue
                if len(row) != width:
                    raise MatrixError(
                        f"Line {line_number} has {len(row)} fields; the header has {width}. "
                        "Check for missing values or stray delimiters."
                    )
                label = _canonicalize_label(row[0])
                if not label:
                    raise MatrixError(f"Line {line_number} has an empty mutation context.")
                for index, (column, cell) in enumerate(zip(columns, row[1:])):
                    if _parse_cell(cell, line_number, column) is None:
                        missing_by_column[index].add(label)
                labels.append(label)
    except UnicodeDecodeError as exc:
        raise MatrixError(f"{name} is not UTF-8 text.") from exc

    context = detect_context(labels)
    keep, warnings = _resolve_missing_columns(
        context, labels, columns, missing_by_column
    )
    return MatrixData(
        context=context,
        labels=_canonical_order(context, labels),
        columns=tuple(columns[index] for index in keep),
        values=None,
        _source_path=path,
        _delimiter=delimiter,
        _source_column_indices=keep,
        warnings=warnings,
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

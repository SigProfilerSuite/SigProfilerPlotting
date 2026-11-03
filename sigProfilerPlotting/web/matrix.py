"""Strict matrix loading and automatic SigProfiler context detection."""

from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Iterable


class MatrixError(ValueError):
    """An actionable error in a matrix supplied to the dashboard."""


@dataclass(frozen=True)
class MatrixData:
    context: str
    labels: tuple[str, ...]
    columns: tuple[str, ...]
    values: tuple[tuple[float, ...], ...]

    def sample_values(self, sample: str) -> list[float]:
        try:
            index = self.columns.index(sample)
        except ValueError as exc:
            raise MatrixError(
                f"Unknown sample/signature {sample!r}. Choose one of: "
                + ", ".join(self.columns)
            ) from exc
        return [row[index] for row in self.values]

    def raw4608(self, sample: str) -> list[dict[str, object]]:
        return [
            {"mutationType": label, "contribution": value}
            for label, value in zip(self.labels, self.sample_values(sample))
        ]


def _references() -> dict[str, tuple[str, ...]]:
    root = files("sigProfilerPlotting.reference_formats")
    result = {}
    for item in root.iterdir():
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
    name = getattr(source, "name", "uploaded matrix")
    if isinstance(source, (str, Path)):
        path = Path(source).expanduser()
        if not path.is_file():
            raise MatrixError(f"Matrix file does not exist: {path}")
        raw = path.read_bytes()
        name = str(path)
    elif isinstance(source, bytes):
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
    first_aliases = {"mutationtype", "mutation types", "classification", "context"}
    if columns and columns[0].lower() in first_aliases:
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


def load_matrix(source) -> MatrixData:
    """Load a file/bytes/file-object or pandas DataFrame and validate all cells."""
    if hasattr(source, "columns") and hasattr(source, "iloc"):
        rows, name = _from_dataframe(source)
    else:
        rows, name = _decode(source)
    if not rows:
        raise MatrixError(f"{name} is empty.")
    width = len(rows[0])
    if width < 2:
        raise MatrixError("The header must contain a context column and at least one sample/signature.")
    header = [str(value).strip() for value in rows[0]]
    aliases = {"mutationtype", "mutation types", "classification", "context"}
    if header[0].lower() not in aliases:
        raise MatrixError(
            "The first column must identify mutation contexts (MutationType, Mutation Types, "
            "classification, or context)."
        )
    columns = header[1:]
    if any(not column for column in columns):
        raise MatrixError("Every sample/signature column needs a non-empty name.")
    if len(set(columns)) != len(columns):
        duplicate = next(c for c in columns if columns.count(c) > 1)
        raise MatrixError(f"Duplicate sample/signature column {duplicate!r}; names must be unique.")

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

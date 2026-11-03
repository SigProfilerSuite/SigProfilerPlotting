import pytest

from sigProfilerPlotting.web.matrix import MatrixError, detect_context, load_matrix

SBS6_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def _write(path, labels, columns=("Sample_1",), values=None):
    values = values or {label: [1.0] * len(columns) for label in labels}
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("MutationType\t" + "\t".join(columns) + "\n")
        for label in labels:
            row = values.get(label, [0.0] * len(columns))
            fh.write(label + "\t" + "\t".join(str(v) for v in row) + "\n")
    return path


def test_loads_a_well_formed_sbs6_matrix(tmp_path):
    path = _write(tmp_path / "ok.tsv", SBS6_LABELS, columns=("S1", "S2"))
    matrix = load_matrix(path)
    assert matrix.context == "SBS6"
    assert matrix.columns == ("S1", "S2")
    assert len(matrix.labels) == 6


def test_duplicate_mutation_context_row_is_rejected(tmp_path):
    path = tmp_path / "dup.tsv"
    path.write_text("MutationType\tS1\nC>A\t1\nC>A\t2\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="Duplicate mutation context"):
        load_matrix(path)


def test_duplicate_sample_column_is_rejected(tmp_path):
    path = tmp_path / "dupcol.tsv"
    path.write_text("MutationType\tS1\tS1\nC>A\t1\t2\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="Duplicate sample/signature column"):
        load_matrix(path)


def test_missing_context_rows_give_an_actionable_message(tmp_path):
    path = _write(tmp_path / "missing.tsv", SBS6_LABELS[:-1])  # drop T>G
    with pytest.raises(MatrixError, match="Could not detect a supported context"):
        load_matrix(path)


def test_unexpected_extra_row_gives_an_actionable_message(tmp_path):
    labels = SBS6_LABELS + ["Z>Z"]
    path = _write(tmp_path / "extra.tsv", labels)
    with pytest.raises(MatrixError, match="Could not detect a supported context"):
        load_matrix(path)


def test_malformed_row_with_wrong_field_count_is_rejected(tmp_path):
    path = tmp_path / "malformed.tsv"
    path.write_text("MutationType\tS1\tS2\nC>A\t1\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="has 2 fields; the header has 3"):
        load_matrix(path)


def test_negative_value_is_rejected(tmp_path):
    path = _write(
        tmp_path / "negative.tsv",
        SBS6_LABELS,
        values={"C>A": [-1.0]},
    )
    with pytest.raises(MatrixError, match="must be non-negative"):
        load_matrix(path)


def test_nonnumeric_value_is_rejected(tmp_path):
    path = tmp_path / "nonnumeric.tsv"
    rows = "\n".join(f"{label}\tabc" if label == "C>A" else f"{label}\t1" for label in SBS6_LABELS)
    path.write_text(f"MutationType\tS1\n{rows}\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="Non-numeric value"):
        load_matrix(path)


def test_nonfinite_value_is_rejected(tmp_path):
    path = tmp_path / "nonfinite.tsv"
    rows = "\n".join(f"{label}\tnan" if label == "C>A" else f"{label}\t1" for label in SBS6_LABELS)
    path.write_text(f"MutationType\tS1\n{rows}\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="must be finite"):
        load_matrix(path)


def test_all_zero_matrix_still_loads(tmp_path):
    path = _write(tmp_path / "zero.tsv", SBS6_LABELS, values={label: [0.0] for label in SBS6_LABELS})
    matrix = load_matrix(path)
    assert all(v == 0.0 for row in matrix.values for v in row)


def test_wide_matrix_with_many_sample_columns(tmp_path):
    columns = tuple(f"Sample_{i}" for i in range(200))
    path = _write(tmp_path / "wide.tsv", SBS6_LABELS, columns=columns)
    matrix = load_matrix(path)
    assert matrix.columns == columns
    assert len(matrix.values[0]) == 200


def test_unusual_sample_names_are_preserved_verbatim(tmp_path):
    columns = ("Sample A", "sig,with,commas", "Sampleα")
    path = _write(tmp_path / "unusual.tsv", SBS6_LABELS, columns=columns)
    matrix = load_matrix(path)
    assert matrix.columns == columns


def test_shuffled_row_order_is_normalized_to_reference_order(tmp_path):
    shuffled = list(reversed(SBS6_LABELS))
    path = _write(tmp_path / "shuffled.tsv", shuffled)
    matrix = load_matrix(path)
    assert matrix.labels == tuple(SBS6_LABELS)


def test_large_sample_count_loads_without_blowing_up(tmp_path):
    columns = tuple(f"S{i}" for i in range(2000))
    path = _write(tmp_path / "large.tsv", SBS6_LABELS, columns=columns)
    matrix = load_matrix(path)
    assert len(matrix.columns) == 2000
    assert len(matrix.values) == 6


def test_detect_context_reports_ambiguous_matches(monkeypatch):
    import sigProfilerPlotting.web.matrix as matrix_module

    monkeypatch.setattr(
        matrix_module,
        "REFERENCE_CONTEXTS",
        {"FAKE_A": ("X", "Y"), "FAKE_B": ("X", "Y")},
    )
    with pytest.raises(MatrixError, match="more than one context"):
        detect_context(["X", "Y"])


def test_empty_matrix_is_rejected(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(MatrixError, match="is empty"):
        load_matrix(path)


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(MatrixError, match="does not exist"):
        load_matrix(tmp_path / "nope.tsv")


def test_wrong_first_header_is_rejected(tmp_path):
    path = tmp_path / "badheader.tsv"
    path.write_text("NotAHeader\tS1\nC>A\t1\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="must identify mutation contexts"):
        load_matrix(path)

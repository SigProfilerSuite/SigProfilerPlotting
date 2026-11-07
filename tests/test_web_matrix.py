import pytest

from sigProfilerPlotting.web.matrix import (
    REFERENCE_CONTEXTS,
    MatrixError,
    detect_context,
    load_matrix,
)

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
    assert all(v == 0.0 for v in matrix.sample_values(matrix.columns[0]))


def test_wide_matrix_with_many_sample_columns(tmp_path):
    columns = tuple(f"Sample_{i}" for i in range(200))
    values = {label: [float(i) for i in range(200)] for label in SBS6_LABELS}
    path = _write(tmp_path / "wide.tsv", SBS6_LABELS, columns=columns, values=values)
    matrix = load_matrix(path)
    assert matrix.columns == columns
    # Every one of the 200 columns is independently addressable and lands
    # on the right value, not just the first/last few.
    assert matrix.sample_values("Sample_0")[0] == 0.0
    assert matrix.sample_values("Sample_199")[0] == 199.0


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
    assert len(matrix.labels) == 6
    assert len(matrix.sample_values("S0")) == 6


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


def test_index_header_is_accepted(tmp_path):
    # pandas' default DataFrame.to_csv() output names an unlabeled index
    # column "index" -- the package's own sigProfilerPlotting.py already
    # treats this as a valid mutation-context header (INDEX_VALS), and the
    # dashboard parser must match that, not reject real-world exports.
    path = _write(tmp_path / "index_header.tsv", SBS6_LABELS, columns=("S1",))
    path.write_text(path.read_text(encoding="utf-8").replace("MutationType", "index", 1))
    matrix = load_matrix(path)
    assert matrix.context == "SBS6"


def test_blank_first_header_cell_is_accepted(tmp_path):
    # DataFrame.to_csv(sep="\t") with no index name set (a very common real
    # -world export, e.g. pd.DataFrame(..., index=mutation_types).to_csv())
    # writes a *blank* first header cell, not "index" or "MutationType" --
    # every data row still carries a real mutation-context label in its
    # first column, so this should not be rejected either.
    path = _write(tmp_path / "blank_header.tsv", SBS6_LABELS, columns=("S1",))
    path.write_text(path.read_text(encoding="utf-8").replace("MutationType", "", 1))
    matrix = load_matrix(path)
    assert matrix.context == "SBS6"


@pytest.mark.parametrize("source_kind", ["path", "bytes"])
def test_r_header_without_top_left_context_cell_is_accepted(tmp_path, source_kind):
    text = "S1\tS2\n" + "\n".join(
        f"{label}\t{index}\t{index + 10}"
        for index, label in enumerate(SBS6_LABELS)
    )
    source = text.encode("utf-8")
    if source_kind == "path":
        source = tmp_path / "r-export.tsv"
        source.write_text(text, encoding="utf-8")

    matrix = load_matrix(source)

    assert matrix.context == "SBS6"
    assert matrix.columns == ("S1", "S2")
    assert matrix.sample_values("S2") == [float(i + 10) for i in range(6)]


def test_r_header_is_not_inferred_when_data_is_more_than_one_field_wider(tmp_path):
    path = tmp_path / "too-wide.tsv"
    path.write_text("S1\nC>A\t1\t2\n", encoding="utf-8")
    with pytest.raises(MatrixError, match="header must contain"):
        load_matrix(path)


@pytest.mark.parametrize("source_kind", ["path", "bytes"])
def test_legacy_id89_repeat_label_is_narrowly_canonicalized(tmp_path, source_kind):
    labels = list(REFERENCE_CONTEXTS["ID89"])
    labels[labels.index("Del(C):R(6,)")] = "Del(C):R(6,9)"
    text = "LegacySignature\n" + "\n".join(
        f"{label}\t{index + 1}" for index, label in enumerate(labels)
    )
    source = text.encode("utf-8")
    if source_kind == "path":
        source = tmp_path / "legacy-id89.tsv"
        source.write_text(text, encoding="utf-8")

    matrix = load_matrix(source)

    assert matrix.context == "ID89"
    assert matrix.labels == REFERENCE_CONTEXTS["ID89"]
    assert "Del(C):R(6,9)" not in matrix.labels
    assert len(matrix.sample_values("LegacySignature")) == 89


def _id476_r_padded_text():
    lines = ["LowerResolution\tComplete476\tEntirelyMissing"]
    for index, label in enumerate(REFERENCE_CONTEXTS["ID476"]):
        lower = str(index + 100) if index < 83 else "NA"
        lines.append(f"{label}\t{lower}\t{index + 1}\tNA")
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("source_kind", ["path", "bytes"])
def test_id476_r_padding_omits_incomplete_profiles_with_warnings(
    tmp_path, source_kind
):
    text = _id476_r_padded_text()
    source = text.encode("utf-8")
    if source_kind == "path":
        source = tmp_path / "id476-r-export.tsv"
        source.write_text(text, encoding="utf-8")

    matrix = load_matrix(source)

    assert matrix.context == "ID476"
    assert matrix.columns == ("Complete476",)
    assert len(matrix.labels) == 476
    assert matrix.sample_values("Complete476") == [float(i + 1) for i in range(476)]
    assert len(matrix.warnings) == 2
    assert "LowerResolution" in matrix.warnings[0]
    assert "83 of 476" in matrix.warnings[0]
    assert "EntirelyMissing" in matrix.warnings[1]


def test_id476_arbitrary_na_pattern_is_rejected(tmp_path):
    labels = REFERENCE_CONTEXTS["ID476"]
    lines = ["MutationType\tBroken"]
    for index, label in enumerate(labels):
        value = "NA" if index == 100 else "1"
        lines.append(f"{label}\t{value}")
    path = tmp_path / "broken-id476.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(MatrixError, match="corrupted matrix"):
        load_matrix(path)


def test_na_is_rejected_outside_documented_id476_layout(tmp_path):
    lines = ["MutationType\tS1"] + [
        f"{label}\t{'NA' if index == 0 else '1'}"
        for index, label in enumerate(SBS6_LABELS)
    ]
    path = tmp_path / "sbs-na.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(MatrixError, match="only supported.*ID476"):
        load_matrix(path)


def test_id476_with_no_complete_profiles_is_rejected(tmp_path):
    text = _id476_r_padded_text().replace("\tComplete476", "").replace(
        "\tEntirelyMissing", ""
    )
    lines = []
    for line in text.splitlines():
        fields = line.split("\t")
        lines.append("\t".join(fields[:2]))
    path = tmp_path / "only-lower-resolution.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(MatrixError, match="no complete profiles"):
        load_matrix(path)

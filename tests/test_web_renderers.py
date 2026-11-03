import pytest

from sigProfilerPlotting.web.matrix import MatrixData, load_matrix
from sigProfilerPlotting.web.renderers import render_sample

SBS6_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def _write(path, labels, columns=("S1",), values=None):
    values = values or {label: [1.0] * len(columns) for label in labels}
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("MutationType\t" + "\t".join(columns) + "\n")
        for label in labels:
            row = values.get(label, [0.0] * len(columns))
            fh.write(label + "\t" + "\t".join(str(v) for v in row) + "\n")
    return path


def _sbs96_matrix(values):
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["SBS96"])
    return MatrixData(
        context="SBS96",
        labels=tuple(labels),
        columns=("S1",),
        values=tuple((values.get(label, 0.0),) for label in labels),
    )


def test_render_sample_counts_mode_reports_raw_values(tmp_path):
    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS, values={"C>A": [4.0], **{l: [1.0] for l in SBS6_LABELS[1:]}})
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="counts")
    assert figure["kind"] == "standard"
    ys = [y for trace in figure["figure"]["traces"] for y in trace["y"]]
    assert sum(ys) == pytest.approx(9.0)


def test_render_sample_percent_mode_sums_to_100(tmp_path):
    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS, values={"C>A": [4.0], **{l: [1.0] for l in SBS6_LABELS[1:]}})
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="percent")
    ys = [y for trace in figure["figure"]["traces"] for y in trace["y"]]
    assert sum(ys) == pytest.approx(100.0)


def test_render_sample_percent_mode_handles_all_zero_total(tmp_path):
    path = _write(tmp_path / "zero.tsv", SBS6_LABELS, values={l: [0.0] for l in SBS6_LABELS})
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="percent")
    ys = [y for trace in figure["figure"]["traces"] for y in trace["y"]]
    assert all(y == 0.0 for y in ys)


def test_render_sample_proportion_mode_sums_to_one(tmp_path):
    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS, values={"C>A": [4.0], **{l: [1.0] for l in SBS6_LABELS[1:]}})
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="proportion")
    ys = [y for trace in figure["figure"]["traces"] for y in trace["y"]]
    assert sum(ys) == pytest.approx(1.0)
    # Same distribution as percent mode, just on a 0-1 scale instead of 0-100.
    percent_ys = [
        y for trace in render_sample(matrix, "S1", unit="percent")["figure"]["traces"] for y in trace["y"]
    ]
    assert ys == pytest.approx([p / 100.0 for p in percent_ys])


def test_render_sample_proportion_mode_handles_all_zero_total(tmp_path):
    path = _write(tmp_path / "zero.tsv", SBS6_LABELS, values={l: [0.0] for l in SBS6_LABELS})
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="proportion")
    ys = [y for trace in figure["figure"]["traces"] for y in trace["y"]]
    assert all(y == 0.0 for y in ys)


def test_unknown_sample_raises_actionable_error(tmp_path):
    from sigProfilerPlotting.web.matrix import MatrixError

    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS)
    matrix = load_matrix(path)
    with pytest.raises(MatrixError, match="Unknown sample/signature"):
        render_sample(matrix, "DoesNotExist")


def test_sbs96_renderer_groups_by_substitution_type():
    values = {label: 1.0 for label in [
        "A[C>A]A", "A[C>G]A", "A[C>T]A", "A[T>A]A", "A[T>C]A", "A[T>G]A",
    ]}
    matrix = _sbs96_matrix(values)
    figure = render_sample(matrix, "S1", unit="counts")
    groups = {trace["name"] for trace in figure["figure"]["traces"]}
    assert groups == {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"}


def test_cnv48_renderer_groups_by_cnv_class(tmp_path):
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["CNV48"])
    path = _write(tmp_path / "cnv.tsv", labels)
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="counts")
    groups = {trace["name"] for trace in figure["figure"]["traces"]}
    assert groups.issubset({"homdel", "LOH", "het"})
    assert len(groups) > 0


def test_generic_renderer_used_for_a_context_without_a_specialized_renderer(tmp_path):
    # SBS18 (T/U/N x 6 substitutions) has no dedicated renderer in renderers.py.
    labels = ["T:C>A", "T:C>G", "T:C>T", "T:T>A", "T:T>C", "T:T>G",
              "U:C>A", "U:C>G", "U:C>T", "U:T>A", "U:T>C", "U:T>G",
              "N:C>A", "N:C>G", "N:C>T", "N:T>A", "N:T>C", "N:T>G"]
    path = _write(tmp_path / "sbs18.tsv", labels)
    matrix = load_matrix(path)
    assert matrix.context == "SBS18"
    figure = render_sample(matrix, "S1", unit="counts")
    assert figure["kind"] == "standard"
    groups = {trace["name"] for trace in figure["figure"]["traces"]}
    assert groups == {"T", "U", "N"}


def test_generic_renderer_does_not_group_by_a_leading_numeric_index(tmp_path):
    # ID28 labels like "1:Del:C:0" -- the leading "1" is a positional index,
    # not a semantic family; grouping on it would produce meaningless "1"-style
    # legend entries instead of "Del"/"Ins"/etc.
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["ID28"])
    path = _write(tmp_path / "id28.tsv", labels)
    matrix = load_matrix(path)
    assert matrix.context == "ID28"
    figure = render_sample(matrix, "S1", unit="counts")
    groups = {trace["name"] for trace in figure["figure"]["traces"]}
    assert "1" not in groups
    assert "Del" in groups and "Ins" in groups


def test_sbs4608_renderer_builds_three_variants():
    labels = []
    for strand in "TUN":
        for sub in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]:
            for left in [a + b for a in "ACGT" for b in "ACGT"]:
                for right in [a + b for a in "ACGT" for b in "ACGT"]:
                    labels.append(f"{strand}:{left}[{sub}]{right}")
    matrix = MatrixData(
        context="SBS4608",
        labels=tuple(labels),
        columns=("S1",),
        values=tuple((1.0,) for _ in labels),
    )
    figure = render_sample(matrix, "S1", unit="counts")
    assert figure["kind"] == "sbs4608"
    assert set(figure["variants"].keys()) == {"default", "strand_bias", "genic_bias"}
    for variant in figure["variants"].values():
        assert "traces" in variant
        assert "layout" in variant


def test_layout_legend_does_not_collide_with_left_aligned_title(tmp_path):
    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS)
    matrix = load_matrix(path)
    figure = render_sample(matrix, "S1", unit="counts")
    layout = figure["figure"]["layout"]
    assert layout["title"]["xanchor"] == "left"
    # The legend must not share the title's left edge at the same height.
    assert not (layout["legend"]["x"] == layout["title"]["x"])


def test_invalid_unit_is_rejected(tmp_path):
    from sigProfilerPlotting.web.matrix import MatrixError

    path = _write(tmp_path / "sbs6.tsv", SBS6_LABELS)
    matrix = load_matrix(path)
    with pytest.raises(MatrixError, match="unit must be"):
        render_sample(matrix, "S1", unit="bogus")

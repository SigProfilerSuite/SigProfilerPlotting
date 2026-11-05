import pytest

from sigProfilerPlotting.web import cosmic
from sigProfilerPlotting.web.matrix import MatrixData, load_matrix
from sigProfilerPlotting.web.renderers import render_sample

SBS6_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

# SigProfilerPlotting does not bundle or redistribute COSMIC's signature
# catalogues -- every test below configures its own small, synthetic,
# user-supplied-style catalogue directory instead of relying on real
# bundled data.
_COSMIC_SIGNATURE_NAMES = ["SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS13"]


@pytest.fixture(autouse=True)
def _cosmic_catalogue(tmp_path):
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["SBS96"])
    directory = tmp_path / "cosmic"
    directory.mkdir()
    header = ["Type"] + _COSMIC_SIGNATURE_NAMES
    lines = ["\t".join(header)]
    for idx, label in enumerate(labels):
        row = [label] + [
            str(((idx * 7 + col * 13 + 3) % 23) / 23.0) for col in range(len(_COSMIC_SIGNATURE_NAMES))
        ]
        lines.append("\t".join(row))
    (directory / "COSMIC_v3.4_SBS_GRCh37.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    cosmic.configure_catalogue_dir(directory)
    yield
    cosmic.configure_catalogue_dir(None)


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


def test_sbs96_renderer_uses_the_classic_display_order_and_colors():
    from sigProfilerPlotting.web.renderers import SBS96_DISPLAY_ORDER, SBS_COLORS

    # One count per label so each bar's position tells us which label landed
    # there: give "A[C>A]A" (the first label in the classic display order)
    # a distinct value to confirm the renderer re-sorts away from
    # REFERENCE_CONTEXTS["SBS96"]'s own (differently ordered) row order.
    values = {label: 1.0 for label in SBS96_DISPLAY_ORDER}
    values["A[C>A]A"] = 7.0
    matrix = _sbs96_matrix(values)
    figure = render_sample(matrix, "S1", unit="counts")
    assert figure["kind"] == "standard"

    traces = figure["figure"]["traces"]
    assert len(traces) == 1  # a single bar trace, not one-trace-per-group
    trace = traces[0]
    assert len(trace["y"]) == 96
    assert trace["y"][0] == pytest.approx(7.0)  # SBS96_DISPLAY_ORDER[0] == "A[C>A]A"
    assert list(trace["customdata"]) == SBS96_DISPLAY_ORDER

    # Each bar's color matches its own substitution-type group.
    for label, color in zip(SBS96_DISPLAY_ORDER, trace["marker"]["color"]):
        substitution = label[2:5]
        assert color == SBS_COLORS[substitution]

    # Six colored header-band shapes (one per substitution type) instead of
    # a legend, and no legend is shown.
    layout = figure["figure"]["layout"]
    assert layout["showlegend"] is False
    assert "legend" not in layout
    assert len(layout["shapes"]) == 6
    header_labels = {a["text"] for a in layout["annotations"] if a["text"].strip("<b>/") in
                     {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"}}
    assert header_labels == {"<b>C>A</b>", "<b>C>G</b>", "<b>C>T</b>", "<b>T>A</b>", "<b>T>C</b>", "<b>T>G</b>"}


@pytest.mark.parametrize("unit, expected_title", [
    ("counts", "Number of Single Base Substitutions"),
    ("percent", "Percentage of Single Base Substitutions"),
    ("proportion", "Proportion of Single Base Substitutions"),
])
def test_sbs96_renderer_uses_the_classic_axis_title_per_unit(unit, expected_title):
    matrix = _sbs96_matrix({"A[C>A]A": 3.0})
    figure = render_sample(matrix, "S1", unit=unit)
    assert figure["figure"]["layout"]["yaxis"]["title"]["text"] == expected_title


def test_sbs96_renderer_percent_and_proportion_sum_correctly():
    matrix = _sbs96_matrix({"A[C>A]A": 3.0, "A[C>G]A": 1.0})
    percent_ys = render_sample(matrix, "S1", unit="percent")["figure"]["traces"][0]["y"]
    proportion_ys = render_sample(matrix, "S1", unit="proportion")["figure"]["traces"][0]["y"]
    assert sum(percent_ys) == pytest.approx(100.0)
    assert sum(proportion_ys) == pytest.approx(1.0)


def test_sbs96_compare_cosmic_attaches_ranked_matches_alongside_the_normal_kind():
    from sigProfilerPlotting.web import cosmic

    reference = cosmic._load_build("GRCh37")["SBS4"]
    matrix = _sbs96_matrix(reference)

    result = render_sample(matrix, "S1", unit="counts", compare_cosmic="GRCh37")

    # The comparison is an addition, not a different response shape -- SBS96
    # still renders through its normal "standard" kind.
    assert result["kind"] == "standard"
    assert "traces" in result["figure"] and "layout" in result["figure"]

    compare = result["cosmic_compare"]
    assert compare["build"] == "GRCh37"
    assert len(compare["matches"]) == 5
    # The query *is* SBS4, so it must rank itself first.
    assert compare["matches"][0]["name"] == "SBS4"
    assert compare["matches"][0]["score"] == pytest.approx(1.0, abs=1e-9)
    scores = [m["score"] for m in compare["matches"]]
    assert scores == sorted(scores, reverse=True)
    for match in compare["matches"]:
        assert "traces" in match["figure"] and "layout" in match["figure"]


@pytest.mark.parametrize("context", ["SBS288", "SBS1536", "SBS4608"])
def test_cosmic_compare_collapses_other_sbs_contexts_to_96(context):
    """A matrix in a wider SBS context (strand-tagged and/or a bigger
    flanking window) still gets a comparison, collapsed down to the 96
    trinucleotide contexts COSMIC's own catalogue uses."""
    import sigProfilerPlotting.web.matrix as matrix_module

    if context == "SBS4608":
        from sigProfilerPlotting.sbs4608.sbs4608_data import DINUCS, SUBS

        labels = [
            f"{strand}:{left}[{sub}]{right}"
            for strand in "TUN"
            for sub in SUBS
            for left in DINUCS
            for right in DINUCS
        ]
    else:
        labels = list(matrix_module.REFERENCE_CONTEXTS[context])
    matrix = MatrixData(
        context=context,
        labels=tuple(labels),
        columns=("S1",),
        values=tuple((1.0,) for _ in labels),
    )

    result = render_sample(matrix, "S1", unit="counts", compare_cosmic="GRCh37")

    compare = result["cosmic_compare"]
    assert compare["build"] == "GRCh37"
    assert len(compare["matches"]) == 5
    scores = [m["score"] for m in compare["matches"]]
    assert scores == sorted(scores, reverse=True)
    # A uniform input collapses to a uniform 96-context vector, which some
    # signature should resemble at least loosely.
    assert scores[0] > 0


def test_cosmic_compare_unit_follows_percent_but_falls_back_from_counts():
    matrix = _sbs96_matrix({"A[C>A]A": 3.0})

    percent_result = render_sample(matrix, "S1", unit="percent", compare_cosmic="GRCh37")
    counts_result = render_sample(matrix, "S1", unit="counts", compare_cosmic="GRCh37")
    proportion_result = render_sample(matrix, "S1", unit="proportion", compare_cosmic="GRCh37")

    def match_ytitle(result):
        return result["cosmic_compare"]["matches"][0]["figure"]["layout"]["yaxis"]["title"]["text"]

    assert match_ytitle(percent_result) == "Percentage of Single Base Substitutions"
    # "Counts" has no meaning for a reference signature -- falls back to
    # proportion rather than, say, erroring or showing a meaningless count.
    assert match_ytitle(counts_result) == "Proportion of Single Base Substitutions"
    assert match_ytitle(proportion_result) == "Proportion of Single Base Substitutions"


def test_sbs96_compare_cosmic_rejects_an_unknown_build():
    from sigProfilerPlotting.web.matrix import MatrixError

    matrix = _sbs96_matrix({"A[C>A]A": 3.0})
    with pytest.raises(MatrixError):
        render_sample(matrix, "S1", unit="counts", compare_cosmic="not-a-real-build")


def test_sbs96_without_compare_cosmic_is_unaffected():
    matrix = _sbs96_matrix({"A[C>A]A": 3.0})
    result = render_sample(matrix, "S1", unit="counts")
    assert result["kind"] == "standard"


def test_sbs96_figure_keeps_consistent_header_spacing_at_any_height():
    # _sbs96_figure's header-band/title offsets are computed in absolute
    # pixels so a short figure (the COSMIC comparison panels) doesn't
    # collide them -- assert the computed pixel gap from the axes top is
    # identical regardless of `height`, the invariant that fix relies on.
    import sigProfilerPlotting.web.matrix as matrix_module
    from sigProfilerPlotting.web.renderers import _sbs96_figure

    raw_by_label = {label: 0.0 for label in matrix_module.REFERENCE_CONTEXTS["SBS96"]}
    raw_by_label["A[C>A]A"] = 3.0
    tall = _sbs96_figure(raw_by_label, "S1", "counts", height=560)
    short = _sbs96_figure(raw_by_label, "S1", "counts", height=300)

    def band_pixels_above_axes(figure):
        margin_t = figure["layout"]["margin"]["t"]
        margin_b = figure["layout"]["margin"]["b"]
        axes_height = figure["layout"]["height"] - margin_t - margin_b
        y0 = figure["layout"]["shapes"][0]["y0"]
        return (y0 - 1) * axes_height

    assert band_pixels_above_axes(tall) == pytest.approx(band_pixels_above_axes(short), abs=0.01)


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


def _sbs1536_matrix(values=None):
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["SBS1536"])
    values = values or {label: float(i % 13 + 1) for i, label in enumerate(labels)}
    return MatrixData(
        context="SBS1536",
        labels=tuple(labels),
        columns=("S1",),
        values=tuple((values.get(label, 0.0),) for label in labels),
    ), sum(values.get(label, 0.0) for label in labels)


def test_sbs1536_renderer_builds_only_a_default_variant_with_no_bias_panel():
    matrix, total = _sbs1536_matrix()
    figure = render_sample(matrix, "S1", unit="counts")
    assert figure["kind"] == "sbs4608"
    assert set(figure["variants"].keys()) == {"default"}
    built = figure["variants"]["default"]
    assert "traces" in built and "layout" in built
    # No transcription-bias traces are ever built for this context.
    bias_start, bias_end = built["sections"]["transcription_bias"]["trace_range"]
    assert bias_start == bias_end == len(built["traces"])
    # Nor is there a yaxis5/yaxis6/legend left behind for it.
    assert "yaxis5" not in built["layout"]
    assert "yaxis6" not in built["layout"]
    assert "legend" not in built["layout"]
    # The reserved bottom band is reclaimed, not left as blank padding.
    assert built["layout"]["height"] < 1200


def test_sbs1536_counts_unit_recovers_raw_totals():
    matrix, total = _sbs1536_matrix()
    built = render_sample(matrix, "S1", unit="counts")["variants"]["default"]
    bar_sum = sum(v for t in built["traces"] if t.get("type") == "bar" for v in t["y"])
    assert bar_sum == pytest.approx(total)
    heatmap_sum = sum(
        cell
        for t in built["traces"]
        if t.get("type") == "heatmap" and t.get("yaxis") == "y4"
        for row in t["z"]
        for cell in row
    )
    assert heatmap_sum == pytest.approx(total)
    assert built["layout"]["yaxis"]["title"]["text"] == "Number of Single Base Substitutions"


def test_sbs1536_total_heatmap_hover_has_no_fabricated_strand_breakdown():
    matrix, _total = _sbs1536_matrix()
    built = render_sample(matrix, "S1", unit="percent")["variants"]["default"]
    total_heatmap_traces = [
        t for t in built["traces"] if t.get("type") == "heatmap" and t.get("yaxis") == "y4"
    ]
    assert total_heatmap_traces
    for trace in total_heatmap_traces:
        assert "T:" not in trace["hovertemplate"]
        assert "customdata.T" not in trace["hovertemplate"]


def _uniform_sbs4608_matrix(value=1.0):
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
        values=tuple((value,) for _ in labels),
    )
    return matrix, float(len(labels)) * value


def _main_bar_sum(built):
    return sum(
        v
        for t in built["traces"]
        if t.get("type") == "bar" and t.get("yaxis") in (None, "y")
        for v in t["y"]
    )


def _total_heatmap_z_sum(built):
    return sum(
        cell
        for t in built["traces"]
        if t.get("type") == "heatmap" and t.get("yaxis") == "y4"
        for row in t["z"]
        for cell in row
    )


# The main plot's own bar sum only covers everything for "default" and
# "genic_bias" (Genic + Intergenic = T + U + N = every row); "strand_bias"
# plots only the T and U bars, so it covers 2/3 of a uniform matrix's total.
MAIN_BAR_FRACTION = {"default": 1.0, "strand_bias": 2.0 / 3.0, "genic_bias": 1.0}


@pytest.mark.parametrize("variant", ["default", "strand_bias", "genic_bias"])
def test_sbs4608_counts_unit_recovers_raw_totals(variant):
    # Every row contributes the same raw count, so the bar chart / total
    # heatmap -- a full partition of all 4608 rows either way -- must sum
    # back to that exact total once "counts" undoes the percent-of-total
    # SBS4608Data/SBS4608PlotBuilder always normalize to internally.
    matrix, total = _uniform_sbs4608_matrix()
    built = render_sample(matrix, "S1", unit="counts")["variants"][variant]
    assert _main_bar_sum(built) == pytest.approx(total * MAIN_BAR_FRACTION[variant])
    assert _total_heatmap_z_sum(built) == pytest.approx(total)
    yaxis = built["layout"]["yaxis"]
    assert yaxis["title"]["text"] == "Number of Single Base Substitutions"
    assert "ticksuffix" not in yaxis


@pytest.mark.parametrize("variant", ["default", "strand_bias", "genic_bias"])
def test_sbs4608_proportion_unit_sums_to_one(variant):
    matrix, _total = _uniform_sbs4608_matrix()
    built = render_sample(matrix, "S1", unit="proportion")["variants"][variant]
    assert _main_bar_sum(built) == pytest.approx(MAIN_BAR_FRACTION[variant])
    assert _total_heatmap_z_sum(built) == pytest.approx(1.0)
    yaxis = built["layout"]["yaxis"]
    assert yaxis["title"]["text"] == "Proportion of Single Base Substitutions"
    assert "ticksuffix" not in yaxis


def test_sbs4608_percent_unit_is_unchanged_by_the_new_scaling_path():
    matrix, _total = _uniform_sbs4608_matrix()
    built = render_sample(matrix, "S1", unit="percent")["variants"]["default"]
    assert _main_bar_sum(built) == pytest.approx(100.0)
    assert built["layout"]["yaxis"]["title"]["text"] == "Percentage of Single Base Substitutions"
    assert built["layout"]["yaxis"]["ticksuffix"] == "%"


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

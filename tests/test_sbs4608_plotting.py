from pathlib import Path

import pytest

import sigProfilerPlotting as plotting_package
from sigProfilerPlotting import plotInteractive, plotStatic
from sigProfilerPlotting.sbs4608.sbs4608_data import (
    DINUCS,
    SUBS,
    load_raw4608,
)


def _write_sbs4608_matrix(
    path,
    columns=("Sample",),
    all_zero=False,
    omit_last=False,
    contributions=None,
):
    contributions = contributions or {}
    rows = []
    for strand in "TUN":
        for substitution in SUBS:
            for left in DINUCS:
                for right in DINUCS:
                    mutation_type = f"{strand}:{left}[{substitution}]{right}"
                    value = (
                        0
                        if all_zero
                        else contributions.get(
                            mutation_type, int(mutation_type == "T:AA[C>A]AA")
                        )
                    )
                    rows.append((mutation_type, *([value] * len(columns))))
    if omit_last:
        rows.pop()

    with open(path, "w", encoding="utf-8") as output:
        output.write("MutationType\t" + "\t".join(columns) + "\n")
        for row in rows:
            output.write("\t".join(map(str, row)) + "\n")
    return path


def test_plot_functions_are_public():
    assert plotting_package.plotInteractive is plotInteractive
    assert plotting_package.plotStatic is plotStatic


def test_load_raw4608_requires_a_complete_matrix(tmp_path):
    matrix = _write_sbs4608_matrix(tmp_path / "incomplete.tsv", omit_last=True)

    with pytest.raises(ValueError, match="found 4607 of 4608"):
        load_raw4608(matrix, "Sample")


@pytest.mark.parametrize("invalid_value", ["nan", "inf", "-1"])
def test_load_raw4608_rejects_invalid_contributions(tmp_path, invalid_value):
    matrix = _write_sbs4608_matrix(tmp_path / "invalid.tsv")
    lines = matrix.read_text(encoding="utf-8").splitlines()
    lines[1] = lines[1].rsplit("\t", 1)[0] + f"\t{invalid_value}"
    matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="finite and non-negative"):
        load_raw4608(matrix, "Sample")


def test_load_raw4608_rejects_duplicate_contexts(tmp_path):
    matrix = _write_sbs4608_matrix(tmp_path / "duplicate.tsv")
    lines = matrix.read_text(encoding="utf-8").splitlines()
    lines[2] = lines[1]
    matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Duplicate SBS4608 MutationType"):
        load_raw4608(matrix, "Sample")


def test_plot_interactive_handles_zero_counts_and_escapes_labels(tmp_path):
    matrix = _write_sbs4608_matrix(
        tmp_path / "zero.tsv", columns=("Sample <zero>",), all_zero=True
    )

    output = Path(
        plotInteractive(
            matrix,
            "unsafe/project <zero>",
            output_path=tmp_path / "plots",
            include_detailed_heatmaps=False,
        )
    )

    assert output.name == "unsafe_project_zero_interactive.html"
    html = output.read_text(encoding="utf-8")
    assert "Sample &lt;zero&gt;" in html
    assert "unsafe/project &lt;zero&gt; -- interactive SBS-4608 plots" in html
    assert 'id="main-plot-type-select"' in html
    assert "applyPlotOptions" in html
    assert '<script src="https://cdn.plot.ly/' not in html
    assert "plotly.js v" in html


def test_count_matrices_are_normalized_to_percentages(tmp_path):
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig

    matrix = _write_sbs4608_matrix(
        tmp_path / "counts.tsv",
        contributions={
            "T:AA[C>A]AA": 400,
            "U:AA[C>A]AA": 200,
            "N:AA[C>A]AA": 200,
        },
    )
    raw = load_raw4608(matrix, "Sample")

    strand_data = SBS4608Data(raw, "Sample", main_plot_type="strand_bias")
    assert sum(row["mutations"] for row in strand_data.get_main_plot_data()) == 75

    genic_data = SBS4608Data(raw, "Sample", main_plot_type="genic_bias")
    assert sum(row["mutations"] for row in genic_data.get_main_plot_data()) == 100

    figure = SBS4608PlotBuilder(
        SBS4608Data(raw, "Sample"),
        SBS4608PlotConfig(include_detailed_heatmaps=True),
    ).build()
    heatmap_values = [
        value
        for trace in figure["traces"]
        if trace.get("type") == "heatmap"
        for row in trace["z"]
        for value in row
    ]
    assert max(heatmap_values) == 100


def test_plot_static_builds_a_valid_multipage_container(tmp_path, monkeypatch):
    import plotly.graph_objects as go
    from pypdf import PdfReader, PdfWriter

    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")

    def fake_write_image(_figure, path, format):
        assert format == "pdf"
        writer = PdfWriter()
        writer.add_blank_page(width=1080, height=1200)
        with open(path, "wb") as output:
            writer.write(output)
        writer.close()

    monkeypatch.setattr(go.Figure, "write_image", fake_write_image)
    output = Path(
        plotStatic(
            matrix,
            "fixture",
            output_path=tmp_path / "plots",
            include_detailed_heatmaps=False,
        )
    )

    assert output.name == "fixture_static.pdf"
    assert len(PdfReader(output).pages) == 1

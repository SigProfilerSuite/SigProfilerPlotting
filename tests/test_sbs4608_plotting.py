from pathlib import Path
import json
import os
import shutil
import subprocess

import pytest

import sigProfilerPlotting as plotting_package
from sigProfilerPlotting import plotInteractive, plotSBS
from sigProfilerPlotting.sbs4608.sbs4608_data import (
    DINUCS,
    SUBS,
    load_raw4608,
)


def _mock_static_page_export(monkeypatch, capture=None):
    """Keep API tests fast while exercising the real multi-page container."""
    from pypdf import PdfWriter
    import sigProfilerPlotting.sbs4608_static as static_module

    def fake_render(figure, _path):
        if capture is not None:
            capture(figure)

    def fake_assemble(_png_path, pdf_path):
        writer = PdfWriter()
        writer.add_blank_page(width=810, height=900)
        with open(pdf_path, "wb") as output:
            writer.write(output)
        writer.close()

    monkeypatch.setattr(static_module, "_render_sbs4608_png", fake_render)
    monkeypatch.setattr(static_module, "_png_to_pdf_page", fake_assemble)


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


def test_heatmap_colorbar_does_not_double_scale_percentages(tmp_path):
    # TraceFormatter.apply_percentage_formatting() already multiplies every
    # heatmap z value by 100 (e.g. a fraction of 0.093 becomes 9.3, meaning
    # "9.3%") before configure_heatmap_colorbar() runs. A colorbar
    # tickformat of ".1%" would multiply by 100 *again* (9.3 -> "930%"),
    # which is exactly the "weird" >100% values a user reported seeing.
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig

    path = _write_sbs4608_matrix(tmp_path / "sbs4608.tsv")
    raw = load_raw4608(str(path), "Sample")
    config = SBS4608PlotConfig(include_detailed_heatmaps=True, main_plot_type="default")
    built = SBS4608PlotBuilder(SBS4608Data(raw, "Sample", "default"), config).build()

    heatmaps_with_colorbars = [
        t for t in built["traces"] if t.get("type") == "heatmap" and t.get("showscale")
    ]
    assert heatmaps_with_colorbars, "expected at least one heatmap trace with a visible colorbar"
    for trace in heatmaps_with_colorbars:
        colorbar = trace["colorbar"]
        assert "tickformat" not in colorbar, (
            "a Plotly '%' tickformat (e.g. '.1%') would re-multiply by 100 on top of "
            "apply_percentage_formatting()'s own *100 -- explicit ticktext is used instead"
        )
        for text in colorbar["ticktext"]:
            assert text.endswith("%")
            float(text[:-1])  # must parse as a plain number, not raise


def test_heatmap_colorbar_labels_every_color_band_boundary(tmp_path):
    # The heatmap colorscale (sbs4608_plot_config.py) is five discrete bands
    # at the 0/0.2/0.4/0.6/0.8/1.0 fractions of [zmin, zmax]. Each boundary
    # should carry its own percentage label instead of Plotly's autotick
    # picking arbitrary round numbers unrelated to the band edges.
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig

    path = _write_sbs4608_matrix(tmp_path / "sbs4608.tsv")
    raw = load_raw4608(str(path), "Sample")
    config = SBS4608PlotConfig(include_detailed_heatmaps=True, main_plot_type="default")
    built = SBS4608PlotBuilder(SBS4608Data(raw, "Sample", "default"), config).build()

    heatmaps_with_colorbars = [
        t for t in built["traces"] if t.get("type") == "heatmap" and t.get("showscale")
    ]
    assert heatmaps_with_colorbars

    for trace in heatmaps_with_colorbars:
        colorbar = trace["colorbar"]
        assert colorbar["tickmode"] == "array"
        assert len(colorbar["tickvals"]) == 6
        assert len(colorbar["ticktext"]) == 6

        zmin, zmax = trace["zmin"], trace["zmax"]
        expected_vals = [zmin + frac * (zmax - zmin) for frac in (0, 0.2, 0.4, 0.6, 0.8, 1.0)]
        for actual, expected in zip(colorbar["tickvals"], expected_vals):
            assert actual == pytest.approx(expected)
        for text, expected in zip(colorbar["ticktext"], expected_vals):
            assert text == f"{expected:.1f}%"


def test_plot_functions_are_public():
    assert plotting_package.plotInteractive is plotInteractive
    assert plotting_package.plotSBS is plotSBS
    # plotStatic was folded into plotSBS(..., "4608", ...) -- there must be
    # no separate public symbol left behind.
    assert not hasattr(plotting_package, "plotStatic")


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


@pytest.mark.parametrize("include_details", [False, True])
def test_interactive_options_compact_and_restore_subplot_layout(
    tmp_path, include_details
):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to exercise client-side layout calculations")
    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")
    output = Path(
        plotInteractive(
            matrix,
            "compact",
            output_path=tmp_path,
            include_detailed_heatmaps=include_details,
        )
    )
    html = output.read_text(encoding="utf-8")
    registry = json.JSONDecoder().raw_decode(html.split("const PLOTS = ", 1)[1])[0]
    function = html.split("function compactPlot(", 1)[1].split(
        "function applyPlotOptions", 1
    )[0]
    script = "const assert = require('node:assert/strict');\n"
    script += "const SECTION_AXES = {total_heatmap:['yaxis4'], detailed_heatmaps:['yaxis2','yaxis3'], transcription_bias:['yaxis5','yaxis6','xaxis2','xaxis3']};\n"
    script += "function compactPlot(" + function
    script += "const variants = " + json.dumps(registry[0]["variants"]) + ";\n"
    script += """
    for (const variant of Object.values(variants)) {
      const original = JSON.stringify(variant);
      const enabled = {total_heatmap:true, detailed_heatmaps:true, transcription_bias:true};
      for (let bits = 0; bits < 8; bits++) {
        const options = {total_heatmap:!!(bits & 1), detailed_heatmaps:!!(bits & 2), transcription_bias:!!(bits & 4)};
        const result = compactPlot(variant, options);
        assert.ok(result.layout.height <= variant.layout.height);
        for (const [name, sec] of Object.entries(variant.sections)) {
          for (let i = sec.trace_range[0]; i < sec.trace_range[1]; i++)
            assert.equal(result.traces[i].visible, options[name]);
        }
        const expectedMainPixels = (variant.layout.height - 180) * (variant.layout.yaxis.domain[1] - variant.layout.yaxis.domain[0]);
        const actualMainPixels = (result.layout.height - 180) * (result.layout.yaxis.domain[1] - result.layout.yaxis.domain[0]);
        assert.ok(Math.abs(expectedMainPixels - actualMainPixels) < 1);
        if (bits === 0) assert.ok(result.layout.height < variant.layout.height * 0.75);
        if (!options.transcription_bias && variant.main_plot_type !== 'default') {
          assert.equal(result.layout.legend.orientation, 'h');
          assert.equal(result.layout.legend.x, 0.5);
          assert.ok(result.layout.legend.y >= 0);
        }
        if (!options.transcription_bias && variant.main_plot_type === 'default')
          assert.equal(result.layout.showlegend, false);
        const restored = compactPlot(variant, enabled);
        assert.equal(restored.layout.height, variant.layout.height);
        assert.deepEqual(restored.layout.yaxis.domain, variant.layout.yaxis.domain);
        assert.equal(JSON.stringify(variant), original);
      }
    }
    """
    # The full-resolution compact/restore registry (3 variants, including
    # detailed heatmaps) serialized into `script` can run into the low
    # hundreds of KB -- past `node -e`'s command-line argument on some CI
    # runners ("OSError: [Errno 7] Argument list too long"). A script file
    # has no such limit.
    script_path = tmp_path / "compact_plot_check.js"
    script_path.write_text(script, encoding="utf-8")
    subprocess.run([node, str(script_path)], check=True, capture_output=True, text=True)


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


def test_plotSBS_4608_builds_a_valid_multipage_container(tmp_path, monkeypatch):
    from pypdf import PdfReader

    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")

    _mock_static_page_export(monkeypatch)
    output = Path(
        plotSBS(str(matrix), str(tmp_path / "plots"), "fixture", "4608")
    )

    assert output.name == "SBS_4608_plots_fixture.pdf"
    page = PdfReader(output).pages[0]
    assert float(page.mediabox.width) == 810
    assert float(page.mediabox.height) == 900


def test_plotSBS_4608_plots_only_the_requested_signatures(tmp_path, monkeypatch):
    from pypdf import PdfReader

    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv", columns=("S1", "S2", "S3"))

    _mock_static_page_export(monkeypatch)
    output = Path(
        plotSBS(
            str(matrix), str(tmp_path / "plots"), "subset", "4608",
            signatures=["S1", "S3"],
        )
    )
    assert len(PdfReader(output).pages) == 2


def test_plotSBS_4608_rejects_an_unknown_signature_name(tmp_path):
    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")
    with pytest.raises(ValueError, match="Signature.*not found"):
        plotSBS(
            str(matrix), str(tmp_path / "plots"), "fixture", "4608",
            signatures=["does-not-exist"],
        )


def test_plotSBS_4608_rejects_non_pdf_format(tmp_path):
    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")
    with pytest.raises(ValueError, match=r"supports only savefig_format='pdf'"):
        plotSBS(
            str(matrix), str(tmp_path / "plots"), "fixture", "4608",
            savefig_format="png",
        )


def test_plotSBS_4608_percentage_defaults_to_true(tmp_path, monkeypatch):
    # SBS4608 input is a signature/profile file, not raw mutation counts --
    # unlike every other plotSBS context, percentage should default on
    # without the caller passing it explicitly.
    matrix = _write_sbs4608_matrix(
        tmp_path / "matrix.tsv",
        contributions={"T:AA[C>A]AA": 400, "U:AA[C>A]AA": 200, "N:AA[C>A]AA": 200},
    )
    captured = {}

    def capture(figure):
        captured["heatmap_max"] = max(
            (value for trace in figure.data if trace.type == "heatmap" for row in trace.z for value in row),
            default=0,
        )

    _mock_static_page_export(monkeypatch, capture)
    plotSBS(str(matrix), str(tmp_path / "plots"), "defaultpct", "4608")  # percentage omitted
    assert captured["heatmap_max"] == pytest.approx(100.0)


def test_plotSBS_4608_percentage_false_recovers_raw_contributions(tmp_path, monkeypatch):
    matrix = _write_sbs4608_matrix(
        tmp_path / "matrix.tsv",
        contributions={"T:AA[C>A]AA": 400, "U:AA[C>A]AA": 200, "N:AA[C>A]AA": 200},
    )
    captured = {}

    def capture(figure):
        captured["heatmap_max"] = max(
            (value for trace in figure.data if trace.type == "heatmap" for row in trace.z for value in row),
            default=0,
        )
        captured["ytitle"] = figure.layout.yaxis.title.text
        captured["colorbars"] = [
            list(trace.colorbar.ticktext)
            for trace in figure.data
            if trace.type == "heatmap" and trace.showscale
        ]

    _mock_static_page_export(monkeypatch, capture)
    plotSBS(str(matrix), str(tmp_path / "plots"), "raw", "4608", percentage=False)
    # The "total" heatmap merges all three strand tags for one pentanucleotide
    # context, so the AA[C>A]AA cell's raw contribution is 400+200+200=800.
    assert captured["heatmap_max"] == pytest.approx(800.0)
    assert captured["ytitle"] == "Raw signature contribution"
    assert captured["colorbars"]
    assert all(text != "0" for texts in captured["colorbars"] for text in texts[1:])


def test_raw_signature_contribution_precision_for_profile_summing_to_one(tmp_path):
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig
    from sigProfilerPlotting.sbs4608_static import _scale_to_raw_contributions

    nonzero = {
        "T:AA[C>A]AA": 0.03125,
        "U:AA[C>A]AA": 0.015625,
        "N:AA[C>A]AA": 0.0078125,
        "T:AC[C>G]GT": 0.00003125,
    }
    # Keep a realistic profile total (~1) while retaining a tiny nonzero cell.
    nonzero["T:TT[T>G]TT"] = 1.0 - sum(nonzero.values())
    matrix = _write_sbs4608_matrix(tmp_path / "profile.tsv", contributions=nonzero)
    raw = load_raw4608(str(matrix), "Sample")
    built = SBS4608PlotBuilder(
        SBS4608Data(raw, "Sample"),
        SBS4608PlotConfig(include_detailed_heatmaps=True),
    ).build()
    _scale_to_raw_contributions(built, sum(row["contribution"] for row in raw))

    assert built["layout"]["yaxis"]["title"]["text"] == "Raw signature contribution"
    assert built["layout"]["yaxis"]["tickformat"] == ".3g"
    colorbars = [
        trace["colorbar"]
        for trace in built["traces"]
        if trace.get("type") == "heatmap" and trace.get("showscale")
    ]
    assert colorbars
    for colorbar in colorbars:
        for value, text in zip(colorbar["tickvals"], colorbar["ticktext"]):
            if value:
                assert text != "0"
                assert float(text) != 0


def test_plotSBS_4608_zero_valued_matrix_does_not_raise(tmp_path, monkeypatch):
    matrix = _write_sbs4608_matrix(tmp_path / "zero.tsv", all_zero=True)
    _mock_static_page_export(monkeypatch)
    plotSBS(str(matrix), str(tmp_path / "plots"), "zero", "4608")


@pytest.mark.parametrize("invalid_value", ["nan", "inf", "-1"])
def test_plotSBS_4608_rejects_invalid_contributions(tmp_path, invalid_value):
    matrix = _write_sbs4608_matrix(tmp_path / "invalid.tsv")
    lines = matrix.read_text(encoding="utf-8").splitlines()
    lines[1] = lines[1].rsplit("\t", 1)[0] + f"\t{invalid_value}"
    matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="finite and non-negative"):
        plotSBS(str(matrix), str(tmp_path / "plots"), "fixture", "4608")


def test_plotSBS_4608_incomplete_matrix_raises(tmp_path):
    matrix = _write_sbs4608_matrix(tmp_path / "incomplete.tsv", omit_last=True)
    with pytest.raises(ValueError, match="found 4607 of 4608"):
        plotSBS(str(matrix), str(tmp_path / "plots"), "fixture", "4608")


def test_plotSBS_4608_custom_text_is_rendered_as_an_annotation(tmp_path, monkeypatch):
    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")
    captured = {}

    def capture(figure):
        captured["texts"] = [a.text for a in figure.layout.annotations]

    _mock_static_page_export(monkeypatch, capture)
    plotSBS(
        str(matrix), str(tmp_path / "plots"), "fixture", "4608",
        custom_text_upper=["Stability: 0.97"],
        custom_text_middle=["Mutations: 12,345"],
    )
    assert any("Stability: 0.97" in t and "Mutations: 12,345" in t for t in captured["texts"])


def test_sbs4608_structural_guard_detects_a_missing_bias_panel(tmp_path):
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig
    from sigProfilerPlotting.sbs4608_static import _assert_sbs4608_figure_structure

    matrix = _write_sbs4608_matrix(tmp_path / "matrix.tsv")
    raw = load_raw4608(str(matrix), "Sample")
    built = SBS4608PlotBuilder(
        SBS4608Data(raw, "Sample"),
        SBS4608PlotConfig(include_detailed_heatmaps=True),
    ).build()
    _assert_sbs4608_figure_structure(built, "Sample")
    built["traces"] = [trace for trace in built["traces"] if trace.get("yaxis") != "y6"]
    with pytest.raises(RuntimeError, match="bias-summary"):
        _assert_sbs4608_figure_structure(built, "Sample")


def test_png_to_pdf_assembly_preserves_full_raster_and_page_ratio(tmp_path):
    from PIL import Image, ImageDraw
    from pypdf import PdfReader
    from sigProfilerPlotting.sbs4608_static import _png_to_pdf_page

    png_path = tmp_path / "page.png"
    pdf_path = tmp_path / "page.pdf"
    image = Image.new("RGB", (2160, 2400), "white")
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, 120, 120), fill="red")
    draw.rectangle((2040, 2280, 2159, 2399), fill="blue")
    image.save(png_path, "PNG")

    _png_to_pdf_page(png_path, pdf_path)
    page = PdfReader(pdf_path).pages[0]
    assert float(page.mediabox.width) == pytest.approx(810)
    assert float(page.mediabox.height) == pytest.approx(900)
    assert page.images
    assert page.images[0].image.size == (2160, 2400)


@pytest.mark.skipif(
    os.environ.get("SPP_RUN_KALEIDO_INTEGRATION") != "1",
    reason="set SPP_RUN_KALEIDO_INTEGRATION=1 for real Chrome/Kaleido release QA",
)
def test_real_kaleido_sbs4608_raster_pdf_integration(tmp_path):
    """Opt-in renderer test; normal unit tests do not require Chrome."""
    from pypdf import PdfReader

    contributions = {}
    index = 1
    total = sum(range(1, 4609))
    for strand in "TUN":
        for substitution in SUBS:
            for left in DINUCS:
                for right in DINUCS:
                    contributions[f"{strand}:{left}[{substitution}]{right}"] = index / total
                    index += 1
    matrix = _write_sbs4608_matrix(tmp_path / "profile.tsv", contributions=contributions)
    output = Path(plotSBS(str(matrix), str(tmp_path), "renderer", "4608"))
    page = PdfReader(output).pages[0]
    assert float(page.mediabox.width) == pytest.approx(810)
    assert float(page.mediabox.height) == pytest.approx(900)
    assert page.images and page.images[0].image.size == (2160, 2400)

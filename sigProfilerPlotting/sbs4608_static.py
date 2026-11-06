"""Static (Plotly -> PDF) SBS-4608 matrix and signature plotting.

Builds a static PDF report with one multi-panel figure per signature: the
96-context bar profile, "front"/"back" pentanucleotide detail heatmaps, a
total heatmap, and stacked Transcribed-vs-Untranscribed + Genic-vs-Intergenic
strand-bias bars -- the same panels as plotInteractive, rendered to one PDF
page per signature instead of an interactive HTML page.

There is no public function in this module. The only public entry point is
``sigPlt.plotSBS(matrix_path, output_path, project, "4608", ...)``, which
calls the private ``_build_sbs4608_pdf`` below. A previous version of this
module exposed a separate public ``plotStatic()`` function with its own
argument order/conventions; it has been folded into ``plotSBS`` so SBS4608
has exactly one static entry point, consistent with every other context.

Unlike plotInteractive(...), which embeds all three main-plot-type variants
(Default / Strand Bias / Genic Bias) in one page with a live dropdown to
switch between them client side, a static PDF page has no client-side
JavaScript to switch views with, so it always renders the "default" view
(the 96-context bar profile) with detailed heatmaps included -- the same
fixed choice ``plotStatic()`` used before. Use this when the plot needs to
go into a manuscript, a print handout, or anywhere an HTML file with
embedded JavaScript isn't practical; use plotInteractive when the viewer
needs hover tooltips or the Main Plot Type / plot-options controls.

This function expects a strand-tagged pentanucleotide (SBS-4608)
signatures/matrix file -- mutationType rows like "T:AA[C>A]AA" (4608 rows
per signature: 6 substitutions x 16 x 16 flanking dinucleotide combos x 3
strand tags T/U/N) -- the same format plotInteractive expects. Unlike every
other plotSBS context, it is file-path only (not DataFrame), since the
underlying SBS4608Data/SBS4608PlotBuilder machinery it shares with
plotInteractive reads directly from a *_Signatures.txt-style file.

Requires the `kaleido` package (for static rendering) and `pypdf` (to merge
per-signature pages into one file), plus a local Chrome/Chromium install.
kaleido finds an existing Chrome/Chromium automatically on most machines; if
it can't, run `plotly_get_chrome` once (works on macOS, Windows, and Linux
x86-64).

Example:
    import sigProfilerPlotting as sigPlt
    sigPlt.plotSBS(
        "SBS4608_S5_Signatures.txt", "plots", "S5", "4608", percentage=True,
    )
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Optional

from PIL import Image, ImageChops

from .sbs4608.sbs4608_data import SBS4608Data, load_raw4608
from .sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
from .sbs4608.sbs4608_plot_config import SBS4608PlotConfig

# Plotly.js silently falls back to its schema default ("auto") for an
# invalid annotation xanchor/yanchor value, which is how plotInteractive's
# HTML renders fine even though some annotations use "bottom"/"top" for
# xanchor (not a valid enum member). The Python graph_objects validator used
# for static rendering is strict, so we normalize before constructing a
# Figure -- this reproduces the same fallback the browser already applies,
# rather than changing how anything looks.
_VALID_ANCHORS = ("auto", "left", "center", "right")

_CUSTOM_TEXT_MAX_LENGTH = 40

# Direct Kaleido vector PDFs can contain malformed clipping state that PDF
# engines interpret differently. Render a complete high-resolution raster
# and wrap it in a simple PDF page instead. At scale 2, 2160x2400 px maps
# exactly to the historical 810x900 pt page at 192 dpi.
_PLOT_WIDTH_PX = 1080
_PLOT_HEIGHT_PX = 1200
_PNG_SCALE = 2
_PDF_WIDTH_PT = 810
_PDF_HEIGHT_PT = 900
_EXPECTED_MUTATION_HEADERS = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


def _safe_filename(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return safe or "SBS4608"


def _custom_text_entry(texts, index):
    """One sample's custom-text string, or None -- applying the same
    <=40-character limit (print a warning and drop it) every other plotSBS
    context's custom_text_upper/middle/bottom already uses, for consistency.
    """
    if not texts or index >= len(texts):
        return None
    value = texts[index]
    if value and len(value) > _CUSTOM_TEXT_MAX_LENGTH:
        print(
            "To add a custom text, please limit the string to <40 characters including spaces."
        )
        return None
    return value


def _custom_text_annotation(custom_text_upper, custom_text_middle, custom_text_bottom, index):
    pieces = [
        text
        for text in (
            _custom_text_entry(custom_text_upper, index),
            _custom_text_entry(custom_text_middle, index),
            _custom_text_entry(custom_text_bottom, index),
        )
        if text
    ]
    if not pieces:
        return []
    return [{
        "text": "<br>".join(pieces),
        "xref": "paper", "yref": "paper",
        "x": 0.99, "y": 0.99,
        "xanchor": "right", "yanchor": "top",
        "showarrow": False,
        "align": "right",
        "font": {"size": 14, "family": "Arial", "color": "black"},
    }]


def _plain_annotation_text(value):
    return re.sub(r"<[^>]+>", "", str(value or "")).strip()


def _assert_sbs4608_figure_structure(built: dict, signature: str) -> None:
    """Require every panel and label that defines a complete static page."""
    traces = built.get("traces", [])
    bars = [
        trace for trace in traces
        if trace.get("type") == "bar" and trace.get("yaxis") in (None, "y")
    ]
    if not set(_EXPECTED_MUTATION_HEADERS).issubset(
        {trace.get("name") for trace in bars}
    ):
        raise RuntimeError("Incomplete SBS4608 main profile: mutation-class bars are missing.")

    heatmap_axes = {
        trace.get("yaxis") for trace in traces if trace.get("type") == "heatmap"
    }
    if not {"y2", "y3", "y4"}.issubset(heatmap_axes):
        raise RuntimeError("Incomplete SBS4608 page: front/back/total heatmaps are required.")

    bias_axes = {
        trace.get("yaxis") for trace in traces if trace.get("type") == "scatter"
    }
    has_signal = any(float(value) != 0 for trace in bars for value in (trace.get("y") or []))
    layout = built.get("layout", {})
    bias_rows_present = {"yaxis5", "yaxis6"}.issubset(layout)
    if (has_signal and not {"y5", "y6"}.issubset(bias_axes)) or not bias_rows_present:
        raise RuntimeError("Incomplete SBS4608 page: both bias-summary rows are required.")

    annotation_text = {
        _plain_annotation_text(annotation.get("text"))
        for annotation in layout.get("annotations", [])
    }
    required = {
        *_EXPECTED_MUTATION_HEADERS,
        signature,
        "Transcribed / Untranscribed strand",
        "Genic / Intergenic",
    }
    missing = sorted(required - annotation_text)
    if missing:
        raise RuntimeError(f"Incomplete SBS4608 page annotations: missing {missing}.")


def _format_raw_contribution(value: float) -> str:
    """Adaptive precision that never formats a nonzero contribution as zero."""
    value = float(value)
    return "0" if value == 0 else f"{value:.3g}"


def _scale_to_raw_contributions(built: dict, total: float) -> None:
    """Undo percent normalization without mislabeling profile values as counts."""
    factor = total / 100.0 if total > 0 else 0.0
    for trace in built["traces"]:
        if trace.get("type") == "heatmap":
            if trace.get("z"):
                trace["z"] = [[cell * factor for cell in row] for row in trace["z"]]
            if "zmin" in trace:
                trace["zmin"] *= factor
            if "zmax" in trace:
                trace["zmax"] *= factor
            colorbar = trace.get("colorbar")
            if colorbar and colorbar.get("tickmode") == "array":
                colorbar["tickvals"] = [value * factor for value in colorbar["tickvals"]]
                colorbar["ticktext"] = [
                    _format_raw_contribution(value) for value in colorbar["tickvals"]
                ]
            trace["hovertemplate"] = (
                "x: %{x}<br>y: %{y}<br>Raw contribution: %{z:.4g}<extra></extra>"
            )
        elif trace.get("type") == "bar" and trace.get("yaxis") in (None, "y"):
            if trace.get("y"):
                trace["y"] = [value * factor for value in trace["y"]]
            trace["hovertemplate"] = "%{x}<br>Raw contribution: %{y:.4g}<extra></extra>"

    yaxis = built["layout"].get("yaxis")
    if yaxis:
        yaxis["title"] = {
            **yaxis.get("title", {}),
            "text": "Raw signature contribution",
        }
        yaxis.pop("ticksuffix", None)
        yaxis["tickformat"] = ".3g"
        if "range" in yaxis:
            yaxis["range"] = [value * factor for value in yaxis["range"]]


def _validate_rendered_png(path: Path) -> None:
    """Reject a truncated raster before it can enter the report."""
    expected_size = (_PLOT_WIDTH_PX * _PNG_SCALE, _PLOT_HEIGHT_PX * _PNG_SCALE)
    with Image.open(path) as image:
        image.load()
        if image.format != "PNG" or image.size != expected_size:
            raise RuntimeError(
                f"Incomplete SBS4608 raster export: expected PNG {expected_size}, "
                f"got {image.format} {image.size}."
            )
        rgb = image.convert("RGB")
        width, height = rgb.size
        bands = {
            "main profile/header": (0, 0, width, round(height * 0.36)),
            "heatmaps": (0, round(height * 0.32), width, round(height * 0.79)),
            "bias summaries": (0, round(height * 0.77), width, height),
        }
        for name, box in bands.items():
            crop = rgb.crop(box)
            diff = ImageChops.difference(crop, Image.new("RGB", crop.size, "white"))
            mask = diff.convert("L").point(lambda value: 255 if value > 8 else 0)
            if mask.histogram()[255] < 500:
                raise RuntimeError(
                    f"Incomplete SBS4608 raster export: {name} band is blank/clipped."
                )


def _render_sbs4608_png(figure, path: Path) -> None:
    figure.write_image(
        str(path),
        format="png",
        width=_PLOT_WIDTH_PX,
        height=_PLOT_HEIGHT_PX,
        scale=_PNG_SCALE,
    )
    _validate_rendered_png(path)


def _png_to_pdf_page(png_path: Path, pdf_path: Path) -> None:
    """Wrap one complete raster in an exact 810x900 pt PDF page."""
    expected_size = (_PLOT_WIDTH_PX * _PNG_SCALE, _PLOT_HEIGHT_PX * _PNG_SCALE)
    with Image.open(png_path) as image:
        image.load()
        if image.size != expected_size:
            raise RuntimeError(
                f"Cannot assemble SBS4608 PDF page from {image.size}; expected {expected_size}."
            )
        resolution = image.width * 72.0 / _PDF_WIDTH_PT
        image.convert("RGB").save(
            pdf_path,
            "PDF",
            resolution=resolution,
            quality=95,
            subsampling=0,
        )

    from pypdf import PdfReader

    reader = PdfReader(str(pdf_path))
    if len(reader.pages) != 1:
        raise RuntimeError("SBS4608 raster-to-PDF assembly did not produce one page.")
    page = reader.pages[0]
    width, height = float(page.mediabox.width), float(page.mediabox.height)
    if abs(width - _PDF_WIDTH_PT) > 0.01 or abs(height - _PDF_HEIGHT_PT) > 0.01:
        raise RuntimeError(
            f"SBS4608 PDF page has size {width:g}x{height:g} pt; "
            f"expected {_PDF_WIDTH_PT}x{_PDF_HEIGHT_PT} pt."
        )
    resources = page.get("/Resources", {}).get_object()
    xobjects = resources.get("/XObject", {}).get_object()
    embedded_sizes = {
        (int(obj.get_object().get("/Width", 0)), int(obj.get_object().get("/Height", 0)))
        for obj in xobjects.values()
        if obj.get_object().get("/Subtype") == "/Image"
    }
    if expected_size not in embedded_sizes:
        raise RuntimeError(
            "SBS4608 PDF page does not contain the full-resolution rendered image."
        )


def _build_sbs4608_pdf(
    matrix_path: str,
    output_path: str,
    project: str,
    percentage: bool = True,
    signatures: Optional[list] = None,
    custom_text_upper: Optional[list] = None,
    custom_text_middle: Optional[list] = None,
    custom_text_bottom: Optional[list] = None,
) -> str:
    """Build the multi-panel SBS4608 PDF, one page per signature/sample.

    Called from sigProfilerPlotting.plotSBS(..., "4608", ...); not public.
    """
    try:
        import plotly.graph_objects as go
        from pypdf import PdfWriter
    except ImportError as error:  # pragma: no cover - dependencies are packaged
        raise RuntimeError(
            "Static SBS4608 plotting requires plotly, kaleido, and pypdf. "
            "Reinstall SigProfilerPlotting with its declared dependencies."
        ) from error

    matrix_file = Path(matrix_path).expanduser().resolve()
    if not matrix_file.is_file():
        raise ValueError(f"matrix_path does not exist: {matrix_file}")

    with open(matrix_file, encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
    available_columns = header[1:] if header and header[0] == "MutationType" else []
    sig_columns = available_columns if signatures is None else list(signatures)
    if not sig_columns:
        raise ValueError("At least one SBS4608 sample/signature column is required.")
    missing = [s for s in sig_columns if s not in available_columns]
    if missing:
        raise ValueError(
            f"Signature(s) not found in {matrix_file}: {missing}. "
            f"Available columns: {available_columns}"
        )

    out_dir = Path(output_path).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"SBS_4608_plots_{_safe_filename(project)}.pdf"

    writer = PdfWriter()
    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            for sig_index, sig in enumerate(sig_columns):
                raw4608 = load_raw4608(str(matrix_file), sig)
                data = SBS4608Data(raw4608, sig, main_plot_type="default")
                config = SBS4608PlotConfig(
                    include_detailed_heatmaps=True, main_plot_type="default",
                )
                builder = SBS4608PlotBuilder(data, config)
                fig_dict = builder.build()

                # The builder always produces percent-of-total. Static
                # signature files may request their original profile
                # contributions, which are not mutation counts.
                total = sum(float(row["contribution"]) for row in raw4608)
                if not percentage:
                    _scale_to_raw_contributions(fig_dict, total)

                _assert_sbs4608_figure_structure(fig_dict, sig)

                layout = dict(fig_dict["layout"])
                annotations = [dict(a) for a in layout.get("annotations", [])]
                for ann in annotations:
                    if ann.get("xanchor") not in _VALID_ANCHORS:
                        ann["xanchor"] = "auto"
                annotations.extend(
                    _custom_text_annotation(
                        custom_text_upper, custom_text_middle, custom_text_bottom, sig_index,
                    )
                )
                layout["annotations"] = annotations

                fig = go.Figure(data=fig_dict["traces"], layout=layout)
                png_path = Path(tmp_dir) / f"page-{sig_index}.png"
                page_path = Path(tmp_dir) / f"page-{sig_index}.pdf"
                try:
                    _render_sbs4608_png(fig, png_path)
                    _png_to_pdf_page(png_path, page_path)
                except Exception as error:
                    raise RuntimeError(
                        "Static SBS4608 PDF export failed while rendering or validating "
                        "a complete raster page. Kaleido requires a local Chrome/Chromium "
                        "installation; run 'plotly_get_chrome' if one is not already "
                        "available."
                    ) from error
                writer.append(str(page_path))

            with open(out_file, "wb") as f:
                writer.write(f)
    finally:
        writer.close()

    return str(out_file.resolve())

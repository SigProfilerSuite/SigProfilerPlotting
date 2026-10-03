"""Static (Plotly -> PDF) SBS-4608 matrix and signature plotting.

Builds a static PDF report with one multi-panel figure per signature: the
96-context bar profile, "front"/"back" pentanucleotide detail heatmaps, a
total heatmap, and stacked Transcribed-vs-Untranscribed + Genic-vs-Intergenic
strand-bias bars -- the same panels as plotInteractive, rendered to one PDF
page per signature instead of an interactive HTML page.

Unlike plotInteractive(...), which embeds all three main-plot-type variants
(Default / Strand Bias / Genic Bias) in one page with a live dropdown to
switch between them client side, this function renders a single fixed view
(main_plot_type) straight to PDF, since a static page has no client-side
JavaScript to switch views with. Use this when the plot needs to go into a
manuscript, a print handout, or anywhere an HTML file with embedded
JavaScript isn't practical; use plotInteractive when the viewer needs hover
tooltips or the Main Plot Type / plot-options controls.

Unlike plotSBS(matrix, output_path, project, "4608", ...), which draws a
single static bar-profile PDF, this function expects a strand-tagged
pentanucleotide (SBS-4608) signatures/matrix file -- mutationType rows like
"T:AA[C>A]AA" (4608 rows per signature: 6 substitutions x 16 x 16 flanking
dinucleotide combos x 3 strand tags T/U/N) -- and renders the richer,
multi-panel view, one page per signature.

Requires the `kaleido` package (for static rendering) and `pypdf` (to merge
per-signature pages into one file), plus a local Chrome/Chromium install.
kaleido finds an existing Chrome/Chromium automatically on most machines; if
it can't, run `plotly_get_chrome` once (works on macOS, Windows, and Linux
x86-64).

Example:
    import sigProfilerPlotting as sigPlt
    sigPlt.plotStatic("SBS4608_S5_Signatures.txt", "S5")
    sigPlt.plotStatic(
        "SBS4608_S5_Signatures.txt", "S5", main_plot_type="strand_bias",
    )
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Optional

from .sbs4608.sbs4608_data import SBS4608Data, load_raw4608
from .sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
from .sbs4608.sbs4608_plot_config import SBS4608PlotConfig

_VALID_MAIN_PLOT_TYPES = ("default", "strand_bias", "genic_bias")

# Plotly.js silently falls back to its schema default ("auto") for an
# invalid annotation xanchor/yanchor value, which is how plotInteractive's
# HTML renders fine even though some annotations use "bottom"/"top" for
# xanchor (not a valid enum member). The Python graph_objects validator used
# for static rendering is strict, so we normalize before constructing a
# Figure -- this reproduces the same fallback the browser already applies,
# rather than changing how anything looks.
_VALID_ANCHORS = ("auto", "left", "center", "right")


def _safe_filename(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return safe or "SBS4608"


def plotStatic(
    matrix_path: str,
    project: str,
    output_path: Optional[str] = None,
    signatures: Optional[list] = None,
    include_detailed_heatmaps: bool = True,
    main_plot_type: str = "default",
) -> str:
    """Build a static multi-panel PDF plot for SBS-4608 (pentanucleotide +
    strand) mutational signatures -- one page per signature.

    Args:
        matrix_path: Path to a *_Signatures.txt (or matrix) file whose
            mutationType rows are strand-tagged pentanucleotide contexts,
            e.g. "T:AA[C>A]AA" (4608 rows per signature). Same format as
            plotInteractive expects.
        project: Name for this run, used in the output filename and the
            per-page signature title.
        output_path: Directory to write the PDF into. Defaults to
            "<project>_plots" in the current working directory.
        signatures: Which signature columns to plot. Defaults to every
            column in the file other than MutationType.
        include_detailed_heatmaps: Include the front/back pentanucleotide
            detail heatmaps (in addition to the total heatmap).
        main_plot_type: Which view the top panel shows -- one of "default"
            (the 96-context bar profile), "strand_bias" (Transcribed vs.
            Untranscribed strand bars), or "genic_bias" (Genic vs.
            Intergenic bars). A static page can't offer the live dropdown
            plotInteractive has, so exactly one variant is rendered.

    Returns:
        The path to the written PDF file.
    """
    try:
        import plotly.graph_objects as go
        from pypdf import PdfWriter
    except ImportError as error:  # pragma: no cover - dependencies are packaged
        raise RuntimeError(
            "Static SBS4608 plotting requires plotly, kaleido, and pypdf. "
            "Reinstall SigProfilerPlotting with its declared dependencies."
        ) from error

    if main_plot_type not in _VALID_MAIN_PLOT_TYPES:
        raise ValueError(
            f"main_plot_type must be one of {_VALID_MAIN_PLOT_TYPES}, got {main_plot_type!r}"
        )

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

    if output_path is None:
        output_path = f"{project}_plots"
    out_dir = Path(output_path).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{_safe_filename(project)}_static.pdf"

    writer = PdfWriter()
    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            for sig_index, sig in enumerate(sig_columns):
                raw4608 = load_raw4608(str(matrix_file), sig)
                data = SBS4608Data(raw4608, sig, main_plot_type=main_plot_type)
                config = SBS4608PlotConfig(
                    include_detailed_heatmaps=include_detailed_heatmaps,
                    main_plot_type=main_plot_type,
                )
                builder = SBS4608PlotBuilder(data, config)
                fig_dict = builder.build()

                layout = dict(fig_dict["layout"])
                annotations = [dict(a) for a in layout.get("annotations", [])]
                for ann in annotations:
                    if ann.get("xanchor") not in _VALID_ANCHORS:
                        ann["xanchor"] = "auto"
                layout["annotations"] = annotations

                fig = go.Figure(data=fig_dict["traces"], layout=layout)
                page_path = Path(tmp_dir) / f"page-{sig_index}.pdf"
                try:
                    fig.write_image(str(page_path), format="pdf")
                except Exception as error:
                    raise RuntimeError(
                        "Static SBS4608 PDF export failed. Kaleido requires a local "
                        "Chrome/Chromium installation; run 'plotly_get_chrome' if one "
                        "is not already available."
                    ) from error
                writer.append(str(page_path))

            with open(out_file, "wb") as f:
                writer.write(f)
    finally:
        writer.close()

    return str(out_file.resolve())

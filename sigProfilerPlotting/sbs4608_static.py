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

    # Deferred: avoids a circular import (web.matrix imports back from the
    # top-level sigProfilerPlotting module, which is what calls in here).
    from .web.renderers import _scale_sbs4608_unit

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

                # SBS4608Data/SBS4608PlotBuilder always build a percent-of-
                # total figure; rescale to raw counts when percentage=False
                # rather than a second, divergent scaling implementation
                # (the same post-hoc rescale the dashboard's counts/percent/
                # proportion unit toggle already uses and tests).
                total = sum(float(row["contribution"]) for row in raw4608)
                _scale_sbs4608_unit(fig_dict, total, "percent" if percentage else "counts")

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

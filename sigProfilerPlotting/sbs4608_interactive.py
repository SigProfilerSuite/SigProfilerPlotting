"""Interactive (Plotly) SBS-4608 matrix and signature plotting.

Builds a single, self-contained HTML report with one multi-panel Plotly
figure per signature: the 96-context bar profile, "front"/"back"
pentanucleotide detail heatmaps, a total heatmap with
Transcribed/Untranscribed/intergenic hover proportions, and stacked
Transcribed-vs-Untranscribed + Genic-vs-Intergenic strand-bias bars.

The page includes a "Main Plot Type" dropdown (Default / Strand Bias /
Genic Bias) that swaps the top panel's chart for every signature at once,
and a "Plot options" panel with checkboxes to toggle the Transcription bias
bars, the 5-base-context total heatmap, and the detailed front/back
heatmaps on or off. Both controls work purely client side -- all three main
plot variants are precomputed and embedded up front, and switching between
them (or toggling a checkbox) just swaps/shows/hides the relevant Plotly
traces and axis labels (Plotly.react / restyle / relayout), no data is
re-fetched or regenerated.

Unlike plotSBS(matrix, output_path, project, "4608", ...), which draws a
single static bar-profile PDF, this function expects a strand-tagged
pentanucleotide (SBS-4608) signatures/matrix file -- mutationType rows like
"T:AA[C>A]AA" (4608 rows per signature: 6 substitutions x 16 x 16 flanking
dinucleotide combos x 3 strand tags T/U/N) -- and renders the richer,
interactive multi-panel view.

Example:
    import sigProfilerPlotting as sigPlt
    sigPlt.plotInteractive("SBS4608_S5_Signatures.txt", "S5")
    sigPlt.plotInteractive(
        "SBS4608_S5_Signatures.txt", "S5", main_plot_type="strand_bias",
    )
"""

from __future__ import annotations

import json
import re
from html import escape
from pathlib import Path
from typing import Optional

from plotly.offline import get_plotlyjs

from .sbs4608.sbs4608_data import SBS4608Data, load_raw4608
from .sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
from .sbs4608.sbs4608_plot_config import SBS_COLOR, SBS4608PlotConfig

_VALID_MAIN_PLOT_TYPES = ("default", "strand_bias", "genic_bias")
_MAIN_PLOT_TYPE_LABELS = {
    "default": "Default",
    "strand_bias": "Strand Bias",
    "genic_bias": "Genic Bias",
}


def _safe_filename(value: str) -> str:
    """Return a filesystem-safe filename component without changing labels."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return safe or "SBS4608"


def _json_for_script(value) -> str:
    """Serialize JSON without allowing user labels to terminate a script tag."""
    return (
        json.dumps(value, allow_nan=False)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


_SBS4608_PAGE_TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<script>{plotly_js}</script>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px 24px 60px; background: #fafafa; line-height: 1.5; }}
h1 {{ font-size: 1.3rem; }}
h2 {{ border-top: 2px solid #ddd; padding-top: 20px; margin-top: 50px; font-size: 1.05rem; }}
nav.toc {{
  position: sticky; top: 0; z-index: 10;
  background: #fff; border: 1px solid #ddd; border-radius: 6px;
  padding: 10px 14px; margin: 16px 0 24px; display: flex; gap: 10px; flex-wrap: wrap;
  align-items: center;
}}
nav.toc span.label {{ font-weight: 600; color: #555; font-size: 0.85rem; margin-right: 4px; }}
nav.toc a {{
  text-decoration: none; color: #1a4fa0; font-size: 0.85rem;
  padding: 4px 10px; border: 1px solid #cddaf0; border-radius: 999px; background: #f2f6ff;
}}
nav.toc a:hover {{ background: #e2ecff; }}
.legend {{
  display: flex; gap: 14px; flex-wrap: wrap; align-items: center;
  background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 10px 14px; margin-bottom: 24px;
}}
.legend .swatch {{ display: inline-flex; align-items: center; gap: 6px; font-size: 0.85rem; color: #333; }}
.legend .dot {{ width: 12px; height: 12px; border-radius: 3px; display: inline-block; border: 1px solid rgba(0,0,0,0.15); }}
.main-plot-type {{
  display: flex; gap: 10px; flex-wrap: wrap; align-items: center;
  background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 10px 14px; margin-bottom: 24px;
}}
.main-plot-type span.label {{ font-weight: 600; color: #555; font-size: 0.85rem; margin-right: 4px; }}
.main-plot-type select {{
  font-size: 0.85rem; padding: 5px 8px; border-radius: 4px; border: 1px solid #ccc; background: #fff; cursor: pointer;
}}
.plot-options {{
  display: flex; gap: 18px; flex-wrap: wrap; align-items: center;
  background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 10px 14px; margin-bottom: 24px;
}}
.plot-options span.label {{ font-weight: 600; color: #555; font-size: 0.85rem; margin-right: 4px; }}
.plot-options label {{ font-size: 0.85rem; color: #333; display: inline-flex; align-items: center; gap: 5px; cursor: pointer; }}
.plot {{ background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<nav class="toc">
<span class="label">Jump to signature:</span>
{toc_links}
</nav>
<div class="legend">
<span class="label" style="margin-right:0;">Substitution type colors:</span>
{legend_swatches}
</div>
{main_plot_type_html}
{plot_options_html}
{sections}
<script>
{plot_calls}

const PLOTS = {plot_registry};

const SECTION_AXIS_OVERRIDES = {{
  total_heatmap: on => ({{'yaxis4.showticklabels': on, 'yaxis4.showline': on}}),
  detailed_heatmaps: on => ({{
    'yaxis2.showticklabels': on, 'yaxis2.showline': on,
    'yaxis3.showticklabels': on, 'yaxis3.showline': on,
  }}),
  transcription_bias: on => ({{'xaxis2.showticklabels': on, 'xaxis3.showticklabels': on}}),
}};

function applyPlotOptions() {{
  const optBias = document.getElementById('opt-bias');
  const optTotal = document.getElementById('opt-total');
  const optDetailed = document.getElementById('opt-detailed');
  const toggles = {{
    transcription_bias: optBias ? optBias.checked : true,
    total_heatmap: optTotal ? optTotal.checked : true,
    detailed_heatmaps: optDetailed ? optDetailed.checked : true,
  }};

  PLOTS.forEach(p => {{
    const visible = new Array(p.numTraces).fill(true);
    let layoutUpdate = {{}};

    for (const [key, sec] of Object.entries(p.sections)) {{
      const on = toggles[key];
      const [start, end] = sec.trace_range;
      for (let i = start; i < end; i++) visible[i] = on;
      layoutUpdate = {{...layoutUpdate, ...SECTION_AXIS_OVERRIDES[key](on)}};
      (sec.annotation_indices || []).forEach(idx => {{
        layoutUpdate[`annotations[${{idx}}].opacity`] = on ? 1 : 0;
      }});
    }}

    Plotly.restyle(p.div, {{visible: visible}});
    Plotly.relayout(p.div, layoutUpdate);
  }});
}}

function applyMainPlotType() {{
  const select = document.getElementById('main-plot-type-select');
  const chosen = select ? select.value : 'default';

  PLOTS.forEach(p => {{
    const variant = p.variants[chosen];
    if (!variant) return;
    Plotly.react(p.div, variant.traces, variant.layout, variant.config);
    p.sections = variant.sections;
    p.numTraces = variant.traces.length;
  }});

  // Re-apply the current checkbox states, since react() resets trace visibility.
  applyPlotOptions();
}}
</script>
</body>
</html>
"""


def plotInteractive(
    matrix_path: str,
    project: str,
    output_path: Optional[str] = None,
    signatures: Optional[list] = None,
    include_detailed_heatmaps: bool = True,
    main_plot_type: str = "default",
) -> str:
    """Build an interactive multi-panel HTML plot for SBS-4608 (pentanucleotide
    + strand) mutational signatures.

    Args:
        matrix_path: Path to a *_Signatures.txt (or matrix) file whose
            mutationType rows are strand-tagged pentanucleotide contexts,
            e.g. "T:AA[C>A]AA" (4608 rows per signature: 6 substitutions x
            16 x 16 flanking-dinucleotide combos x 3 strand tags T/U/N).
            Not the same as a plain 96/288/384/1536-context matrix -- use
            plotSBS(..., "4608") for a static PDF of that instead.
        project: Name for this run, used in the output filename and page title.
        output_path: Directory to write the HTML into. Defaults to
            "<project>_plots" in the current working directory.
        signatures: Which signature columns to plot. Defaults to every
            column in the file other than MutationType.
        include_detailed_heatmaps: Include the front/back pentanucleotide
            detail heatmaps (in addition to the total heatmap). When True,
            the page's "Plot options" panel also gets a checkbox to toggle
            them on/off.
        main_plot_type: Which view the top panel shows when the page first
            loads -- one of "default" (the 96-context bar profile),
            "strand_bias" (Transcribed vs. Untranscribed strand bars), or
            "genic_bias" (Genic vs. Intergenic bars). All three are always
            computed and embedded in the page; a "Main Plot Type" dropdown
            lets the viewer switch between them live, for every signature
            in this report at once, without regenerating anything.

    Returns:
        The path to the written HTML file.
    """
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
    out_file = out_dir / f"{_safe_filename(project)}_interactive.html"

    sections = []
    plot_calls = []
    plot_registry = []
    section_ids = []
    for sig_index, sig in enumerate(sig_columns):
        raw4608 = load_raw4608(str(matrix_file), sig)
        section_id = f"signature-{sig_index}"
        div_id = f"plot-{sig_index}"
        section_ids.append(section_id)

        variants = {}
        for mpt in _VALID_MAIN_PLOT_TYPES:
            data = SBS4608Data(raw4608, sig, main_plot_type=mpt)
            config = SBS4608PlotConfig(
                include_detailed_heatmaps=include_detailed_heatmaps, main_plot_type=mpt
            )
            builder = SBS4608PlotBuilder(data, config)
            fig = builder.build()
            variants[mpt] = {
                "traces": fig["traces"],
                "layout": fig["layout"],
                "config": fig["config"],
                "sections": fig["sections"],
            }

        initial = variants[main_plot_type]
        sections.append(
            f'<h2 id="{section_id}">{escape(sig)}</h2>\n'
            f'<div class="plot" id="{div_id}"></div>'
        )
        plot_calls.append(
            f'Plotly.newPlot("{div_id}", {_json_for_script(initial["traces"])}, '
            f'{_json_for_script(initial["layout"])}, '
            f'{_json_for_script(initial["config"])});'
        )
        plot_registry.append(
            {
                "div": div_id,
                "numTraces": len(initial["traces"]),
                "sections": initial["sections"],
                "variants": variants,
            }
        )

    legend_swatches = "\n".join(
        f'<span class="swatch"><span class="dot" style="background:{color};"></span>{mutation}</span>'
        for mutation, color in SBS_COLOR.items()
    )
    toc_links = "\n".join(
        f'<a href="#{section_id}">{escape(sig)}</a>'
        for section_id, sig in zip(section_ids, sig_columns)
    )
    title = escape(f"{project} -- interactive SBS-4608 plots")

    main_plot_type_options = "\n".join(
        f'<option value="{mpt}"{" selected" if mpt == main_plot_type else ""}>'
        f"{_MAIN_PLOT_TYPE_LABELS[mpt]}</option>"
        for mpt in _VALID_MAIN_PLOT_TYPES
    )
    main_plot_type_html = (
        '<div class="main-plot-type">\n'
        '<span class="label">Main Plot Type:</span>\n'
        '<select id="main-plot-type-select" onchange="applyMainPlotType()">\n'
        f"{main_plot_type_options}\n"
        "</select>\n"
        "</div>"
    )

    detailed_checkbox = (
        '<label><input type="checkbox" id="opt-detailed" checked onchange="applyPlotOptions()"> '
        "Detailed heatmaps (front/back)</label>"
        if include_detailed_heatmaps
        else ""
    )
    plot_options_html = (
        '<div class="plot-options">\n'
        '<span class="label">Plot options:</span>\n'
        '<label><input type="checkbox" id="opt-bias" checked onchange="applyPlotOptions()"> '
        "Transcription bias</label>\n"
        '<label><input type="checkbox" id="opt-total" checked onchange="applyPlotOptions()"> '
        "5-base context (total heatmap)</label>\n"
        f"{detailed_checkbox}\n"
        "</div>"
    )

    page_html = _SBS4608_PAGE_TEMPLATE.format(
        plotly_js=get_plotlyjs(),
        title=title,
        toc_links=toc_links,
        legend_swatches=legend_swatches,
        main_plot_type_html=main_plot_type_html,
        plot_options_html=plot_options_html,
        sections="\n".join(sections),
        plot_calls="\n".join(plot_calls),
        plot_registry=_json_for_script(plot_registry),
    )
    out_file.write_text(page_html, encoding="utf-8")

    return str(out_file.resolve())

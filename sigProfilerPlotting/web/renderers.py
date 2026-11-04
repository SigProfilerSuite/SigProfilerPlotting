"""Plotly JSON renderers used by the local dashboard."""

from __future__ import annotations

import re
from collections import OrderedDict

from .matrix import MatrixData, MatrixError

SBS_COLORS = {
    "C>A": "#03BCEE", "C>G": "#010101", "C>T": "#E32926",
    "T>A": "#CAC9C9", "T>C": "#A1CE63", "T>G": "#EBC6C4",
}
PALETTE = ["#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B", "#EECA3B", "#B279A2", "#FF9DA6"]

SBS96_SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
# The classic SigProfilerPlotting SBS96 display order groups all 16 contexts
# of one substitution type together (5' flank, then 3' flank) -- this is
# `seq96` in sigProfilerPlotting.py's plotSBS "96" branch, reproduced here as
# a generator rather than copying that 96-line literal. It differs from
# reference_formats/SBS96.txt's own row order (grouped by 5' flank first),
# which is why matrix.labels can't be used directly for this layout.
SBS96_DISPLAY_ORDER = [
    f"{flank5}[{sub}]{flank3}"
    for sub in SBS96_SUBSTITUTIONS
    for flank5 in "ACGT"
    for flank3 in "ACGT"
]


def _unit(values: list[float], unit: str) -> tuple[list[float], str, str]:
    if unit not in ("counts", "percent", "proportion"):
        raise MatrixError("unit must be 'counts', 'percent', or 'proportion'")
    if unit == "counts":
        return values, "Mutations", ".4~g"
    total = sum(values)
    if unit == "percent":
        return ([100.0 * value / total if total else 0.0 for value in values],
                "Percent of sample", ".3f")
    # "proportion": the same fraction-of-total as percent, kept on a 0-1
    # scale (no *100) rather than converted to a 0-100 percentage.
    return ([value / total if total else 0.0 for value in values],
            "Proportion of sample", ".4f")


def _layout(sample: str, context: str, ytitle: str, height: int = 520) -> dict:
    return {
        "title": {"text": f"{sample} — {context}", "x": 0.01, "xanchor": "left"},
        "height": height,
        "margin": {"l": 68, "r": 25, "t": 58, "b": 135},
        "paper_bgcolor": "#ffffff", "plot_bgcolor": "#ffffff",
        "hovermode": "closest", "barmode": "stack",
        "xaxis": {"tickangle": -90, "automargin": True, "showgrid": False},
        "yaxis": {"title": ytitle, "rangemode": "tozero", "gridcolor": "#e7ecf3"},
        # Right-aligned at the top-right corner, clear of the left-aligned
        # title above -- a horizontal legend at x=0 sat directly on top of it.
        "legend": {"orientation": "h", "y": 1.12, "x": 1, "xanchor": "right"},
    }


def _grouped(labels, values, grouping, colors, sample, context, ytitle, fmt):
    groups: OrderedDict[str, tuple[list[str], list[float]]] = OrderedDict()
    for label, value in zip(labels, values):
        group = grouping(label)
        xs, ys = groups.setdefault(group, ([], []))
        xs.append(label)
        ys.append(value)
    traces = []
    for index, (group, (xs, ys)) in enumerate(groups.items()):
        color = colors.get(group, PALETTE[index % len(PALETTE)])
        traces.append({
            "type": "bar", "name": group, "x": xs, "y": ys,
            "marker": {"color": color},
            "hovertemplate": f"%{{x}}<br>{sample}: %{{y:{fmt}}}<extra>{group}</extra>",
        })
    return {"traces": traces, "layout": _layout(sample, context, ytitle),
            "config": {"responsive": True, "displaylogo": False, "toImageButtonOptions": {"format": "svg"}}}


def _sbs96_group(label: str) -> str:
    match = re.search(r"\[([CT]>[ACGT])\]", label)
    return match.group(1) if match else "Other"


def _sbs96_tick_text(label: str, color: str) -> str:
    # Trinucleotide only (5' flank, mutated base, 3' flank), flanks in gray
    # and the mutated base bold in its group color -- e.g. "A[C>A]A" -> ACA
    # with the middle C colored/bold. Plotly tick labels support this
    # limited inline HTML.
    flank5, flank3 = label[0], label[-1]
    ref_base = label[2]
    return (
        f"<span style='color:#808080'>{flank5}</span>"
        f"<span style='color:{color}'><b>{ref_base}</b></span>"
        f"<span style='color:#808080'>{flank3}</span>"
    )


def _sbs96(matrix: MatrixData, sample: str, unit: str) -> dict:
    """The classic SigProfilerPlotting SBS96 look: six colored mutation-type
    header bands in place of a legend, trinucleotide-only x labels colored
    per group, and a bold sample title -- matching the static plotSBS "96"
    layout (sigProfilerPlotting.py's make_pickle_file/plotSBS) rather than
    the generic grouped-bar-with-legend style used for other contexts.
    """
    raw_by_label = dict(zip(matrix.labels, matrix.sample_values(sample)))
    raw = [raw_by_label[label] for label in SBS96_DISPLAY_ORDER]
    total = sum(raw)
    if unit == "counts":
        values = raw
        ytitle = "Number of Single Base Substitutions"
    elif unit == "percent":
        values = [100.0 * v / total if total else 0.0 for v in raw]
        ytitle = "Percentage of Single Base Substitutions"
    else:  # proportion
        values = [v / total if total else 0.0 for v in raw]
        ytitle = "Proportion of Single Base Substitutions"

    groups = [_sbs96_group(label) for label in SBS96_DISPLAY_ORDER]
    bar_colors = [SBS_COLORS[group] for group in groups]
    tick_text = [_sbs96_tick_text(label, color) for label, color in zip(SBS96_DISPLAY_ORDER, bar_colors)]

    top = max(values) * 1.1 if values and max(values) > 0 else 1.0
    tick_values = [top * i / 4 for i in range(5)]
    if unit == "counts":
        tick_labels = ["0"] + [f"{v:,.0f}" for v in tick_values[1:]]
    elif unit == "percent":
        tick_labels = ["0"] + [f"{v:.1f}%" for v in tick_values[1:]]
    else:
        tick_labels = ["0"] + [f"{v:.3f}" for v in tick_values[1:]]

    n = len(SBS96_DISPLAY_ORDER)
    per_group = n // len(SBS96_SUBSTITUTIONS)  # 16
    shapes = []
    annotations = []
    for index, substitution in enumerate(SBS96_SUBSTITUTIONS):
        x0, x1 = index * per_group - 0.5, (index + 1) * per_group - 0.5
        color = SBS_COLORS[substitution]
        shapes.append({
            "type": "rect", "xref": "x", "yref": "paper",
            "x0": x0, "x1": x1, "y0": 1.03, "y1": 1.1,
            "fillcolor": color, "line": {"width": 0},
        })
        annotations.append({
            "xref": "x", "yref": "paper", "x": (x0 + x1) / 2, "y": 1.17,
            "text": f"<b>{substitution}</b>", "showarrow": False,
            "font": {"family": "Arial, sans-serif", "size": 15, "color": "#000000"},
        })
    annotations.append({
        # Inside the plotting area (below the header bands at paper y>1),
        # near its top-left -- not at y=1 exactly, which the bands' own
        # shapes (y0=1.03) would otherwise clip into.
        "xref": "paper", "yref": "paper", "x": 0.01, "y": 0.96, "yanchor": "top",
        "xanchor": "left", "text": f"<b>{sample}</b>", "showarrow": False,
        "font": {"family": "Arial, sans-serif", "size": 18, "color": "#000000"},
    })

    layout = {
        "height": 560,
        "margin": {"l": 70, "r": 25, "t": 100, "b": 70},
        "paper_bgcolor": "#ffffff", "plot_bgcolor": "#ffffff",
        "showlegend": False,
        "shapes": shapes,
        "annotations": annotations,
        "xaxis": {
            "tickmode": "array",
            "tickvals": list(range(n)),
            "ticktext": tick_text,
            "tickfont": {"family": "Courier New, monospace", "size": 11},
            "showgrid": False,
            "range": [-0.5, n - 0.5],
        },
        "yaxis": {
            "title": {"text": ytitle, "font": {"family": "Times New Roman, serif", "size": 15}},
            "tickmode": "array",
            "tickvals": tick_values,
            "ticktext": tick_labels,
            "range": [0, top],
            "gridcolor": "#ededed",
            "zeroline": False,
        },
    }
    traces = [{
        "type": "bar", "x": list(range(n)), "y": values,
        "marker": {"color": bar_colors},
        "customdata": SBS96_DISPLAY_ORDER,
        "hovertemplate": f"%{{customdata}}<br>{sample}: %{{y}}<extra></extra>",
        "showlegend": False,
    }]
    return {
        "traces": traces, "layout": layout,
        "config": {"responsive": True, "displaylogo": False, "toImageButtonOptions": {"format": "svg"}},
    }


def _sbs288(matrix, sample, values, ytitle, fmt):
    colors = {"T": "#2878B5", "U": "#D9534F", "N": "#8C8C8C"}
    figure = _grouped(matrix.labels, values, lambda label: label.split(":", 1)[0], colors,
                      sample, matrix.context, ytitle, fmt)
    figure["layout"]["barmode"] = "group"
    return figure


def _dbs78(matrix, sample, values, ytitle, fmt):
    return _grouped(matrix.labels, values, lambda label: label.split(">", 1)[0], {},
                    sample, matrix.context, ytitle, fmt)


def _id83(matrix, sample, values, ytitle, fmt):
    def group(label):
        parts = label.split(":")
        return ":".join(parts[1:3]) if len(parts) >= 3 else parts[0]
    colors = {"Del:C": "#F4A582", "Del:T": "#D6604D", "Ins:C": "#92C5DE",
              "Ins:T": "#4393C3", "Del:R": "#B2182B", "Ins:R": "#2166AC", "Del:M": "#7B3294"}
    return _grouped(matrix.labels, values, group, colors, sample, matrix.context, ytitle, fmt)


def _cnv48(matrix, sample, values, ytitle, fmt):
    colors = {"homdel": "#4D4D4D", "LOH": "#377EB8", "het": "#4DAF4A"}
    return _grouped(matrix.labels, values, lambda label: label.split(":")[1], colors,
                    sample, matrix.context, ytitle, fmt)


def _sv32(matrix, sample, values, ytitle, fmt):
    colors = {"clustered": "#E45756", "non-clustered": "#4C78A8"}
    figure = _grouped(matrix.labels, values,
                      lambda label: "non-clustered" if label.startswith("non-clustered_") else "clustered",
                      colors, sample, matrix.context, ytitle, fmt)
    figure["layout"]["barmode"] = "group"
    return figure


def _generic(matrix, sample, values, ytitle, fmt):
    def family(label):
        # A leading segment that's purely numeric (e.g. "1:Del:C:0" for
        # ID28) is a positional index, not a semantic family -- skip it and
        # use the next colon-delimited token instead.
        segments = label.split(":") if ":" in label else []
        semantic = next((segment for segment in segments if not segment.isdigit()), None)
        if semantic is not None:
            return semantic
        if ">" in label:
            return label.split(">", 1)[0]
        match = re.search(r"\[([^]]+)\]", label)
        if match:
            return match.group(1)
        # No usable family token anywhere in this label set: fall back to one
        # ungrouped bucket rather than grouping by a meaningless digit/index.
        return matrix.context

    figure = _grouped(matrix.labels, values, family, {}, sample, matrix.context, ytitle, fmt)
    figure["layout"]["height"] = 570
    return figure


def _sbs4608(matrix: MatrixData, sample: str, unit: str) -> dict:
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig

    # `unit` is intentionally unused here: SBS4608Data/SBS4608PlotBuilder
    # already renormalize every contribution to a percentage of its own
    # total internally (matching plotInteractive, which has no raw-counts
    # mode at all), so pre-dividing by total here would be a no-op -- the
    # dashboard UI hides the counts/percentage/proportion control for this
    # context rather than offering a choice that can't change the output.
    raw = matrix.raw4608(sample)
    variants = {}
    for main_type in ("default", "strand_bias", "genic_bias"):
        config = SBS4608PlotConfig(include_detailed_heatmaps=True, main_plot_type=main_type)
        built = SBS4608PlotBuilder(SBS4608Data(raw, sample, main_type), config).build()
        total_top = config.total_heatmap_domain[1]
        built["sections"]["total_heatmap"]["paper_band"] = [0.19, total_top]
        built["sections"]["detailed_heatmaps"]["paper_band"] = [total_top, config.main_plot_domain[0]]
        built["sections"]["transcription_bias"]["paper_band"] = [0, 0.19]
        built["config"].update({"responsive": True, "displayModeBar": True, "displaylogo": False})
        built["main_plot_type"] = main_type  # read by the ported compactPlot() in app.js
        variants[main_type] = built
    return {"kind": "sbs4608", "variants": variants}


def render_sample(matrix: MatrixData, sample: str, unit: str = "counts") -> dict:
    """Render one sample/signature as JSON accepted by Plotly.newPlot."""
    if unit not in ("counts", "percent", "proportion"):
        raise MatrixError("unit must be 'counts', 'percent', or 'proportion'")
    if matrix.context == "SBS4608":
        return _sbs4608(matrix, sample, unit)
    if matrix.context == "SBS96":
        # Bypasses the generic _unit()/_grouped() path entirely: SBS96 gets
        # its own classic SigProfilerPlotting layout (header bands, colored
        # trinucleotide ticks) rather than the shared legend-based style.
        return {"kind": "standard", "figure": _sbs96(matrix, sample, unit)}
    raw = matrix.sample_values(sample)
    values, ytitle, fmt = _unit(raw, unit)
    renderer = {
        "SBS288": _sbs288, "DBS78": _dbs78,
        "ID83": _id83, "CNV48": _cnv48, "SV32": _sv32,
    }.get(matrix.context, _generic)
    return {"kind": "standard", "figure": renderer(matrix, sample, values, ytitle, fmt)}

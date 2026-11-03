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


def _sbs96(matrix, sample, values, ytitle, fmt):
    def group(label):
        match = re.search(r"\[([CT]>[ACGT])\]", label)
        return match.group(1) if match else "Other"
    return _grouped(matrix.labels, values, group, SBS_COLORS, sample, matrix.context, ytitle, fmt)


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
        variants[main_type] = built
    return {"kind": "sbs4608", "variants": variants}


def render_sample(matrix: MatrixData, sample: str, unit: str = "counts") -> dict:
    """Render one sample/signature as JSON accepted by Plotly.newPlot."""
    raw = matrix.sample_values(sample)
    values, ytitle, fmt = _unit(raw, unit)
    if matrix.context == "SBS4608":
        return _sbs4608(matrix, sample, unit)
    renderer = {
        "SBS96": _sbs96, "SBS288": _sbs288, "DBS78": _dbs78,
        "ID83": _id83, "CNV48": _cnv48, "SV32": _sv32,
    }.get(matrix.context, _generic)
    return {"kind": "standard", "figure": renderer(matrix, sample, values, ytitle, fmt)}

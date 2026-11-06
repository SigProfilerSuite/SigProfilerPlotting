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


def _sbs96_figure(
    raw_by_label: dict, title: str, unit: str, height: int = 560, title_font_size: int = 18
) -> dict:
    """The classic SigProfilerPlotting SBS96 look: six colored mutation-type
    header bands in place of a legend, trinucleotide-only x labels colored
    per group, and a bold title -- matching the static plotSBS "96" layout
    (sigProfilerPlotting.py's make_pickle_file/plotSBS) rather than the
    generic grouped-bar-with-legend style used for other contexts.

    Takes a plain {96-context label: value} mapping rather than a
    MatrixData/sample pair so it can draw any SBS96-shaped vector -- a
    dashboard sample/signature, or (smaller, for the COSMIC comparison
    panels) a reference signature -- through the same code path.
    """
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
    margin_top, margin_bottom = 100, 70
    # yref="paper" fractions above 1.0 (or below 1.0, for the title) scale
    # with the AXES height, not with the figure's total height -- fixed
    # fractions tuned for the default 560px figure (390px of axes) put the
    # header bands and title comfortably apart there, but at a much shorter
    # figure (the COSMIC comparison panels, ~300px) the same fractions
    # resolve to far fewer pixels while the text stays the same size,
    # visually overlapping. Anchor these to fixed pixel offsets instead, so
    # the spacing is identical regardless of `height`.
    axes_height = max(height - margin_top - margin_bottom, 1)

    def above_axes(pixels):
        return 1 + pixels / axes_height

    def inside_axes(pixels):
        return 1 - pixels / axes_height

    shapes = []
    annotations = []
    for index, substitution in enumerate(SBS96_SUBSTITUTIONS):
        x0, x1 = index * per_group - 0.5, (index + 1) * per_group - 0.5
        color = SBS_COLORS[substitution]
        shapes.append({
            "type": "rect", "xref": "x", "yref": "paper",
            "x0": x0, "x1": x1, "y0": above_axes(12), "y1": above_axes(39),
            "fillcolor": color, "line": {"width": 0},
        })
        annotations.append({
            "xref": "x", "yref": "paper", "x": (x0 + x1) / 2, "y": above_axes(66),
            "text": f"<b>{substitution}</b>", "showarrow": False,
            "font": {"family": "Arial, sans-serif", "size": 15, "color": "#000000"},
        })
    annotations.append({
        # Inside the plotting area (below the header bands at paper y>1),
        # near its top-left -- not at y=1 exactly, which the bands' own
        # shapes would otherwise clip into.
        "xref": "paper", "yref": "paper", "x": 0.01, "y": inside_axes(16), "yanchor": "top",
        "xanchor": "left", "text": f"<b>{title}</b>", "showarrow": False,
        "font": {"family": "Arial, sans-serif", "size": title_font_size, "color": "#000000"},
    })

    layout = {
        "height": height,
        "margin": {"l": 70, "r": 25, "t": margin_top, "b": margin_bottom},
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
        "hovertemplate": f"%{{customdata}}<br>{title}: %{{y}}<extra></extra>",
        "showlegend": False,
    }]
    return {
        "traces": traces, "layout": layout,
        "config": {"responsive": True, "displaylogo": False, "toImageButtonOptions": {"format": "svg"}},
    }


def _sbs96(matrix: MatrixData, sample: str, unit: str) -> dict:
    raw_by_label = dict(zip(matrix.labels, matrix.sample_values(sample)))
    return _sbs96_figure(raw_by_label, sample, unit)


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


# ID89/ID476 (Koh et al. 2025) group the same way the static plotID()
# contexts do in sigProfilerPlotting.py -- duplicated here, in pure-Python
# string matching with no matplotlib dependency, to keep this module's
# only import light (see module docstring); the two implementations are
# kept in sync by tests/test_web_renderers.py comparing group names/sizes
# against reference_formats/ID89.txt and ID476.txt directly.
_ID89_T_BUCKETS = (
    ("R(1,4)", "(1-4)"), ("R(5,7)", "(5-7)"), ("R(8,)", "(8+)"), ("R(0,4)", "(0-4)"),
)


def _id89_group(label: str) -> str:
    if label == "Complex":
        return "Complex"
    if "Del(C)" in label:
        return "Del 1 C"
    if "Ins(C)" in label:
        return "Ins 1 C"
    for token, bucket in _ID89_T_BUCKETS:
        if token in label:
            if "Del(T)" in label:
                return f"Del 1 T {bucket}"
            if "Ins(T)" in label:
                return f"Ins 1 T {bucket}"
    if re.search(r":M\d", label) or re.search(r":M\(\d", label):
        return "Del Mh"
    return "Del >=2" if label.startswith("Del") else "Ins >=2"


def _id476_group(label: str) -> str:
    if label.startswith("Complex"):
        return "Complex"
    match = re.search(r"\[(Del|Ins)\((C|T)\)", label)
    if match:
        return f"{match.group(1)} 1 {match.group(2)}"
    if re.search(r":M\d", label) or re.search(r":M\(\d", label):
        return "Mh"
    return "Del >=2" if label.startswith("Del") else "Ins >=2"


def _id89(matrix, sample, values, ytitle, fmt):
    colors = {
        "Del 1 C": "#FDBE6F", "Del 1 T (1-4)": "#FD8D3C", "Del 1 T (5-7)": "#E6550D",
        "Del 1 T (8+)": "#A63603", "Ins 1 C": "#C7E9C0", "Ins 1 T (0-4)": "#A1D99B",
        "Ins 1 T (5-7)": "#41AB5D", "Ins 1 T (8+)": "#006D2C", "Del >=2": "#DE2D26",
        "Ins >=2": "#3182BD", "Del Mh": "#756BB1", "Complex": "#969696",
    }
    return _grouped(matrix.labels, values, _id89_group, colors, sample, matrix.context, ytitle, fmt)


def _id476(matrix, sample, values, ytitle, fmt):
    colors = {
        "Del 1 C": "#FDAE61", "Del 1 T": "#E6550D", "Ins 1 C": "#C7E9C0",
        "Ins 1 T": "#41AB5D", "Del >=2": "#DE2D26", "Ins >=2": "#3182BD",
        "Mh": "#756BB1", "Complex": "#969696",
    }
    figure = _grouped(matrix.labels, values, _id476_group, colors, sample, matrix.context, ytitle, fmt)
    figure["layout"]["height"] = 570
    return figure


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


def _scale_sbs4608_unit(built: dict, total: float, unit: str) -> None:
    """Rescale an already-built SBS4608 figure out of its internal percent-
    of-`total` representation.

    SBS4608Data/SBS4608PlotBuilder (shared with plotInteractive(), which
    keeps its own percent-only COSMIC-style look unchanged) always
    renormalize every contribution to a percentage of `total` first. Every
    number this function touches is ultimately a sum of those per-row
    percentages, i.e. `contribution * (100 / total)` for one fixed `total`
    -- so undoing that scale factor is a uniform linear transform, safe to
    apply after the fact rather than re-deriving each view from raw counts.

    Traces left untouched on purpose: the "Strand Bias (%)"/"Genic (%)"
    sub-panel (yaxis5/yaxis6, scatter traces) and its text labels are
    ratios between two sibling bars (e.g. Transcribed vs Untranscribed), not
    shares of `total` -- they were never percent-of-total to begin with, so
    "counts"/"proportion" must not change them.
    """
    if unit == "percent" or total <= 0:
        return
    factor = (total / 100.0) if unit == "counts" else (1.0 / 100.0)
    if unit == "counts":
        value_text = "%{z:,.0f}"
        tick_text = lambda v: f"{v:,.0f}"  # noqa: E731
    else:
        value_text = "%{z:.4f}"
        tick_text = lambda v: f"{v:.4f}"  # noqa: E731

    for trace in built["traces"]:
        if trace.get("type") == "heatmap":
            if trace.get("z"):
                trace["z"] = [[cell * factor for cell in row] for row in trace["z"]]
            if "zmin" in trace:
                trace["zmin"] = trace["zmin"] * factor
            if "zmax" in trace:
                trace["zmax"] = trace["zmax"] * factor
            colorbar = trace.get("colorbar")
            if colorbar and colorbar.get("tickmode") == "array":
                colorbar["tickvals"] = [v * factor for v in colorbar["tickvals"]]
                colorbar["ticktext"] = [tick_text(v) for v in colorbar["tickvals"]]
            trace["hovertemplate"] = f"x: %{{x}}<br>y: %{{y}}<br>Value: {value_text}<extra></extra>"
        elif trace.get("type") == "bar" and trace.get("yaxis") in (None, "y"):
            # The main plot's own bar trace(s): either the default 96-bar
            # chart, or the strand/genic-bias stacked bars -- both are
            # percent-of-total (not the ratio-only transcription-bias panel,
            # whose traces are type "scatter" on yaxis5/yaxis6).
            if trace.get("y"):
                trace["y"] = [v * factor for v in trace["y"]]

    yaxis = built["layout"].get("yaxis")
    if yaxis:
        title = (
            "Number of Single Base Substitutions"
            if unit == "counts"
            else "Proportion of Single Base Substitutions"
        )
        yaxis["title"] = {**yaxis.get("title", {}), "text": title}
        yaxis.pop("ticksuffix", None)
        yaxis["tickformat"] = "~s" if unit == "counts" else ".3f"
        if "range" in yaxis:
            yaxis["range"] = [v * factor for v in yaxis["range"]]


def _reclaim_bottom_band(built: dict, band_hi: float):
    """Permanently remap every paper-fraction position out of a [0, band_hi]
    band, the same way the client's compactPlot() remaps a toggled-off
    section -- for a section that is never built at all (no checkbox to
    toggle it back on), the remap only ever needs to happen once, here,
    rather than shipping an always-empty section for the client to reclaim.

    Returns the same `map_y` used internally, so callers can also remap any
    of SBS4608PlotConfig's own (pre-reclaim, now-stale) domain fractions --
    e.g. the paper_band boundaries of sections that remain real and
    toggleable, like SBS1536's total/detailed heatmaps.
    """
    remaining = 1.0 - band_hi

    def map_y(y):
        return (y - min(y, band_hi)) / remaining

    layout = built["layout"]
    for key, axis in layout.items():
        if re.fullmatch(r"yaxis\d*", key) and isinstance(axis, dict) and "domain" in axis:
            axis["domain"] = [map_y(v) for v in axis["domain"]]
    for annotation in layout.get("annotations", []):
        if annotation.get("yref") == "paper":
            annotation["y"] = map_y(annotation["y"])
    for shape in layout.get("shapes", []):
        if shape.get("yref") == "paper":
            shape["y0"] = map_y(shape["y0"])
            shape["y1"] = map_y(shape["y1"])
    for trace in built["traces"]:
        colorbar = trace.get("colorbar")
        if colorbar:
            if "y" in colorbar:
                colorbar["y"] = map_y(colorbar["y"])
            if "len" in colorbar and colorbar.get("lenmode") != "pixels":
                colorbar["len"] = colorbar["len"] / remaining
    legend = layout.get("legend")
    if legend and "y" in legend:
        legend["y"] = map_y(legend["y"])
    margin = layout.get("margin") or {}
    vertical_margins = margin.get("t", 100) + margin.get("b", 80)
    layout["height"] = round(vertical_margins + (layout["height"] - vertical_margins) * remaining)
    return map_y


def _sbs4608(matrix: MatrixData, sample: str, unit: str) -> dict:
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig

    raw = matrix.raw4608(sample)
    total = sum(float(row["contribution"]) for row in raw)
    variants = {}
    for main_type in ("default", "strand_bias", "genic_bias"):
        config = SBS4608PlotConfig(include_detailed_heatmaps=True, main_plot_type=main_type)
        built = SBS4608PlotBuilder(SBS4608Data(raw, sample, main_type), config).build()
        _scale_sbs4608_unit(built, total, unit)
        total_top = config.total_heatmap_domain[1]
        built["sections"]["total_heatmap"]["paper_band"] = [0.19, total_top]
        built["sections"]["detailed_heatmaps"]["paper_band"] = [total_top, config.main_plot_domain[0]]
        built["sections"]["transcription_bias"]["paper_band"] = [0, 0.19]
        built["config"].update({"responsive": True, "displayModeBar": True, "displaylogo": False})
        built["main_plot_type"] = main_type  # read by the ported compactPlot() in app.js
        variants[main_type] = built
    return {"kind": "sbs4608", "variants": variants}


def _sbs1536(matrix: MatrixData, sample: str, unit: str) -> dict:
    # SBS1536's 1536 pentanucleotide contexts are exactly what SBS4608Data
    # collapses its strand-tagged rows down to once strand is summed away --
    # so tag every row with one fixed, made-up strand ("T:") and SBS4608Data's
    # own strand-sum becomes a no-op (nothing else to sum over), letting
    # SBS4608PlotBuilder's bar chart / total heatmap / front+back detailed
    # heatmaps run completely unmodified on real SBS1536 data.
    #
    # What's genuinely different -- there is no real strand/transcription
    # data for this context -- is handled explicitly below: the
    # transcription-bias traces are never built (only the fabricated T tag
    # exists, so a real Transcribed/Untranscribed/Genic/Intergenic split
    # would be meaningless), and the total heatmap's hover text has its
    # (equally fabricated) T/U/N breakdown stripped.
    from sigProfilerPlotting.sbs4608.annotation_builder import AnnotationBuilder
    from sigProfilerPlotting.sbs4608.layout_builder import LayoutBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_data import SBS4608Data
    from sigProfilerPlotting.sbs4608.sbs4608_plot_builder import SBS4608PlotBuilder
    from sigProfilerPlotting.sbs4608.sbs4608_plot_config import SBS4608PlotConfig
    from sigProfilerPlotting.sbs4608.tick_label_formatter import TickLabelFormatter
    from sigProfilerPlotting.sbs4608.trace_formatter import TraceFormatter

    raw = matrix.raw4608(sample)
    total = sum(float(row["contribution"]) for row in raw)
    fake_raw4608 = [
        {"mutationType": f"T:{row['mutationType']}", "contribution": row["contribution"]}
        for row in raw
    ]
    data = SBS4608Data(fake_raw4608, sample, "default")
    config = SBS4608PlotConfig(include_detailed_heatmaps=True, main_plot_type="default")
    builder = SBS4608PlotBuilder(data, config)

    traces = builder._build_main_plot_traces()
    heatmap_start = len(traces)
    total_traces = builder.heatmap_builder.build_total_heatmap(
        data.group_by_outer_with_tun_proportions(), data.get_total_mutations()
    )
    for trace in total_traces:
        trace["hovertemplate"] = "%{customdata.context}: %{z}<extra></extra>"
    traces.extend(total_traces)
    total_heatmap_range = [heatmap_start, len(traces)]

    detailed_start = len(traces)
    traces.extend(builder._build_detailed_heatmap_traces())
    detailed_heatmaps_range = [detailed_start, len(traces)]

    TraceFormatter.apply_percentage_formatting(traces)
    TraceFormatter.configure_heatmap_colorbar(traces, config.should_include_detailed_heatmaps())

    bar_data = data.get_96_bar_chart_data()
    annotations = builder._build_annotations(bar_data)
    shapes = builder._build_shapes(bar_data)
    layout = builder.layout_builder.build_base_layout(
        bar_data["flatSorted"], bar_data["maxValue"], annotations, shapes
    )
    builder._configure_axes(layout, None)
    layout["barmode"] = "group"
    layout["plot_bgcolor"] = "white"
    layout["paper_bgcolor"] = "white"
    layout["height"] = config.plot_height
    layout["width"] = config.plot_width
    LayoutBuilder.apply_percentage_layout_modifications(
        layout, include_strand_bias=False,
        include_detailed_heatmaps=config.should_include_detailed_heatmaps(),
    )
    TickLabelFormatter.apply_to_layout(layout, data.get_96_data(), config.get_colors())
    layout["annotations"] = AnnotationBuilder.remove_detailed_trinucleotide_labels(layout["annotations"])
    AnnotationBuilder.add_missing_ca_label(layout["annotations"])

    built = {
        "traces": traces,
        "layout": layout,
        "config": {"responsive": True, "displayModeBar": True, "displaylogo": False},
        "sections": {
            "total_heatmap": {"trace_range": total_heatmap_range},
            "detailed_heatmaps": {"trace_range": detailed_heatmaps_range},
            "transcription_bias": {"trace_range": [len(traces), len(traces)], "annotation_indices": []},
        },
        "main_plot_type": "default",
    }
    # There is no transcription-bias panel for SBS1536 at all -- reclaim the
    # paper-fraction band SBS4608PlotConfig always reserves for it, instead
    # of shipping an always-empty, never-toggleable section that just leaves
    # blank space at the bottom of every plot.
    map_y = _reclaim_bottom_band(built, 0.19)
    _scale_sbs4608_unit(built, total, unit)
    # SBS4608PlotConfig's own domain fractions are the pre-reclaim ones;
    # remap them through the same transform so the remaining two sections'
    # paper_band stays consistent with where compactPlot() will find them.
    total_top = config.total_heatmap_domain[1]
    built["sections"]["total_heatmap"]["paper_band"] = [map_y(0.19), map_y(total_top)]
    built["sections"]["detailed_heatmaps"]["paper_band"] = [
        map_y(total_top), map_y(config.main_plot_domain[0])
    ]
    return {"kind": "sbs4608", "variants": {"default": built}}


# Contexts the COSMIC comparison is offered for -- any SBS context with a
# well-defined trinucleotide collapse (cosmic.collapse_to_sbs96 strips a
# strand tag and/or wider pentanucleotide flanks down to the 96 contexts
# COSMIC's own catalogue is expressed in).
COSMIC_COMPATIBLE_CONTEXTS = ("SBS96", "SBS288", "SBS1536", "SBS4608")


def _cosmic_compare_unit(unit: str) -> str:
    # A reference signature has no "raw count" of its own -- it is already a
    # normalized profile -- so "counts" falls back to "proportion" rather
    # than being meaningless. "percent" and "proportion" both make sense and
    # pass straight through, matching whatever the reader has the main
    # Values control set to.
    return unit if unit in ("percent", "proportion") else "proportion"


def _cosmic_compare(matrix: MatrixData, sample: str, unit: str, build: str) -> dict:
    """The 5 COSMIC signatures (from `build`) most cosine-similar to this
    sample/signature, each drawn through the same classic SBS96 renderer as
    the query (at full size, so they're directly, visually comparable) and
    ranked highest-similarity first.

    Works for any context in COSMIC_COMPATIBLE_CONTEXTS: raw_by_label is
    collapsed to the 96 trinucleotide contexts first (a no-op for an
    already-96-context SBS96 matrix).
    """
    from . import cosmic

    raw_by_label = dict(zip(matrix.labels, matrix.sample_values(sample)))
    collapsed = cosmic.collapse_to_sbs96(raw_by_label)
    matches = cosmic.top_matches(collapsed, build)
    match_unit = _cosmic_compare_unit(unit)
    match_figures = [
        {
            "name": match["name"],
            "score": match["score"],
            "figure": _sbs96_figure(
                match["values"], f"{match['name']} (cos {match['score']:.3f})", match_unit
            ),
        }
        for match in matches
    ]
    return {"build": build, "matches": match_figures}


def render_sample(
    matrix: MatrixData, sample: str, unit: str = "counts", compare_cosmic: str | None = None
) -> dict:
    """Render one sample/signature as JSON accepted by Plotly.newPlot."""
    if unit not in ("counts", "percent", "proportion"):
        raise MatrixError("unit must be 'counts', 'percent', or 'proportion'")
    if matrix.context == "SBS4608":
        response = _sbs4608(matrix, sample, unit)
    elif matrix.context == "SBS1536":
        response = _sbs1536(matrix, sample, unit)
    elif matrix.context == "SBS96":
        # Bypasses the generic _unit()/_grouped() path entirely: SBS96 gets
        # its own classic SigProfilerPlotting layout (header bands, colored
        # trinucleotide ticks) rather than the shared legend-based style.
        response = {"kind": "standard", "figure": _sbs96(matrix, sample, unit)}
    else:
        raw = matrix.sample_values(sample)
        values, ytitle, fmt = _unit(raw, unit)
        renderer = {
            "SBS288": _sbs288, "DBS78": _dbs78,
            "ID83": _id83, "ID89": _id89, "ID476": _id476,
            "CNV48": _cnv48, "SV32": _sv32,
        }.get(matrix.context, _generic)
        response = {"kind": "standard", "figure": renderer(matrix, sample, values, ytitle, fmt)}

    if compare_cosmic and matrix.context in COSMIC_COMPATIBLE_CONTEXTS:
        # An addition alongside whatever the context's own "kind" already
        # is, not a separate response shape -- so every context keeps its
        # normal rendering (SBS4608/SBS1536's multi-panel view included)
        # with the comparison attached underneath it.
        response["cosmic_compare"] = _cosmic_compare(matrix, sample, unit, compare_cosmic)
    return response

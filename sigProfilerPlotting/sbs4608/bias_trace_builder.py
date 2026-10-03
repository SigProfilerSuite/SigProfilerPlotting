"""Python port of BiasTraceBuilder.js."""

import re

MUTATION_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def _get_mutation(ctx):
    return re.search(r"\[(.*)\]", ctx).group(1)


def _sort_by_mutation_group(item):
    ctx = item["context96"]
    mut = _get_mutation(ctx)
    order = MUTATION_ORDER.index(mut) if mut in MUTATION_ORDER else len(MUTATION_ORDER)
    return (order, ctx)


class BiasTraceBuilder:
    @staticmethod
    def build_strand_bias_traces(data_192_strand, colors):
        transcribed_data = []
        untranscribed_data = []

        for item in data_192_strand:
            strand, context96 = item["mutationType"].split(":", 1)
            if strand == "T":
                transcribed_data.append(
                    {"context96": context96, "value": item["mutations"]}
                )
            elif strand == "U":
                untranscribed_data.append(
                    {"context96": context96, "value": item["mutations"]}
                )

        transcribed_data.sort(key=_sort_by_mutation_group)
        untranscribed_data.sort(key=_sort_by_mutation_group)

        transcribed_x = list(range(len(transcribed_data)))
        transcribed_y = [d["value"] for d in transcribed_data]
        untranscribed_x = list(range(len(untranscribed_data)))
        untranscribed_y = [d["value"] for d in untranscribed_data]

        return [
            {
                "name": "Transcribed",
                "type": "bar",
                "marker": {"color": colors["transcribed"]},
                "x": transcribed_x,
                "y": transcribed_y,
                "hoverinfo": "x+y",
                "showlegend": True,
                "legendgroup": "strand",
            },
            {
                "name": "Untranscribed",
                "type": "bar",
                "marker": {"color": colors["untranscribed"]},
                "x": untranscribed_x,
                "y": untranscribed_y,
                "hoverinfo": "x+y",
                "showlegend": True,
                "legendgroup": "strand",
            },
        ]

    @staticmethod
    def build_genic_bias_traces(data_192_genic, colors):
        genic_data = []
        intergenic_data = []

        for item in data_192_genic:
            genic_status, context96 = item["mutationType"].split(":", 1)
            if genic_status == "Genic":
                genic_data.append({"context96": context96, "value": item["mutations"]})
            elif genic_status == "Intergenic":
                intergenic_data.append(
                    {"context96": context96, "value": item["mutations"]}
                )

        genic_data.sort(key=_sort_by_mutation_group)
        intergenic_data.sort(key=_sort_by_mutation_group)

        genic_x = list(range(len(genic_data)))
        genic_y = [d["value"] for d in genic_data]
        intergenic_x = list(range(len(intergenic_data)))
        intergenic_y = [d["value"] for d in intergenic_data]

        return [
            {
                "name": "Genic",
                "type": "bar",
                "marker": {"color": colors["genic"]},
                "x": genic_x,
                "y": genic_y,
                "hoverinfo": "x+y",
                "showlegend": True,
                "legendgroup": "genic",
            },
            {
                "name": "Intergenic",
                "type": "bar",
                "marker": {"color": colors["intergenic"]},
                "x": intergenic_x,
                "y": intergenic_y,
                "hoverinfo": "x+y",
                "showlegend": True,
                "legendgroup": "genic",
            },
        ]

    @staticmethod
    def build_stacked_transcription_bias_traces(data, mutation_colors):
        mutations = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        bar_border_color = "#CDCDCD"
        label_font_size = 10

        def get_raw_data(mutation, abbreviation):
            for d in data:
                if (
                    d["mutationType"][0] == abbreviation
                    and mutation in d["mutationType"]
                ):
                    return d["contribution"]
            return 0

        traces = []

        for mut_index, mutation in enumerate(mutations):
            T = get_raw_data(mutation, "T")
            U = get_raw_data(mutation, "U")
            N = get_raw_data(mutation, "N")

            total = T + U + N
            if total <= 0:
                continue

            genic_prop = (T + U) / total
            intergenic_prop = N / total
            transcribed_prop = T / (T + U) if (T + U) else 0
            untranscribed_prop = U / (T + U) if (T + U) else 0

            genic_perc = genic_prop * 100
            intergenic_perc = intergenic_prop * 100
            transcribed_perc = transcribed_prop * 100
            untranscribed_perc = untranscribed_prop * 100

            bar_y = 0.5
            total_groups = 6
            padding_per_group = 0.015
            total_padding = total_groups * padding_per_group * 2
            available_width = 1 - total_padding
            bar_width = available_width / total_groups
            x_start = (
                mut_index * (bar_width + 2 * padding_per_group) + padding_per_group
            )
            bar_height = 0.25

            base_color = mutation_colors[mutation]["base"]
            pale_color = mutation_colors[mutation]["pale"]
            transcribed_color = base_color
            untranscribed_color = pale_color
            genic_color = base_color
            intergenic_color = pale_color

            base_text_color = "white"
            pale_text_color = "black"

            # Transcriptional strand bias plot (y5)
            transcribed_bar_width = transcribed_prop * bar_width
            traces.append(
                {
                    "x": [
                        x_start,
                        x_start + transcribed_bar_width,
                        x_start + transcribed_bar_width,
                        x_start,
                        x_start,
                    ],
                    "y": [
                        bar_y - bar_height,
                        bar_y - bar_height,
                        bar_y + bar_height,
                        bar_y + bar_height,
                        bar_y - bar_height,
                    ],
                    "fill": "toself",
                    "type": "scatter",
                    "mode": "lines",
                    "line": {"color": bar_border_color, "width": 1},
                    "fillcolor": transcribed_color,
                    "name": "",
                    "legendgroup": "transcription",
                    "legend": "legend",
                    "xaxis": "x2",
                    "yaxis": "y5",
                    "showlegend": False,
                    "text": f"{mutation} - Transcribed: {transcribed_perc:.1f}%",
                    "hoverinfo": "text",
                    "hovertemplate": f"{mutation} - Transcribed: {transcribed_perc:.1f}%<extra></extra>",
                }
            )

            untranscribed_bar_width = untranscribed_prop * bar_width
            untranscribed_start = x_start + transcribed_bar_width
            traces.append(
                {
                    "x": [
                        untranscribed_start,
                        untranscribed_start + untranscribed_bar_width,
                        untranscribed_start + untranscribed_bar_width,
                        untranscribed_start,
                        untranscribed_start,
                    ],
                    "y": [
                        bar_y - bar_height,
                        bar_y - bar_height,
                        bar_y + bar_height,
                        bar_y + bar_height,
                        bar_y - bar_height,
                    ],
                    "fill": "toself",
                    "type": "scatter",
                    "mode": "lines",
                    "line": {"color": bar_border_color, "width": 1},
                    "fillcolor": untranscribed_color,
                    "name": "",
                    "legendgroup": "transcription",
                    "legend": "legend",
                    "xaxis": "x2",
                    "yaxis": "y5",
                    "showlegend": False,
                    "text": f"{mutation} - Untranscribed: {untranscribed_perc:.1f}%",
                    "hoverinfo": "text",
                    "hovertemplate": f"{mutation} - Untranscribed: {untranscribed_perc:.1f}%<extra></extra>",
                }
            )

            if transcribed_bar_width > 0.02:
                traces.append(
                    {
                        "x": [x_start + transcribed_bar_width / 2],
                        "y": [bar_y],
                        "mode": "text",
                        "type": "scatter",
                        "text": [f"{transcribed_perc:.1f}%"],
                        "textfont": {"color": base_text_color, "size": label_font_size},
                        "showlegend": False,
                        "xaxis": "x2",
                        "yaxis": "y5",
                        "hoverinfo": "skip",
                    }
                )

            if untranscribed_bar_width > 0.02:
                traces.append(
                    {
                        "x": [untranscribed_start + untranscribed_bar_width / 2],
                        "y": [bar_y],
                        "mode": "text",
                        "type": "scatter",
                        "text": [f"{untranscribed_perc:.1f}%"],
                        "textfont": {"color": pale_text_color, "size": label_font_size},
                        "showlegend": False,
                        "xaxis": "x2",
                        "yaxis": "y5",
                        "hoverinfo": "skip",
                    }
                )

            # Genic/intergenic bias plot (y6)
            genic_bar_width = genic_prop * bar_width
            traces.append(
                {
                    "x": [
                        x_start,
                        x_start + genic_bar_width,
                        x_start + genic_bar_width,
                        x_start,
                        x_start,
                    ],
                    "y": [
                        bar_y - bar_height,
                        bar_y - bar_height,
                        bar_y + bar_height,
                        bar_y + bar_height,
                        bar_y - bar_height,
                    ],
                    "fill": "toself",
                    "type": "scatter",
                    "mode": "lines",
                    "line": {"color": bar_border_color, "width": 1},
                    "fillcolor": genic_color,
                    "name": "",
                    "legendgroup": "genic",
                    "legend": "legend2",
                    "xaxis": "x3",
                    "yaxis": "y6",
                    "showlegend": False,
                    "text": f"{mutation} - Genic: {genic_perc:.1f}%",
                    "hoverinfo": "text",
                    "hovertemplate": f"{mutation} - Genic: {genic_perc:.1f}%<extra></extra>",
                }
            )

            intergenic_bar_width = intergenic_prop * bar_width
            intergenic_start = x_start + genic_bar_width
            traces.append(
                {
                    "x": [
                        intergenic_start,
                        intergenic_start + intergenic_bar_width,
                        intergenic_start + intergenic_bar_width,
                        intergenic_start,
                        intergenic_start,
                    ],
                    "y": [
                        bar_y - bar_height,
                        bar_y - bar_height,
                        bar_y + bar_height,
                        bar_y + bar_height,
                        bar_y - bar_height,
                    ],
                    "fill": "toself",
                    "type": "scatter",
                    "mode": "lines",
                    "line": {"color": bar_border_color, "width": 1},
                    "fillcolor": intergenic_color,
                    "name": "",
                    "legendgroup": "genic",
                    "legend": "legend2",
                    "xaxis": "x3",
                    "yaxis": "y6",
                    "showlegend": False,
                    "text": f"{mutation} - Intergenic: {intergenic_perc:.1f}%",
                    "hoverinfo": "text",
                    "hovertemplate": f"{mutation} - Intergenic: {intergenic_perc:.1f}%<extra></extra>",
                }
            )

            if genic_bar_width > 0.02:
                traces.append(
                    {
                        "x": [x_start + genic_bar_width / 2],
                        "y": [bar_y],
                        "mode": "text",
                        "type": "scatter",
                        "text": [f"{genic_perc:.1f}%"],
                        "textfont": {"color": base_text_color, "size": label_font_size},
                        "showlegend": False,
                        "xaxis": "x3",
                        "yaxis": "y6",
                        "hoverinfo": "skip",
                    }
                )

            if intergenic_bar_width > 0.02:
                traces.append(
                    {
                        "x": [intergenic_start + intergenic_bar_width / 2],
                        "y": [bar_y],
                        "mode": "text",
                        "type": "scatter",
                        "text": [f"{intergenic_perc:.1f}%"],
                        "textfont": {"color": pale_text_color, "size": label_font_size},
                        "showlegend": False,
                        "xaxis": "x3",
                        "yaxis": "y6",
                        "hoverinfo": "skip",
                    }
                )

        return traces

"""Python port of LayoutBuilder.js."""


class LayoutBuilder:
    def __init__(self, config):
        self.config = config

    @staticmethod
    def apply_percentage_layout_modifications(
        layout,
        include_strand_bias=False,
        max_strand_bias_percentage=0,
        include_detailed_heatmaps=False,
        use_stacked_bars=False,
    ):
        if layout.get("yaxis"):
            layout["yaxis"] = {
                **layout["yaxis"],
                "title": {
                    "text": "Percentage of Single Base Substitutions",
                    "font": {"family": "Verdana", "size": 10},
                },
                "ticksuffix": "%",
            }

        if include_strand_bias:
            if use_stacked_bars:
                layout["yaxis5"] = {
                    "domain": [0.12, 0.165],
                    "title": {
                        "text": "",
                        "font": {"family": "Courier New, monospace", "size": 10},
                    },
                    "autorange": False,
                    "range": [0.25, 1],
                    "tickmode": "array",
                    "tickvals": [],
                    "ticktext": [],
                    "showline": False,
                    "showticklabels": False,
                    "ticks": "",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                    "showgrid": False,
                    "zeroline": False,
                }

                padding_per_group = 0.015
                bar_width = (1 - 6 * padding_per_group * 2) / 6
                positions = []
                for i in range(6):
                    x_start = (
                        i * (bar_width + 2 * padding_per_group) + padding_per_group
                    )
                    positions.append(x_start + bar_width / 2)

                layout["xaxis2"] = {
                    "domain": [0, 1],
                    "range": [0, 1],
                    "showticklabels": True,
                    "tickmode": "array",
                    "tickvals": positions,
                    "ticktext": ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
                    "showline": False,
                    "ticks": "",
                    "anchor": "y5",
                    "side": "bottom",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                    "showgrid": False,
                    "gridcolor": "#E0E0E0",
                    "gridwidth": 1,
                    "zeroline": False,
                }

                layout["yaxis6"] = {
                    "domain": [0.015, 0.06],
                    "title": {
                        "text": "",
                        "font": {"family": "Courier New, monospace", "size": 10},
                    },
                    "autorange": False,
                    "range": [0.25, 1],
                    "tickmode": "array",
                    "tickvals": [],
                    "ticktext": [],
                    "showline": False,
                    "showticklabels": False,
                    "ticks": "",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                    "showgrid": False,
                    "zeroline": False,
                }

                layout.setdefault("annotations", [])
                layout["annotations"].append(
                    {
                        "text": "<b>Transcribed / Untranscribed strand</b>",
                        "xref": "paper",
                        "yref": "paper",
                        "x": 0.5,
                        "y": 0.155,
                        "xanchor": "center",
                        "yanchor": "bottom",
                        "showarrow": False,
                        "font": {"family": "Verdana", "size": 10, "color": "black"},
                    }
                )
                layout["annotations"].append(
                    {
                        "text": "<b>Genic / Intergenic</b>",
                        "xref": "paper",
                        "yref": "paper",
                        "x": 0.5,
                        "y": 0.05,
                        "xanchor": "center",
                        "yanchor": "bottom",
                        "showarrow": False,
                        "font": {"family": "Verdana", "size": 10, "color": "black"},
                    }
                )

                layout["legend"] = {
                    "x": 1,
                    "xanchor": "left",
                    "y": 0,
                    "yanchor": "bottom",
                    "font": {"size": 10},
                    "bordercolor": "#E0E0E0",
                    "borderwidth": 1,
                }

                layout["xaxis3"] = {
                    "domain": [0, 1],
                    "range": [0, 1],
                    "showticklabels": True,
                    "tickmode": "array",
                    "tickvals": positions,
                    "ticktext": ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
                    "showline": False,
                    "ticks": "",
                    "anchor": "y6",
                    "side": "bottom",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                    "showgrid": False,
                    "gridcolor": "#E0E0E0",
                    "gridwidth": 1,
                    "zeroline": False,
                }
            else:
                layout["yaxis5"] = {
                    "domain": [0, 0.17],
                    "title": {
                        "text": "Strand Bias (%)",
                        "font": {"family": "Courier New, monospace", "size": 10},
                    },
                    "autorange": False,
                    "range": [
                        0,
                        max_strand_bias_percentage + max_strand_bias_percentage * 0.1,
                    ],
                    "linecolor": "#D3D3D3",
                    "linewidth": 1,
                    "mirror": "all",
                    "tickformat": "d",
                    "ticksuffix": "%",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                    "showgrid": True,
                    "gridcolor": "#F5F5F5",
                }
                layout["xaxis2"] = {
                    "domain": [0, 1],
                    "range": [0, 1],
                    "showticklabels": True,
                    "tickmode": "array",
                    "tickvals": [1 / 12, 3 / 12, 5 / 12, 7 / 12, 9 / 12, 11 / 12],
                    "ticktext": ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
                    "showline": True,
                    "linecolor": "#D3D3D3",
                    "linewidth": 1,
                    "mirror": "all",
                    "anchor": "y5",
                    "side": "bottom",
                    "tickfont": {"family": "Courier New, monospace", "size": 8},
                }

        return layout

    def build_base_layout(self, flat_sorted, max_val, annotations, shapes):
        y_max = max_val * 1.2 if max_val > 0 else 1.0
        return {
            "title": "",
            "hoverlabel": {"bgcolor": "#FFF"},
            "height": 800,
            "width": 1080,
            "grid": {"rows": 4, "columns": 1},
            "xaxis": {
                "showticklabels": False,
                "showline": True,
                "tickangle": -90,
                "tickfont": {"family": "Courier New, monospace", "size": 8},
                "tickmode": "array",
                "tickvals": list(range(len(flat_sorted))),
                "ticktext": [e["mutationType"] for e in flat_sorted],
                "linecolor": "#D3D3D3",
                "linewidth": 1,
                "mirror": "all",
                "tickformat": "~s" if max_val > 1000 else "",
                "ticks": "",
            },
            "yaxis": {
                "autorange": False,
                "range": [0, y_max],
                "linecolor": "#D3D3D3",
                "linewidth": 1,
                "mirror": "all",
                "domain": [0.72, 1],
                "tickformat": "~s" if max_val > 1000 else "",
                "tickfont": {"family": "Verdana", "size": 8},
                "showgrid": True,
                "gridcolor": "#F5F5F5",
            },
            "yaxis2": {
                "autorange": True,
                "linecolor": "#D3D3D3",
                "linewidth": 1,
                "ticks": "",
                "mirror": "all",
                "anchor": "x",
                "domain": [0.54, 0.715],
                "tickfont": {"family": "Courier New, monospace", "size": 8},
            },
            "yaxis3": {
                "autorange": True,
                "linecolor": "#D3D3D3",
                "linewidth": 1,
                "ticks": "",
                "mirror": "all",
                "anchor": "x",
                "tickfont": {"family": "Courier New, monospace", "size": 8},
                "domain": [0.36, 0.535],
            },
            "yaxis4": {
                "autorange": True,
                "linecolor": "#D3D3D3",
                "linewidth": 1,
                "ticks": "",
                "mirror": "all",
                "anchor": "x",
                "tickfont": {"family": "Courier New, monospace", "size": 8},
                "domain": [0, 0.35],
            },
            "shapes": shapes,
            "annotations": annotations,
        }

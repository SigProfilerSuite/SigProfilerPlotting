"""Python port of TraceFormatter.js."""


class TraceFormatter:
    @staticmethod
    def apply_percentage_formatting(traces):
        for trace in traces:
            if trace.get("type") == "heatmap" and trace.get("yaxis") == "y4":
                if trace.get("z"):
                    trace["z"] = [[val * 100 for val in row] for row in trace["z"]]
                if trace.get("zmax"):
                    trace["zmax"] = trace["zmax"] * 100
                if (
                    not trace.get("hovertemplate")
                    or "customdata" not in trace["hovertemplate"]
                ):
                    trace["hovertemplate"] = (
                        "x: %{x}<br>y: %{y}<br>Value: %{z:.2g}%<extra></extra>"
                    )
                else:
                    trace["hovertemplate"] = trace["hovertemplate"].replace(
                        ": %{z}", ": %{z:.2g}%"
                    )

            if trace.get("type") == "heatmap" and trace.get("yaxis") in ("y2", "y3"):
                if trace.get("z"):
                    trace["z"] = [[val * 100 for val in row] for row in trace["z"]]
                if trace.get("zmax"):
                    trace["zmax"] = trace["zmax"] * 100
                trace["hovertemplate"] = (
                    "x: %{x}<br>y: %{y}<br>Value: %{z:.2g}%<extra></extra>"
                )

            if trace.get("type") == "bar" and trace.get("yaxis") == "y":
                trace["hovertemplate"] = "%{x}<br>%{y:.2g}%<extra></extra>"

    @staticmethod
    def configure_heatmap_colorbar(traces, include_detailed_heatmaps=False):
        colorbar_configs = (
            {
                "y2": {"len": 0.13, "y": 0.63, "x": 1.02, "show": True},
                "y3": {"len": 0.13, "y": 0.505, "x": 1.02, "show": True},
                "y4": {"len": 0.19, "y": 0.31, "x": 1.02, "show": True},
            }
            if include_detailed_heatmaps
            else {
                "y4": {"len": 0.31, "y": 0.345, "x": 1.02, "show": True},
            }
        )

        shown_colorbars = {}

        for trace in traces:
            if trace.get("type") == "heatmap" and trace.get("yaxis"):
                config = colorbar_configs.get(trace["yaxis"])
                if config:
                    if not shown_colorbars.get(trace["yaxis"]):
                        # apply_percentage_formatting() (called just before
                        # this, sbs4608_plot_builder.py:53-54) already
                        # multiplies every z value by 100 -- e.g. a fraction
                        # of 0.093 becomes 9.3, meaning "9.3%". Plotly's own
                        # ".1%" tickformat would multiply by 100 *again*
                        # (turning 9.3 into "930%"), so use a plain numeric
                        # format plus a literal "%" suffix instead.
                        trace["colorbar"] = {
                            "len": config["len"],
                            "y": config["y"],
                            "x": 1.01,
                            "xanchor": "left",
                            "xpad": 5,
                            "thickness": 20,
                            "tickformat": ".1f",
                            "ticksuffix": "%",
                        }
                        trace["showscale"] = True
                        shown_colorbars[trace["yaxis"]] = True
                    else:
                        trace["showscale"] = False

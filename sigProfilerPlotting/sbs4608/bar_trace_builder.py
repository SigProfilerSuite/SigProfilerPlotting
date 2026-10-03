"""Python port of BarTraceBuilder.js."""


class BarTraceBuilder:
    def __init__(self, colors):
        self.colors = colors

    @staticmethod
    def get_max_y_value(traces):
        ys = [y for trace in traces for y in trace.get("y", [])]
        return max(ys) if ys else 0

    def build_bar_traces(self, grouped_by_mutation):
        items = list(grouped_by_mutation.items())
        traces = []
        for group_index, (mutation, signatures) in enumerate(items):
            prefix_len = sum(len(sigs) for _, sigs in items[:group_index])
            x = [prefix_len + i for i in range(len(signatures))]
            y = [e["contribution"] for e in signatures]
            traces.append(
                {
                    "name": mutation,
                    "type": "bar",
                    "marker": {"color": self.colors[mutation]},
                    "x": x,
                    "y": y,
                    "hoverinfo": "x+y",
                    "showlegend": False,
                }
            )
        return traces

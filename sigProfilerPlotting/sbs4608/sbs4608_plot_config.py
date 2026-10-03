"""Python port of SBS4608PlotConfig.js.

Uses the standard COSMIC SBS96 mutation-type palette (the field convention
used by SigProfilerPlotting / the COSMIC signature browser).
"""

SBS_COLOR = {
    "C>A": "#03BCEE",  # light blue
    "C>G": "#000000",  # black
    "C>T": "#E32926",  # red
    "T>A": "#CAC9C9",  # gray
    "T>C": "#A1CE63",  # green
    "T>G": "#EBC6C4",  # pink
}


def _lighten(hex_color, amount=0.6):
    """Blend a hex color toward white by `amount` (0=no change, 1=white).
    Used to derive the paler bias-bar tint from the same dark base color."""
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    r = round(r + (255 - r) * amount)
    g = round(g + (255 - g) * amount)
    b = round(b + (255 - b) * amount)
    return f"#{r:02x}{g:02x}{b:02x}"


class SBS4608PlotConfig:
    """Centralized configuration with pre-computed layout properties.
    Acts as a dataclass with constants and computed domains. Read-only after construction.
    """

    def __init__(self, include_detailed_heatmaps=False, main_plot_type="default"):
        self.include_detailed_heatmaps = include_detailed_heatmaps
        self.main_plot_type = main_plot_type

        self.plot_height = 1200 if include_detailed_heatmaps else 1050
        self.plot_width = 1080
        self.chunk_size = 16
        self.colors = SBS_COLOR

        self.heatmap_colorscale = [
            [0, "rgb(56,56,156)"],
            [0.2, "rgb(56,56,156)"],
            [0.2, "rgb(106,106,128)"],
            [0.4, "rgb(106,106,128)"],
            [0.4, "rgb(155,146,98)"],
            [0.6, "rgb(155,146,98)"],
            [0.6, "rgb(205,186,69)"],
            [0.8, "rgb(205,186,69)"],
            [0.8, "rgb(255,255,39)"],
            [1, "rgb(255,255,39)"],
        ]

        self.bias_colors = {
            "strand": {
                "transcribed": "#406ae3",
                "untranscribed": "#9acd32",
            },
            "genic": {
                "genic": "#0afaff",
                "intergenic": "#808080",
            },
            "mutations": {
                mutation: {"base": color, "pale": _lighten(color)}
                for mutation, color in SBS_COLOR.items()
            },
        }

        self.heatmap_positions = {
            "front": {"yaxis": "y2", "colorbarY": 0.625, "colorbarLen": 0.2},
            "back": {"yaxis": "y3", "colorbarY": 0.44, "colorbarLen": 0.2},
            "total": {"yaxis": "y4", "colorbarY": 0.17, "colorbarLen": 0.38},
        }

        self._domains = self._compute_domains()

    @property
    def main_plot_domain(self):
        return self._domains["mainPlot"]

    @property
    def front_heatmap_domain(self):
        return self._domains["frontHeatmap"]

    @property
    def back_heatmap_domain(self):
        return self._domains["backHeatmap"]

    @property
    def total_heatmap_domain(self):
        return self._domains["totalHeatmap"]

    @property
    def transcription_bias_domain(self):
        return self._domains["transcriptionBias"]

    @property
    def genic_bias_domain(self):
        return self._domains["genicBias"]

    def _compute_domains(self):
        domains = {}
        if self.include_detailed_heatmaps:
            domains["mainPlot"] = [0.74, 1]
            domains["frontHeatmap"] = [0.57, 0.69]
            domains["backHeatmap"] = [0.445, 0.565]
            domains["totalHeatmap"] = [0.22, 0.44]
        else:
            domains["mainPlot"] = [0.58, 1]
            domains["frontHeatmap"] = None
            domains["backHeatmap"] = None
            domains["totalHeatmap"] = [0.19, 0.50]

        domains["transcriptionBias"] = [0.105, 0.165]
        domains["genicBias"] = [0, 0.06]
        return domains

    def should_include_detailed_heatmaps(self):
        return self.include_detailed_heatmaps

    def is_strand_bias(self):
        return self.main_plot_type == "strand_bias"

    def is_genic_bias(self):
        return self.main_plot_type == "genic_bias"

    def is_default_mode(self):
        return self.main_plot_type == "default"

    def get_colors(self):
        return self.colors

    def get_colorscale(self):
        return self.heatmap_colorscale

    def get_bias_colors(self):
        return self.bias_colors

    def get_heatmap_config(self, type_):
        config = self.heatmap_positions.get(type_)
        if config is None:
            raise ValueError(f"Unknown heatmap type: {type_}")
        return {
            **config,
            "colorscale": self.heatmap_colorscale,
            "chunkSize": self.chunk_size,
        }

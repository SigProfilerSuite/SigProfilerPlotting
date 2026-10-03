"""Python port of SBS4608PlotBuilder.js -- orchestrates plot construction."""

from .heatmap_trace_builder import HeatmapTraceBuilder
from .bar_trace_builder import BarTraceBuilder
from .bias_trace_builder import BiasTraceBuilder
from .annotation_builder import AnnotationBuilder
from .trace_formatter import TraceFormatter
from .tick_label_formatter import TickLabelFormatter
from .layout_builder import LayoutBuilder


class SBS4608PlotBuilder:
    def __init__(self, data, config):
        self.data = data
        self.config = config
        self.heatmap_builder = HeatmapTraceBuilder(config)
        self.bar_builder = BarTraceBuilder(config.get_colors())
        self.layout_builder = LayoutBuilder(config)

    def build(self):
        traces, sections = self._build_all_traces()
        layout = self._build_layout()
        sections["transcription_bias"]["annotation_indices"] = (
            self._find_bias_annotation_indices(layout)
        )
        return {
            "traces": traces,
            "layout": layout,
            "config": {"displayModeBar": False},
            "sections": sections,
        }

    def _build_all_traces(self):
        traces = []
        traces.extend(self._build_main_plot_traces())

        heatmap_start = len(traces)
        total_traces = self.heatmap_builder.build_total_heatmap(
            self.data.group_by_outer_with_tun_proportions(),
            self.data.get_total_mutations(),
        )
        traces.extend(total_traces)
        total_heatmap_range = [heatmap_start, len(traces)]

        detailed_start = len(traces)
        traces.extend(self._build_detailed_heatmap_traces())
        detailed_heatmaps_range = [detailed_start, len(traces)]

        bias_start = len(traces)
        traces.extend(self._build_bias_traces())
        transcription_bias_range = [bias_start, len(traces)]

        TraceFormatter.apply_percentage_formatting(traces)
        TraceFormatter.configure_heatmap_colorbar(
            traces, self.config.should_include_detailed_heatmaps()
        )

        sections = {
            "total_heatmap": {"trace_range": total_heatmap_range},
            "detailed_heatmaps": {"trace_range": detailed_heatmaps_range},
            "transcription_bias": {"trace_range": transcription_bias_range},
        }
        return traces, sections

    def _build_main_plot_traces(self):
        bias_colors = self.config.get_bias_colors()
        if self.data.is_strand_bias():
            return BiasTraceBuilder.build_strand_bias_traces(
                self.data.get_main_plot_data(), bias_colors["strand"]
            )
        elif self.data.is_genic_bias():
            return BiasTraceBuilder.build_genic_bias_traces(
                self.data.get_main_plot_data(), bias_colors["genic"]
            )
        else:
            bar_data = self.data.get_96_bar_chart_data()
            return self.bar_builder.build_bar_traces(bar_data["groupedByMutation"])

    def _build_detailed_heatmap_traces(self):
        traces = []
        if self.config.should_include_detailed_heatmaps():
            detail_heatmap_z_max = self.data.calculate_detail_heatmap_z_max()
            traces.extend(
                self.heatmap_builder.build_front_heatmap(
                    self.data.get_1536_data(),
                    self.data.get_total_mutations(),
                    detail_heatmap_z_max,
                )
            )
            traces.extend(
                self.heatmap_builder.build_back_heatmap(
                    self.data.get_1536_data(),
                    self.data.get_total_mutations(),
                    detail_heatmap_z_max,
                )
            )
        return traces

    def _find_bias_annotation_indices(self, layout):
        wanted = [
            "<b>Transcribed / Untranscribed strand</b>",
            "<b>Genic / Intergenic</b>",
        ]
        indices = []
        for text in wanted:
            for i, ann in enumerate(layout.get("annotations", [])):
                if ann.get("text") == text:
                    indices.append(i)
                    break
        return indices

    def _build_bias_traces(self):
        bias_colors = self.config.get_bias_colors()
        return BiasTraceBuilder.build_stacked_transcription_bias_traces(
            self.data.get_bias_data(),
            bias_colors["mutations"],
        )

    def _build_layout(self):
        bar_data = self.data.get_96_bar_chart_data()
        max_y_value = self._calculate_max_y_value()

        annotations = self._build_annotations(bar_data)
        shapes = self._build_shapes(bar_data)

        layout = self.layout_builder.build_base_layout(
            bar_data["flatSorted"], bar_data["maxValue"], annotations, shapes
        )

        self._configure_axes(layout, max_y_value)

        layout["barmode"] = "group"
        layout["plot_bgcolor"] = "white"
        layout["paper_bgcolor"] = "white"
        layout["height"] = self.config.plot_height
        layout["width"] = self.config.plot_width

        LayoutBuilder.apply_percentage_layout_modifications(
            layout,
            include_strand_bias=True,
            max_strand_bias_percentage=0,
            include_detailed_heatmaps=self.config.should_include_detailed_heatmaps(),
            use_stacked_bars=True,
        )

        TickLabelFormatter.apply_to_layout(
            layout, self.data.get_96_data(), self.config.get_colors()
        )

        layout["annotations"] = AnnotationBuilder.remove_detailed_trinucleotide_labels(
            layout["annotations"]
        )
        AnnotationBuilder.add_missing_ca_label(layout["annotations"])

        return layout

    def _calculate_max_y_value(self):
        if self.data.is_strand_bias() or self.data.is_genic_bias():
            main_traces = self._build_main_plot_traces()
            return BarTraceBuilder.get_max_y_value(main_traces)
        return None

    def _build_annotations(self, bar_data):
        annotations = []
        annotations.extend(
            AnnotationBuilder.build_mutation_annotations(bar_data["groupedByMutation"])
        )

        sample_annotation = AnnotationBuilder.build_sample_annotation(
            self.data.get_1536_data(), "", 0.95
        )
        annotations.append(sample_annotation)

        right_hand_label = self.data.get_right_hand_label()
        if right_hand_label:
            AnnotationBuilder.add_right_hand_label(annotations, right_hand_label)

        return annotations

    def _build_shapes(self, bar_data):
        shapes = []
        shapes.extend(
            AnnotationBuilder.build_shapes(
                bar_data["groupedByMutation"], self.config.get_colors()
            )
        )
        shapes.extend(self._build_vertical_separators(len(bar_data["flatSorted"])))
        return shapes

    def _build_vertical_separators(self, total_mutations):
        n = total_mutations // self.config.chunk_size
        separators = []
        for i in range(n):
            x = (i + 1) * self.config.chunk_size - 0.5
            separators.append(
                {
                    "type": "line",
                    "xref": "x",
                    "yref": "paper",
                    "x0": x,
                    "x1": x,
                    "y0": 0,
                    "y1": 1,
                    "line": {"color": "white", "width": 1},
                }
            )
        return separators

    def _configure_axes(self, layout, max_y_value):
        layout["xaxis"] = {
            **layout["xaxis"],
            "domain": [0, 1],
            "anchor": "y",
            "side": "bottom",
            "showline": True,
            "showticklabels": True,
            "tickangle": -90,
            "linecolor": "#E0E0E0",
            "linewidth": 1,
            "mirror": "all",
        }

        layout["yaxis"] = {
            **layout["yaxis"],
            "domain": self.config.main_plot_domain,
            "tickfont": {"family": "Verdana", "size": 8},
        }
        if max_y_value is not None:
            layout["yaxis"]["range"] = [
                0,
                max_y_value * 1.1 if max_y_value > 0 else 1.0,
            ]

        layout["yaxis4"] = {
            **layout["yaxis4"],
            "domain": self.config.total_heatmap_domain,
            "tickfont": {"family": "Courier New, monospace", "size": 8},
        }

        if self.config.should_include_detailed_heatmaps():
            layout["yaxis2"] = {
                **layout["yaxis2"],
                "domain": self.config.front_heatmap_domain,
                "tickfont": {"family": "Courier New, monospace", "size": 8},
            }
            layout["yaxis3"] = {
                **layout["yaxis3"],
                "domain": self.config.back_heatmap_domain,
                "tickfont": {"family": "Courier New, monospace", "size": 8},
            }

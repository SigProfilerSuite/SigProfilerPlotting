"""Python port of TickLabelFormatter.js."""

from .mutation_data import MutationData
from .annotation_builder import AnnotationBuilder


class TickLabelFormatter:
    @staticmethod
    def generate_x_axis_ticks(data96, colors):
        order = list(colors.keys())

        def mutation_group_sort(group):
            mutation = group["mutation"]
            return order.index(mutation) if mutation in order else len(order)

        grouped_data96 = MutationData.group_by_mutation(
            data96, r"\[(.*)]", mutation_group_sort=mutation_group_sort
        )

        mutation_type_names96 = []
        for group in grouped_data96:
            for e in group["data"]:
                mutation_type_names96.append(
                    {"mutation": group["mutation"], "mutationType": e["mutationType"]}
                )

        tickvals = list(range(len(mutation_type_names96)))
        ticktext = [
            AnnotationBuilder.format_tick_label(
                e["mutation"], e["mutationType"], colors
            )
            for e in mutation_type_names96
        ]
        return tickvals, ticktext

    @staticmethod
    def apply_to_layout(layout, data96, colors):
        tickvals, ticktext = TickLabelFormatter.generate_x_axis_ticks(data96, colors)
        layout["xaxis"] = {
            **layout["xaxis"],
            "tickmode": "array",
            "tickvals": tickvals,
            "ticktext": ticktext,
            "tickfont": {"family": "Courier New, monospace", "color": "#A0A0A0"},
        }

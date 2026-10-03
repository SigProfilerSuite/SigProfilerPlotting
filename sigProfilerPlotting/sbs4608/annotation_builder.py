"""Python port of AnnotationBuilder.js."""

import re
from .mutation_data import MutationData


class AnnotationBuilder:
    @staticmethod
    def build_mutation_annotations(grouped_by_mutation):
        items = list(grouped_by_mutation.items())
        annotations = []
        for group_index, (mutation, signatures) in enumerate(items):
            prefix_len = sum(len(sigs) for _, sigs in items[:group_index])
            annotations.append(
                {
                    "xref": "x",
                    "yref": "paper",
                    "xanchor": "bottom",
                    "yanchor": "bottom",
                    "x": prefix_len + (len(signatures) - 1) * 0.5,
                    "y": 1.04,
                    "text": f"<b>{mutation}</b>",
                    "showarrow": False,
                    "font": {"size": 18},
                    "align": "center",
                }
            )
        return annotations

    @staticmethod
    def build_shapes(grouped_by_mutation, colors):
        items = list(grouped_by_mutation.items())
        shapes = []
        for group_index, (mutation, _) in enumerate(items):
            x0 = sum(len(sigs) for _, sigs in items[:group_index]) - 0.4
            x1 = sum(len(sigs) for _, sigs in items[: group_index + 1]) - 0.6
            shapes.append(
                {
                    "type": "rect",
                    "xref": "x",
                    "yref": "paper",
                    "x0": x0,
                    "x1": x1,
                    "y0": 1.04,
                    "y1": 1.01,
                    "fillcolor": colors[mutation],
                    "line": {"width": 0},
                }
            )
        return shapes

    @staticmethod
    def build_x_axis_annotations(flat_sorted):
        annotations = []
        for index, num in enumerate(flat_sorted):
            annotations.append(
                {
                    "xref": "x",
                    "yref": "paper",
                    "xanchor": "bottom",
                    "yanchor": "bottom",
                    "x": index,
                    "y": -0.04,
                    "text": re.sub(r"\[(.*)\]", "-", num["mutationType"]),
                    "showarrow": False,
                    "font": {"size": 7.5, "family": "Courier New, monospace"},
                    "align": "center",
                    "num": num,
                    "index": index,
                    "textangle": -90,
                }
            )
        return annotations

    @staticmethod
    def build_y_axis_label():
        return {
            "xref": "paper",
            "yref": "paper",
            "xanchor": "top",
            "yanchor": "top",
            "x": -0.045,
            "y": 1.02,
            "text": "<b>Number of Single Base Substitutions</b>",
            "showarrow": False,
            "font": {"size": 10, "family": "Times New Roman"},
            "align": "center",
            "textangle": -90,
        }

    @staticmethod
    def update_y_axis_label(annotations):
        for a in annotations:
            if a.get("text") == "<b>Number of Single Base Substitutions</b>":
                a["text"] = "<b>Percentage of Single Base Substitutions</b>"
                break

    @staticmethod
    def build_sample_annotation(api_data, text="", y_pos=0.88):
        total_mutations = MutationData.get_total_mutations(api_data)
        first = api_data[0]
        sample = first.get("sample")
        if sample and round(total_mutations, 2) > 1:
            label = f"<b>{sample}: {total_mutations:,.0f} {text or ('Indels' if first.get('profile') == 'ID' else 'Substitutions')}</b>"
        elif sample and total_mutations <= 1.1:
            label = f"<b>{sample}</b>"
        else:
            label = f"<b>{first.get('signatureName')}</b>"
        return {
            "xref": "paper",
            "yref": "paper",
            "xanchor": "bottom",
            "yanchor": "bottom",
            "x": 0.01,
            "y": y_pos,
            "text": label,
            "showarrow": False,
            "font": {"size": 24, "family": "Arial"},
            "align": "center",
        }

    @staticmethod
    def format_tick_label(mutation, mutation_type, colors):
        bracket_start = mutation_type.index("[")
        bracket_end = mutation_type.index("]")

        prefix = mutation_type[:bracket_start]
        suffix = mutation_type[bracket_end + 1 :]
        change = mutation_type[bracket_start + 1 : bracket_end]

        if len(prefix) == 3 and len(suffix) == 3:
            context = prefix
        elif len(prefix) == 1 and len(suffix) == 1:
            middle_base = change.split(">")[0]
            context = prefix + middle_base + suffix
        else:
            return mutation_type

        color = colors[change]
        return (
            f"{context[0]}<span style=\"color: {color}; font-family: 'Courier New', monospace;\">"
            f"<b>{context[1]}</b></span>{context[2]}"
        )

    @staticmethod
    def remove_detailed_trinucleotide_labels(annotations):
        result = []
        for a in annotations:
            if (
                a.get("y") == -0.04
                and a.get("text")
                and re.fullmatch(r"\w-\w", a["text"])
            ):
                continue
            result.append(a)
        return result

    @staticmethod
    def add_missing_ca_label(annotations):
        annotations.append(
            {
                "xref": "x",
                "yref": "paper",
                "xanchor": "bottom",
                "yanchor": "bottom",
                "x": 7.5,
                "y": 1.04,
                "text": "<b>C>A</b>",
                "showarrow": False,
                "font": {"size": 18},
                "align": "center",
            }
        )

    @staticmethod
    def add_right_hand_label(annotations, label_text):
        annotations.append(
            {
                "text": f"{label_text}",
                "font": {"size": 18, "family": "Arial"},
                "x": 0.99,
                "xref": "paper",
                "y": 0.98,
                "yref": "paper",
                "xanchor": "right",
                "yanchor": "top",
                "showarrow": False,
            }
        )

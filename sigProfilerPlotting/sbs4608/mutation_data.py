"""Python port of MutationData.js -- line-by-line, same semantics."""

import re
import math


class MutationData:
    def __init__(self, raw_data):
        # deep copy (list of dicts) to avoid mutating caller's data
        self.raw = [dict(e) for e in raw_data]
        self.total_mutations = self._convert_to_percentages()

    @staticmethod
    def convert_to_percentages(data):
        total_mutations = sum(float(e["contribution"]) for e in data)
        for e in data:
            contribution = float(e["contribution"])
            e["mutations"] = (
                (contribution / total_mutations) * 100 if total_mutations else 0.0
            )
            e["originalContribution"] = e["contribution"]
        return total_mutations

    @staticmethod
    def convert_to_percentages_without_normalizing(data):
        for e in data:
            e["mutations"] = float(e["contribution"]) * 100
            e["originalContribution"] = e["contribution"]

    @staticmethod
    def get_total_mutations(api_data):
        total = 0.0
        for e in api_data:
            mutations = float(e.get("mutations") or 0)
            contribution = float(e.get("contribution") or 0)
            total += mutations + contribution
        return total

    @staticmethod
    def group_by_mutation(
        api_data, group_regex, mutation_group_sort=None, mutation_type_sort=None
    ):
        group_by = {}
        for e in api_data:
            mutation = re.search(group_regex, e["mutationType"]).group(1)
            group_by.setdefault(mutation, []).append(e)

        grouped_data = []
        for mutation, data in group_by.items():
            if mutation_type_sort:
                data = sorted(data, key=mutation_type_sort)
            grouped_data.append({"mutation": mutation, "data": data})

        if mutation_group_sort:
            grouped_data.sort(key=mutation_group_sort)
        return grouped_data

    @staticmethod
    def group_by_total(data):
        groups = {}
        for signature in data:
            mutation = re.search(r"\[(.*)\]", signature["mutationType"]).group(1)
            groups.setdefault(mutation, []).append(signature)
        return groups

    @staticmethod
    def chunks(array, size):
        n = math.ceil(len(array) / size) if size else 0
        return [array[i * size : i * size + size] for i in range(n)]

    def _convert_to_percentages(self):
        total_mutations = sum(float(e["contribution"]) for e in self.raw)
        for e in self.raw:
            contribution = float(e["contribution"])
            e["mutations"] = (
                (contribution / total_mutations) * 100 if total_mutations else 0.0
            )
            e["originalContribution"] = e["contribution"]
        return total_mutations

    def get_data(self):
        return self.raw

    def get_total(self):
        return self.total_mutations

    def group_by_outer(self):
        groups = {}
        for entry in self.raw:
            mt = entry["mutationType"]
            mutation = mt[0] + mt[-1]
            signature = {"mutationType": mt, "contribution": entry["mutations"]}
            groups.setdefault(mutation, []).append(signature)
        return groups

    def group_by(self, extractor_fn):
        groups = {}
        for entry in self.raw:
            key = extractor_fn(entry)
            signature = {
                "mutationType": entry["mutationType"],
                "contribution": entry["mutations"],
            }
            groups.setdefault(key, []).append(signature)
        return groups

    def _collapse_to_trinucleotides(self):
        """Collapse 1536-context data to 96 trinucleotide contexts.
        Example: AAAA[C>A]AAAA + CAAA[C>A]AAAA + ... -> AAA[C>A]AAA (sum of 16 variations).
        """
        collapsed = {}
        for entry in self.raw:
            mt = entry["mutationType"]
            trinucleotide = mt[1:-1]
            if trinucleotide not in collapsed:
                collapsed[trinucleotide] = {
                    "mutationType": trinucleotide,
                    "contribution": 0.0,
                    "mutations": 0.0,
                }
            collapsed[trinucleotide]["mutations"] += float(entry["mutations"])
            collapsed[trinucleotide]["contribution"] += float(entry["mutations"])
        return list(collapsed.values())

    def get_bar_chart_data(self):
        data96 = self._collapse_to_trinucleotides()
        grouped_by_mutation = MutationData.group_by_total(data96)
        flat_sorted = [item for group in grouped_by_mutation.values() for item in group]
        max_value = max(o["contribution"] for o in flat_sorted)
        return {
            "groupedByMutation": grouped_by_mutation,
            "flatSorted": flat_sorted,
            "maxValue": max_value,
        }

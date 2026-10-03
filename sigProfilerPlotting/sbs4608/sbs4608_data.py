"""Python port of SBS4608Data.js.

Derives ALL views (1536-context, 96-context bar data, strand/genic bias,
T/U/N hover proportions) from the raw strand-tagged pentanucleotide rows a
signatures/matrix file provides -- nothing is precomputed or invented.
"Genic" = transcribed (T) + untranscribed (U); "Intergenic" = N.
"""

import math
import re
from .mutation_data import MutationData

SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
DINUCS = [a + b for a in "ACGT" for b in "ACGT"]  # 16, alphabetical: AA, AC, ... TT
MUT_RE = re.compile(r"^([TUN]):([ACGT]{2})\[([ACGT]>[ACGT])\]([ACGT]{2})$")


def is_sbs4608_row(mutation_type):
    """True if a mutationType string matches the strand-tagged pentanucleotide
    format this module expects, e.g. 'T:AA[C>A]AA'."""
    return bool(MUT_RE.match(mutation_type))


def load_raw4608(path, column):
    """Read one signature column from a SigProfilerExtractor-style
    *_Signatures.txt (or matrix) file whose mutationType rows are strand-tagged
    pentanucleotide contexts, e.g. 'T:AA[C>A]AA'. Returns
    list of {mutationType: str, contribution: float}."""
    rows = []
    seen = set()
    expected = {
        f"{strand}:{left}[{sub}]{right}"
        for strand in "TUN"
        for sub in SUBS
        for left in DINUCS
        for right in DINUCS
    }
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        if not header or header[0] != "MutationType":
            raise ValueError(
                f"The first column in {path!r} must be named 'MutationType'."
            )
        if column not in header:
            raise ValueError(
                f"Column {column!r} not found in {path!r}. Available columns: "
                f"{', '.join(c for c in header if c != 'MutationType')}"
            )
        col_idx = header.index(column)
        for line_number, line in enumerate(f, start=2):
            parts = line.rstrip("\r\n").split("\t")
            if not parts or not parts[0]:
                continue
            mutation_type = parts[0]
            if mutation_type not in expected:
                raise ValueError(
                    f"Invalid SBS4608 MutationType {mutation_type!r} at line "
                    f"{line_number} in {path!r}. Expected a T/U/N-tagged "
                    "pentanucleotide context such as 'T:AA[C>A]AA'."
                )
            if mutation_type in seen:
                raise ValueError(
                    f"Duplicate SBS4608 MutationType {mutation_type!r} at line "
                    f"{line_number} in {path!r}."
                )
            if len(parts) <= col_idx:
                raise ValueError(
                    f"Missing value for column {column!r} at line {line_number} "
                    f"in {path!r}."
                )
            try:
                contribution = float(parts[col_idx])
            except ValueError as error:
                raise ValueError(
                    f"Non-numeric value {parts[col_idx]!r} for column {column!r} "
                    f"at line {line_number} in {path!r}."
                ) from error
            if not math.isfinite(contribution) or contribution < 0:
                raise ValueError(
                    f"SBS4608 contribution must be finite and non-negative; got "
                    f"{parts[col_idx]!r} at line {line_number} in {path!r}."
                )
            seen.add(mutation_type)
            rows.append({"mutationType": mutation_type, "contribution": contribution})
    missing = expected - seen
    if missing:
        examples = ", ".join(sorted(missing)[:3])
        raise ValueError(
            f"Incomplete SBS4608 matrix in {path!r}: found {len(seen)} of 4608 "
            f"required MutationType rows. Missing examples: {examples}."
        )
    return rows


class SBS4608Data:
    def __init__(self, raw4608, signature_name, main_plot_type="default"):
        self.main_plot_type = main_plot_type
        self.signature_name = signature_name
        self.raw4608 = raw4608
        self.right_hand_label = None

        self.raw1536 = self._build_raw1536()
        self.raw96 = self._build_raw96()

        self._mutation_data = MutationData(self.raw1536)
        self._bar_chart_data = self._mutation_data.get_bar_chart_data()

        self._prepared_bias_data = None
        if main_plot_type == "strand_bias":
            self._prepared_bias_data = self._prepare_bias_data(
                self._build_raw192_strand()
            )
        elif main_plot_type == "genic_bias":
            self._prepared_bias_data = self._prepare_bias_data(
                self._build_raw192_genic()
            )

    def _build_raw1536(self):
        """Sum over strand (T+U+N) for each pentanucleotide context. Iteration
        order is sub-outer, left-outer, right-inner -- required for the
        front/back detailed-heatmap Y-axis label sampling in
        HeatmapTraceBuilder to land on the correct nucleotide blocks."""
        sums = {}
        for row in self.raw4608:
            m = MUT_RE.match(row["mutationType"])
            _, left, sub, right = m.groups()
            key = (sub, left, right)
            sums[key] = sums.get(key, 0.0) + row["contribution"]

        rows = []
        for sub in SUBS:
            for left in DINUCS:
                for right in DINUCS:
                    contribution = sums.get((sub, left, right), 0.0)
                    rows.append(
                        {
                            "mutationType": f"{left}[{sub}]{right}",
                            "contribution": contribution,
                            "signatureName": self.signature_name,
                        }
                    )
        return rows

    def _build_raw96(self):
        """Trinucleotide-context labels only (used solely for X-axis tick
        generation, which needs mutationType strings, not values)."""
        seen = set()
        rows = []
        for sub in SUBS:
            for l1 in "ACGT":
                for r1 in "ACGT":
                    tri = f"{l1}[{sub}]{r1}"
                    if tri not in seen:
                        seen.add(tri)
                        rows.append({"mutationType": tri})
        return rows

    def _build_raw192_strand(self):
        """96-context x {T,U} bias data, contribution = raw (unnormalized) sum
        over the dropped dimension (the 5'/3' outer bases), matching the
        96-trinucleotide resolution BiasTraceBuilder expects."""
        sums = {}
        for row in self.raw4608:
            m = MUT_RE.match(row["mutationType"])
            strand, left, sub, right = m.groups()
            if strand not in ("T", "U"):
                continue
            context96 = f"{left[1]}[{sub}]{right[0]}"
            key = (strand, context96)
            sums[key] = sums.get(key, 0.0) + row["contribution"]
        return [
            {"mutationType": f"{strand}:{ctx}", "contribution": val}
            for (strand, ctx), val in sums.items()
        ]

    def _build_raw192_genic(self):
        """96-context x {Genic, Intergenic} bias data. Genic = T+U, Intergenic = N."""
        sums = {}
        for row in self.raw4608:
            m = MUT_RE.match(row["mutationType"])
            strand, left, sub, right = m.groups()
            status = "Genic" if strand in ("T", "U") else "Intergenic"
            context96 = f"{left[1]}[{sub}]{right[0]}"
            key = (status, context96)
            sums[key] = sums.get(key, 0.0) + row["contribution"]
        return [
            {"mutationType": f"{status}:{ctx}", "contribution": val}
            for (status, ctx), val in sums.items()
        ]

    def _prepare_bias_data(self, raw_data):
        data = [dict(e) for e in raw_data]
        total = sum(float(entry["contribution"]) for entry in self.raw4608)
        for entry in data:
            contribution = float(entry["contribution"])
            entry["mutations"] = (contribution / total) * 100 if total else 0.0
            entry["originalContribution"] = entry["contribution"]
        return data

    def get_1536_data(self):
        return self._mutation_data.get_data()

    def get_96_data(self):
        return self.raw96

    def get_96_bar_chart_data(self):
        return self._bar_chart_data

    def get_main_plot_data(self):
        if self.main_plot_type in ("strand_bias", "genic_bias"):
            return self._prepared_bias_data
        return self._bar_chart_data["groupedByMutation"]

    def get_bias_data(self):
        """18 = 6 mutation types x 3 strands (T/U/N), each the RAW (non-percentage)
        sum of contributions across all contexts. Feeds the bottom stacked bias bars."""
        sums = {}
        for row in self.raw4608:
            m = MUT_RE.match(row["mutationType"])
            strand, left, sub, right = m.groups()
            key = (strand, sub)
            sums[key] = sums.get(key, 0.0) + row["contribution"]
        return [
            {"mutationType": f"{strand}:XX[{sub}]XX", "contribution": val}
            for (strand, sub), val in sums.items()
        ]

    def group_by_outer(self):
        return self._mutation_data.group_by_outer()

    def get_total_mutations(self):
        return self._mutation_data.get_total()

    def _build_tun_map(self):
        tun_map = {}
        for entry in self.raw4608:
            mt = entry["mutationType"]
            prefix = mt[0]
            context1536 = mt[2:]
            tun_map.setdefault(context1536, {"T": 0.0, "U": 0.0, "N": 0.0})
            tun_map[context1536][prefix] = entry["contribution"]
        return tun_map

    def group_by_outer_with_tun_proportions(self):
        tun_map = self._build_tun_map()
        grouped = self._mutation_data.group_by_outer()

        enhanced = {}
        for outer_pair, signatures in grouped.items():
            enhanced_sigs = []
            for sig in signatures:
                context1536 = sig["mutationType"]
                tun_data = tun_map.get(context1536, {"T": 0.0, "U": 0.0, "N": 0.0})
                total = tun_data["T"] + tun_data["U"] + tun_data["N"]
                enhanced_sigs.append(
                    {
                        "mutationType": sig["mutationType"],
                        "contribution": sig["contribution"],
                        "tunProportions": {
                            "T": (tun_data["T"] / total) if total > 0 else 0,
                            "U": (tun_data["U"] / total) if total > 0 else 0,
                            "N": (tun_data["N"] / total) if total > 0 else 0,
                        },
                    }
                )
            enhanced[outer_pair] = enhanced_sigs
        return enhanced

    def get_max_bar_value(self):
        return self._bar_chart_data["maxValue"]

    def get_right_hand_label(self):
        return self.right_hand_label

    def is_strand_bias(self):
        return self.main_plot_type == "strand_bias"

    def is_genic_bias(self):
        return self.main_plot_type == "genic_bias"

    def is_default_mode(self):
        return self.main_plot_type == "default"

    def calculate_detail_heatmap_z_max(self):
        data = self.get_1536_data()

        def calculate_max_z(group_by, regroup_by, sort_grouped=None):
            grouped = {}
            for entry in data:
                key = group_by(entry)
                signature = {
                    "mutationType": entry["mutationType"],
                    "contribution": entry["mutations"],
                }
                grouped.setdefault(key, []).append(signature)

            summed = [
                {
                    "mutationType": group[0]["mutationType"],
                    "contribution": sum(e["contribution"] for e in group),
                }
                for group in grouped.values()
            ]

            if sort_grouped:
                summed.sort(key=sort_grouped)

            regrouped = {}
            for entry in summed:
                mutation = regroup_by(entry)
                regrouped.setdefault(mutation, []).append(entry)

            z_values = [
                [e["contribution"] / 100 for e in group] for group in regrouped.values()
            ]
            flat = [v for group in z_values for v in group]
            return max(flat) if flat else 0

        front_max_z = calculate_max_z(
            lambda e: e["mutationType"][:-1],
            lambda e: re.search(r"\[(.*)\]", e["mutationType"]).group(1),
        )
        back_max_z = calculate_max_z(
            lambda e: e["mutationType"][1:],
            lambda e: re.search(r"\[(.*)\]", e["mutationType"]).group(1),
            sort_grouped=lambda e: e["mutationType"][-1],
        )
        return max(front_max_z, back_max_z)

"""Python port of HeatmapTraceBuilder.js."""

import re
import math


class HeatmapTraceBuilder:
    def __init__(self, config):
        self.config = config

    def build_front_heatmap(self, data, total_mutations, shared_z_max=None):
        front_config = self.config.get_heatmap_config("front")

        def build_y_labels(y_labels):
            return [
                y_labels[0][0][0] + "--N",
                y_labels[0][16][0] + "--N",
                y_labels[0][32][0] + "--N",
                y_labels[0][48][0] + "--N",
            ]

        return self.build_heatmap(
            data=data,
            group_by=lambda e: e["mutationType"][:-1],
            regroup_by=lambda e: re.search(r"\[(.*)\]", e["mutationType"]).group(1),
            build_y_labels=build_y_labels,
            yaxis=front_config["yaxis"],
            colorbar_y=front_config["colorbarY"],
            colorbar_len=front_config["colorbarLen"],
            total_mutations=total_mutations,
            chunk_size=front_config["chunkSize"],
            colorscale=front_config["colorscale"],
            shared_z_max=shared_z_max,
        )

    def build_back_heatmap(self, data, total_mutations, shared_z_max=None):
        back_config = self.config.get_heatmap_config("back")

        def build_y_labels(y_labels):
            return [
                "N--" + y_labels[0][0][-1],
                "N--" + y_labels[0][16][-1],
                "N--" + y_labels[0][32][-1],
                "N--" + y_labels[0][48][-1],
            ]

        def sort_grouped(entry):
            return entry["mutationType"][-1]

        return self.build_heatmap(
            data=data,
            group_by=lambda e: e["mutationType"][1:],
            regroup_by=lambda e: re.search(r"\[(.*)\]", e["mutationType"]).group(1),
            sort_grouped=sort_grouped,
            build_y_labels=build_y_labels,
            yaxis=back_config["yaxis"],
            colorbar_y=back_config["colorbarY"],
            colorbar_len=back_config["colorbarLen"],
            total_mutations=total_mutations,
            chunk_size=back_config["chunkSize"],
            colorscale=back_config["colorscale"],
            shared_z_max=shared_z_max,
        )

    def build_total_heatmap(self, group_by_mutation_outer, total_mutations):
        heatmap_y = []
        heatmap_z = []
        heatmap_customdata = []

        for key, value in group_by_mutation_outer.items():
            value = sorted(value, key=lambda v: v["mutationType"][3:6])

            heatmap_y.append(key[0] + "--" + key[-1])
            heatmap_z.append([v["contribution"] / 100 for v in value])
            heatmap_customdata.append(
                [
                    {
                        "T": v["tunProportions"]["T"],
                        "U": v["tunProportions"]["U"],
                        "N": v["tunProportions"]["N"],
                        "context": v["mutationType"],
                    }
                    for v in value
                ]
            )

        # Transpose: 16 rows of 96 items -> 6 columns of 16x16 chunks
        heat_map_z_final = [
            [row[chunk_index * 16 : (chunk_index + 1) * 16] for row in heatmap_z]
            for chunk_index in range(6)
        ]
        heatmap_customdata_final = [
            [
                row[chunk_index * 16 : (chunk_index + 1) * 16]
                for row in heatmap_customdata
            ]
            for chunk_index in range(6)
        ]

        total_config = self.config.get_heatmap_config("total")
        max_z = max(v for chunk in heat_map_z_final for row in chunk for v in row)
        zmax = max_z * 1.1 if max_z > 0 else 1.0

        traces = []
        for index, num in enumerate(heat_map_z_final):
            # JS computed this via `array.slice(0,index).reduce((x0,[_,sigs])=>x0+sigs.length,0) + i`;
            # since every prior chunk's destructured "sigs.length" is 16, this reduces to index*16 + i.
            x_values = [index * 16 + i for i in range(len(num[0]))] if num else []
            traces.append(
                {
                    "colorbar": {
                        "len": total_config["colorbarLen"],
                        "y": total_config["colorbarY"],
                    },
                    "colorscale": total_config["colorscale"],
                    "zmin": 0,
                    "zmax": zmax,
                    "z": num,
                    "y": heatmap_y,
                    "type": "heatmap",
                    "hoverongaps": False,
                    "xaxis": "x",
                    "yaxis": total_config["yaxis"],
                    "x": x_values,
                    "xgap": 0.1,
                    "ygap": 0.1,
                    "customdata": heatmap_customdata_final[index],
                    "hovertemplate": (
                        "%{customdata.context}: %{z}<br>"
                        "T: %{customdata.T:.0%}<br>"
                        "U: %{customdata.U:.0%}<br>"
                        "N: %{customdata.N:.0%}<extra></extra>"
                    ),
                }
            )

        return traces

    def build_heatmap(
        self,
        data,
        group_by,
        regroup_by,
        build_y_labels=None,
        sort_grouped=None,
        yaxis=None,
        colorbar_y=None,
        colorbar_len=0.2,
        total_mutations=None,
        chunk_size=16,
        colorscale=None,
        shared_z_max=None,
    ):
        grouped = self._group_data(data, group_by)
        summed = self._sum_contributions(grouped)
        if sort_grouped:
            summed.sort(key=sort_grouped)
        regrouped = self._regroup_by_mutation(summed, regroup_by)
        y_labels, z_values = self._build_yz(regrouped, total_mutations)
        return self._create_traces(
            z_values=z_values,
            y_labels=y_labels,
            yaxis=yaxis,
            colorbar_y=colorbar_y,
            colorbar_len=colorbar_len,
            chunk_size=chunk_size,
            colorscale=colorscale,
            build_y_labels_fn=build_y_labels,
            shared_z_max=shared_z_max,
        )

    def _group_data(self, data, group_by_fn):
        groups = {}
        for entry in data:
            key = group_by_fn(entry)
            signature = {
                "mutationType": entry["mutationType"],
                "contribution": entry.get("mutations", entry["contribution"]),
            }
            groups.setdefault(key, []).append(signature)
        return groups

    def _sum_contributions(self, grouped):
        return [
            {
                "mutationType": key,
                "contribution": sum(float(e["contribution"]) for e in value),
            }
            for key, value in grouped.items()
        ]

    def _regroup_by_mutation(self, summed, regroup_by_fn):
        groups = {}
        for entry in summed:
            mutation = regroup_by_fn(entry)
            signature = {
                "mutationType": entry["mutationType"],
                "contribution": entry["contribution"],
            }
            groups.setdefault(mutation, []).append(signature)
        return groups

    def _build_yz(self, regrouped, total_mutations):
        y_arrays = []
        z_arrays = []
        for key, value in regrouped.items():
            y_arrays.append([v["mutationType"] for v in value])
            z_arrays.append([v["contribution"] / 100 for v in value])
        return y_arrays, z_arrays

    def _create_traces(
        self,
        z_values,
        y_labels,
        yaxis,
        colorbar_y,
        colorbar_len,
        chunk_size,
        colorscale,
        build_y_labels_fn=None,
        shared_z_max=None,
    ):
        def chunks(arr, size):
            n = math.ceil(len(arr) / size) if size else 0
            return [arr[i * size : i * size + size] for i in range(n)]

        if build_y_labels_fn:
            y_labels_final = build_y_labels_fn(y_labels)
        else:
            y_labels_final = [
                y_labels[0][i * chunk_size]
                for i in range(math.ceil(len(y_labels[0]) / chunk_size))
            ]

        z_chunked = [chunks(z, chunk_size) for z in z_values]

        flat_all = [v for group in z_chunked for chunk in group for v in chunk]
        max_z = (
            shared_z_max
            if shared_z_max is not None
            else (max(flat_all) if flat_all else 0)
        )
        zmax = max_z * 1.1 if max_z > 0 else 1.0

        traces = []
        for mutation_group_index, chunked_array in enumerate(z_chunked):
            start_x = mutation_group_index * chunk_size
            x_values = [
                start_x + chunk_index for chunk_index in range(len(chunked_array))
            ]
            traces.append(
                {
                    "colorbar": {"len": colorbar_len, "y": colorbar_y},
                    "colorscale": colorscale,
                    "zmin": 0,
                    "zmax": zmax,
                    "z": chunked_array,
                    "y": y_labels_final,
                    "type": "heatmap",
                    "hoverongaps": False,
                    "xaxis": "x",
                    "yaxis": yaxis,
                    "x": x_values,
                    "xgap": 0.1,
                    "ygap": 0.1,
                    "hovertemplate": "x: %{x}<br>y: %{y}<br>Value: %{z}<extra></extra>",
                }
            )
        return traces

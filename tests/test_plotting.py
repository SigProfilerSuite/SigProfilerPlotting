import os
import numpy as np
from PIL import Image
import sigProfilerPlotting as sigPlt
import pytest
import pandas as pd

from tests.image_comparison import image_similarities

current_script_path = os.path.abspath(__file__)

SPP_PATH = os.path.dirname(current_script_path)
def plotSV_wrapper(
    matrix_path, output_path, project, context, savefig_format="png", **kwargs
):
    # Call the actual plotSV function with the correct parameters
    sigPlt.plotSV(
        matrix_path=matrix_path,
        output_path=output_path,
        project=project,
        savefig_format=savefig_format,
        percentage=kwargs.get("percentage", False),
        aggregate=kwargs.get("aggregate", False),
        dpi=kwargs.get("dpi", 100),
    )


def plotCNV_wrapper(
    matrix_path, output_path, project, context, savefig_format="png", **kwargs
):
    if type(matrix_path) == str:
        read_from_file = True
    else:
        read_from_file = False
    # Call the actual plotCNV function with the correct parameters
    sigPlt.plotCNV(
        matrix_path=matrix_path,
        output_path=output_path,
        project=project,
        savefig_format=savefig_format,
        read_from_file=read_from_file,
        percentage=kwargs.get("percentage", False),
        aggregate=kwargs.get("aggregate", False),
        dpi=kwargs.get("dpi", 100),
    )


test_configs = {
    "SBS96": {
        "type": "SBS",
        "context": "96",
        "function": sigPlt.plotSBS,
        "example_file": "example.SBS96.all",
    },
    "SBS288": {
        "type": "SBS",
        "context": "288",
        "function": sigPlt.plotSBS,
        "example_file": "example.SBS288.all",
    },
    "DBS78": {
        "type": "DBS",
        "context": "78",
        "function": sigPlt.plotDBS,
        "example_file": "example.DBS78.all",
    },
    "ID83": {
        "type": "ID",
        "context": "83",
        "function": sigPlt.plotID,
        "example_file": "example.ID83.all",
    },
    "CNV48": {
        "type": "CNV",
        "context": "48",
        "function": plotCNV_wrapper,
        "example_file": "example.CNV48.tsv",
    },
    "SV32": {
        "type": "SV",
        "context": "32",
        "function": plotSV_wrapper,
        "example_file": "example.SV32.tsv",
    },
}


# Compare two input routes rendered on the same host. A static PNG is not a
# reliable oracle across OS font rasterizers and Matplotlib versions.
@pytest.mark.parametrize("config_key", test_configs.keys())
def test_plot_generation(config_key, tmp_path):
    config = test_configs[config_key]
    matrix_path = os.path.join(
        SPP_PATH, "input", config["type"], "unordered", config["example_file"]
    )
    matrix_frame = pd.read_csv(matrix_path, sep="\t")
    outputs = []

    for input_name, matrix in (("file", matrix_path), ("dataframe", matrix_frame)):
        output_directory = tmp_path / input_name
        output_directory.mkdir()
        config["function"](
            matrix,
            str(output_directory),
            "test",
            config["context"],
            savefig_format="png",
            percentage=False,
        )
        output = output_directory / (
            f"{config['type']}_{config['context']}_plots_Random.png"
        )
        assert output.is_file() and output.stat().st_size > 0
        outputs.append(output)

    structure, chroma = image_similarities(outputs[0], outputs[1])
    assert structure >= 0.99, f"{config_key} file/DataFrame SSIM={structure:.6f}"
    assert chroma >= 0.99, f"{config_key} file/DataFrame chroma={chroma:.6f}"

    with Image.open(outputs[0]) as file_image, Image.open(outputs[1]) as frame_image:
        assert file_image.size == frame_image.size
        assert file_image.width > 100 and file_image.height > 100
        pixels = np.asarray(file_image.convert("RGB"))
        height, width = pixels.shape[:2]
        # Exclude the colored headings and x-axis labels: this region must
        # contain actual bars, not just a well-formed but empty plot frame.
        chart = pixels[
            int(0.20 * height) : int(0.85 * height),
            int(0.15 * width) : int(0.95 * width),
        ]
        assert np.mean(chart.mean(axis=2) < 245) > 0.025
        assert np.mean(chart.max(axis=2) - chart.min(axis=2) > 10) > 0.02

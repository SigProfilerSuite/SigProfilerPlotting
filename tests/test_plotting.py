import os
from PIL import Image
import sigProfilerPlotting as sigPlt
import pytest
import pandas as pd

from tests.image_comparison import (
    MINIMUM_CHROMA_SIMILARITY,
    MINIMUM_STRUCTURAL_SIMILARITY,
    image_similarities,
)

current_script_path = os.path.abspath(__file__)

SPP_PATH = os.path.dirname(current_script_path)
SPP_SBS = os.path.join(SPP_PATH, "input", "SBS")
SPP_DBS = os.path.join(SPP_PATH, "input", "DBS")
SPP_ID = os.path.join(SPP_PATH, "input", "ID")
SPP_CNV = os.path.join(SPP_PATH, "input", "CNV")
SPP_SV = os.path.join(SPP_PATH, "input", "SV")
SPP_STANDARD_PNG = os.path.join(SPP_PATH, "standard_png")


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
        "crop_dimensions": {
            "Random": None,
            "xaxis": (170, 905, 6000, 2000),
            "bars": (170, 300, 4340, 910),
        },
    },
    "SBS288": {
        "type": "SBS",
        "context": "288",
        "function": sigPlt.plotSBS,
        "example_file": "example.SBS288.all",
        "crop_dimensions": {
            "Random": None,
            "xaxis": (170, 905, 3250, 2000),
            "bars": (170, 300, 3250, 910),
        },
    },
    "DBS78": {
        "type": "DBS",
        "context": "78",
        "function": sigPlt.plotDBS,
        "example_file": "example.DBS78.all",
        "crop_dimensions": {
            "Random": None,
            "xaxis": (170, 905, 4350, 990),
            "bars": (170, 280, 4350, 910),
        },
    },
    "ID83": {
        "type": "ID",
        "context": "83",
        "function": sigPlt.plotID,
        "example_file": "example.ID83.all",
        "crop_dimensions": {
            "Random": None,
            "xaxis": (170, 1063, 4250, 1100),
            "bars": (190, 350, 4250, 1003),
        },
    },
    "CNV48": {
        "type": "CNV",
        "context": "48",
        "function": plotCNV_wrapper,
        "example_file": "example.CNV48.tsv",
        "crop_dimensions": {
            "Random": None,
            "xaxis": (0, 878, 1372, 1021),
            "bars": (120, 95, 1372, 878),
        },
    },
    "SV32": {
        "type": "SV",
        "context": "32",
        "function": plotSV_wrapper,
        "example_file": "example.SV32.tsv",
        "crop_dimensions": {
            "Random": None,
            "xaxis": (0, 740, 1360, 865),
            "bars": (120, 120, 1340, 725),
        },
    },
}


@pytest.fixture
def config(request):
    config_key = request.param
    return test_configs[config_key]


@pytest.fixture(params=["file", "dataframe"])
def input_data(request):
    def _input_data(config):
        example_file_path = os.path.join(
            SPP_PATH, "input", config["type"], "unordered", config["example_file"]
        )
        if request.param == "dataframe":
            df = pd.read_csv(example_file_path, sep="\t")
            return df
        else:
            return example_file_path

    return _input_data


# Main test function tests plots for SBS, DBS, ID, CNV, and SV for dataframes and files
@pytest.mark.parametrize("config_key", test_configs.keys())
def test_plot_generation(config_key, input_data):
    config = test_configs[config_key]
    mutation_type = config["type"]
    mutation_path = os.path.join(SPP_PATH, "input", mutation_type)
    output_subdir = os.path.join(mutation_path, "output")
    output_directory = os.path.join(output_subdir, f"{config_key}_full_image{os.sep}")

    os.makedirs(output_directory, exist_ok=True)

    data = input_data(config)  # Pass config to input_data fixture

    # Modify the wrapper functions to accept data and input_type
    # and call the plotting function accordingly
    config["function"](
        data,
        output_directory,
        "test",
        config["context"],
        savefig_format="png",
        percentage=False,
    )

    for test_case, crop_area in config["crop_dimensions"].items():
        test_image_path = os.path.join(
            output_directory, f"{config['type']}_{config['context']}_plots_Random.png"
        )
        cropped_test_image_path = os.path.join(
            output_directory,
            f"{config['type']}_{config['context']}_plots_{test_case}.png",
        )
        standard_image_name = (
            f"{config['type']}_{config['context']}_plots_{test_case}.png"
        )
        standard_image_path = os.path.join(SPP_STANDARD_PNG, standard_image_name)

        if crop_area:
            with Image.open(test_image_path) as img:
                img = img.crop(crop_area)
                img.save(cropped_test_image_path)
        else:
            # If there's no crop_area, use the original test image for comparison
            # This handles the case for comparing the entire plot
            os.rename(test_image_path, cropped_test_image_path)

        structure, chroma = image_similarities(
            cropped_test_image_path, standard_image_path
        )
        assert structure >= MINIMUM_STRUCTURAL_SIMILARITY, (
            f"Images for {config_key}, {test_case} have a structural mismatch: "
            f"SSIM={structure:.6f}, minimum={MINIMUM_STRUCTURAL_SIMILARITY:.6f}."
        )
        assert chroma >= MINIMUM_CHROMA_SIMILARITY, (
            f"Images for {config_key}, {test_case} have a color mismatch: "
            f"chroma={chroma:.6f}, minimum={MINIMUM_CHROMA_SIMILARITY:.6f}."
        )

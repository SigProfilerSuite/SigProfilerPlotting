import pandas as pd
import sigProfilerPlotting as sigPlt
import pytest
import os
from sigProfilerPlotting import process_input, get_context_reference


# Path to the tests directory
SPP_TEST_PATH = os.path.dirname(os.path.abspath(__file__))

SPP_SBS = os.path.join(SPP_TEST_PATH, "input/SBS/")
SPP_DBS = os.path.join(SPP_TEST_PATH, "input/DBS/")
SPP_ID = os.path.join(SPP_TEST_PATH, "input/ID/")
SPP_CNV = os.path.join(SPP_TEST_PATH, "input/CNV/")
SPP_SV = os.path.join(SPP_TEST_PATH, "input/SV/")

# List of plot types to test
plot_types_SBS = ["96", "288"]
plot_types_DBS = ["78"]
plot_types_ID = ["83"]
plot_types_CNV = ["48"]
plot_types_SV = ["32"]

path_options = ["unordered", "ordered", "wrong_rows"]


################SBS#############
@pytest.mark.parametrize("plot_type", plot_types_SBS)
@pytest.mark.parametrize("path", path_options)
def test_get_context_reference_SBS(plot_type, path):
    file_path = os.path.join(SPP_SBS, path, f"example.SBS{plot_type}.all")

    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)
    input_data_index = input_data.iloc[:, 0].tolist()
    expected_index = get_context_reference(plot_type)

    if path == "ordered":
        assert input_data_index == expected_index
    else:
        assert input_data_index != expected_index


@pytest.mark.parametrize("plot_type", plot_types_SBS)
@pytest.mark.parametrize("path", path_options)
def test_process_input_SBS(plot_type, path):
    file_path = os.path.join(SPP_SBS, path, f"example.SBS{plot_type}.all")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    if path == "wrong_rows":
        with pytest.raises(
            ValueError, match=f"Input matrix file should have {plot_type} rows"
        ):
            process_input(input_data, plot_type)
    else:
        ordered_input_data = process_input(input_data, plot_type)
        ordered_input_data_index = ordered_input_data.index.tolist()
        expected_index = get_context_reference(plot_type)
        assert ordered_input_data_index == expected_index


################DBS#############
@pytest.mark.parametrize("plot_type", plot_types_DBS)
@pytest.mark.parametrize("path", path_options)
def test_get_context_reference_DBS(plot_type, path):
    file_path = os.path.join(SPP_DBS, path, f"example.DBS{plot_type}.all")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    input_data_index = input_data.iloc[:, 0].tolist()
    expected_index = get_context_reference(plot_type)

    if path == "ordered":
        assert input_data_index == expected_index
    else:
        assert input_data_index != expected_index


@pytest.mark.parametrize("plot_type", plot_types_DBS)
@pytest.mark.parametrize("path", path_options)
def test_process_input_DBS(plot_type, path):
    file_path = os.path.join(SPP_DBS, path, f"example.DBS{plot_type}.all")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    if path == "wrong_rows":
        with pytest.raises(
            ValueError, match=f"Input matrix file should have {plot_type} rows"
        ):
            process_input(input_data, plot_type)
    else:
        ordered_input_data = process_input(input_data, plot_type)
        ordered_input_data_index = ordered_input_data.index.tolist()
        expected_index = get_context_reference(plot_type)
        assert ordered_input_data_index == expected_index


################ID#############
@pytest.mark.parametrize("plot_type", plot_types_ID)
@pytest.mark.parametrize("path", path_options)
def test_get_context_reference_ID(plot_type, path):
    file_path = os.path.join(SPP_ID, path, f"example.ID{plot_type}.all")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    input_data_index = input_data.iloc[:, 0].tolist()
    expected_index = get_context_reference(plot_type)

    if path == "ordered":
        assert input_data_index == expected_index
    else:
        assert input_data_index != expected_index


@pytest.mark.parametrize("plot_type", plot_types_ID)
@pytest.mark.parametrize("path", path_options)
def test_process_input_ID(plot_type, path):
    file_path = os.path.join(SPP_ID, path, f"example.ID{plot_type}.all")

    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    if path == "wrong_rows":
        with pytest.raises(
            ValueError, match=f"Input matrix file should have {plot_type} rows"
        ):
            process_input(input_data, plot_type)
    else:
        ordered_input_data = process_input(input_data, plot_type)
        ordered_input_data_index = ordered_input_data.index.tolist()
        expected_index = get_context_reference(plot_type)
        assert ordered_input_data_index == expected_index


################CNV#############
@pytest.mark.parametrize("plot_type", plot_types_CNV)
@pytest.mark.parametrize("path", path_options)
def test_get_context_reference_CNV(plot_type, path):
    file_path = os.path.join(SPP_CNV, path, f"example.CNV{plot_type}.tsv")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)
    input_data_index = input_data.iloc[:, 0].tolist()
    expected_index = get_context_reference(plot_type)

    if path == "ordered":
        assert input_data_index == expected_index
    else:
        assert input_data_index != expected_index


@pytest.mark.parametrize("plot_type", plot_types_CNV)
@pytest.mark.parametrize("path", path_options)
def test_process_input_CNV(plot_type, path):
    file_path = os.path.join(SPP_CNV, path, f"example.CNV{plot_type}.tsv")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    if path == "wrong_rows":
        with pytest.raises(
            ValueError, match=f"Input matrix file should have {plot_type} rows"
        ):
            process_input(input_data, plot_type)
    else:
        ordered_input_data = process_input(input_data, plot_type)
        ordered_input_data_index = ordered_input_data.index.tolist()
        expected_index = get_context_reference(plot_type)
        assert ordered_input_data_index == expected_index


# ###############SV#############
@pytest.mark.parametrize("plot_type", plot_types_SV)
@pytest.mark.parametrize("path", path_options)
def test_get_context_reference_SV(plot_type, path):
    file_path = os.path.join(SPP_SV, path, f"example.SV{plot_type}.tsv")

    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)
    input_data_index = input_data.iloc[:, 0].tolist()
    expected_index = get_context_reference(plot_type)

    if path == "ordered":
        assert input_data_index == expected_index
    else:
        assert input_data_index != expected_index


@pytest.mark.parametrize("plot_type", plot_types_SV)
@pytest.mark.parametrize("path", path_options)
def test_process_input_SV(plot_type, path):
    file_path = os.path.join(SPP_SV, path, f"example.SV{plot_type}.tsv")
    input_data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)

    if path == "wrong_rows":
        with pytest.raises(
            ValueError, match=f"Input matrix file should have {plot_type} rows"
        ):
            process_input(input_data, plot_type)
    else:
        ordered_input_data = process_input(input_data, plot_type)
        ordered_input_data_index = ordered_input_data.index.tolist()
        expected_index = get_context_reference(plot_type)
        assert ordered_input_data_index == expected_index

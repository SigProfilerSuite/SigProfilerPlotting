import importlib
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from PIL import Image

import sigProfilerPlotting as sigplt
from tests.image_comparison import (
    MINIMUM_CHROMA_SIMILARITY,
    MINIMUM_STRUCTURAL_SIMILARITY,
    image_similarities,
)


plotting = importlib.import_module("sigProfilerPlotting.sigProfilerPlotting")
TESTS_PATH = os.path.dirname(os.path.abspath(__file__))
REPOSITORY_ROOT = os.path.dirname(TESTS_PATH)
SPP_STANDARD_PNG = os.path.join(TESTS_PATH, "standard_png")
CNV48_PATH = os.path.join(TESTS_PATH, "input", "CNV", "unordered", "example.CNV48.tsv")
SBS_INPUT_PATH = os.path.join(TESTS_PATH, "input", "SBS", "ordered")


@pytest.mark.parametrize("as_dataframe", [False, True])
@pytest.mark.parametrize("percentage", [False, True])
def test_plot_cnv_aggregate_accepts_file_and_dataframe(as_dataframe, percentage):
    matrix = pd.read_csv(CNV48_PATH, sep="\t") if as_dataframe else CNV48_PATH

    images = sigplt.plotCNV(
        matrix,
        output_path="unused-for-pil-images",
        project="aggregate",
        percentage=percentage,
        aggregate=True,
        savefig_format="PIL_Image",
    )

    assert list(images) == [""]
    assert images[""].size[0] > 0
    images[""].close()


def test_plot_cnv_aggregate_requires_a_sample_column():
    matrix = pd.DataFrame({"MutationType": plotting.get_context_reference("48")})

    with pytest.raises(ValueError, match="at least one sample column"):
        sigplt.plotCNV(
            matrix,
            output_path="unused-for-pil-images",
            project="aggregate",
            aggregate=True,
            read_from_file=False,
            savefig_format="PIL_Image",
        )


def test_plot_cnv_aggregate_averages_across_numeric_sample_columns(monkeypatch):
    matrix = pd.DataFrame(
        {
            "MutationType": plotting.get_context_reference("48"),
            "sample-a": 2,
            "sample-b": 4,
        }
    )
    monkeypatch.setattr(
        plotting,
        "output_results",
        lambda savefig_format, output_path, project, figs, context_type, dpi: figs,
    )

    figures = sigplt.plotCNV(
        matrix,
        output_path="unused-for-captured-figures",
        project="aggregate",
        aggregate=True,
        read_from_file=False,
        savefig_format="PIL_Image",
    )

    bar_heights = [patch.get_height() for patch in figures[""].axes[0].patches[:48]]
    assert bar_heights == [3] * 48
    plt.close(figures[""])


def test_png_output_path_does_not_require_trailing_separator(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    plotting.output_results("png", str(tmp_path), "project", {"sample": fig}, "SBS_96")

    assert (tmp_path / "SBS_96_plots_sample.png").is_file()
    assert not (tmp_path.parent / f"{tmp_path.name}SBS_96_plots_sample.png").exists()


def test_process_input_rejects_wrong_channel_set_with_correct_row_count():
    contexts = plotting.get_context_reference("96")
    matrix = pd.DataFrame(
        {"MutationType": ["not-a-channel", *contexts[1:]], "sample": 1}
    )

    with pytest.raises(ValueError, match=r"missing:.*unexpected:"):
        plotting.process_input(matrix, "96")


def test_process_input_rejects_duplicate_channels():
    contexts = plotting.get_context_reference("96")
    matrix = pd.DataFrame(
        {"MutationType": [contexts[0], contexts[0], *contexts[2:]], "sample": 1}
    )

    with pytest.raises(ValueError, match="duplicate mutation contexts"):
        plotting.process_input(matrix, "96")


@pytest.mark.parametrize("plot_type", ["384", "1536", "4608"])
def test_legacy_sbs_matrix_is_canonicalized_before_plotting(plot_type, tmp_path):
    source = os.path.join(SBS_INPUT_PATH, f"example.SBS{plot_type}.all")
    matrix = pd.read_csv(source, sep="\t", index_col=0)
    shuffled_path = tmp_path / f"shuffled-SBS{plot_type}.tsv"
    matrix.sample(frac=1, random_state=17).to_csv(shuffled_path, sep="\t")

    with plotting._open_canonical_legacy_matrix(
        str(shuffled_path), plot_type
    ) as canonical_file:
        canonical = pd.read_csv(canonical_file, sep="\t", index_col=0)

    assert canonical.index.tolist() == plotting.get_context_reference(plot_type)
    pd.testing.assert_frame_equal(canonical, matrix)


@pytest.mark.parametrize("plot_type", ["384", "1536", "4608"])
def test_legacy_sbs_invalid_count_raises_value_error(plot_type, tmp_path):
    source = os.path.join(SBS_INPUT_PATH, f"example.SBS{plot_type}.all")
    matrix = pd.read_csv(source, sep="\t", index_col=0).astype(object)
    matrix.iloc[0, 0] = "not-a-count"
    invalid_path = tmp_path / f"invalid-SBS{plot_type}.tsv"
    matrix.to_csv(invalid_path, sep="\t")

    with pytest.raises(ValueError, match="Invalid matrix value.*not-a-count"):
        sigplt.plotSBS(
            str(invalid_path),
            str(tmp_path),
            "invalid",
            plot_type,
        )


@pytest.mark.parametrize(
    ("function", "plot_type"),
    [
        (sigplt.plotSBS, "96"),
        (sigplt.plotSBS, "288"),
        (sigplt.plotID, "83"),
        (sigplt.plotDBS, "78"),
    ],
)
def test_zero_event_percentage_plot_returns_an_image(function, plot_type):
    matrix = pd.DataFrame(
        {"zero-sample": 0}, index=plotting.get_context_reference(plot_type)
    )
    matrix.index.name = plotting.MUTTYPE

    images = function(
        matrix,
        "unused-for-pil-images",
        "zero-events",
        plot_type,
        percentage=True,
        savefig_format="PIL_Image",
    )

    assert list(images) == ["zero-sample"]
    assert images["zero-sample"].size[0] > 0
    images["zero-sample"].close()


def test_default_template_is_cached_outside_the_package(monkeypatch, tmp_path):
    monkeypatch.delenv("SIGPROFILERPLOTTING_VOLUME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg-cache"))

    cache_dir = plotting._default_template_cache_dir()
    template = plotting.make_pickle_file("SBS96", return_plot_template=True)

    assert os.path.commonpath([cache_dir, str(tmp_path)]) == str(tmp_path)
    assert f"spp-{plotting.SPP_VERSION}" in cache_dir
    assert f"matplotlib-{plotting.matplotlib.__version__}" in cache_dir
    assert os.path.isfile(os.path.join(cache_dir, "SBS96.pkl"))
    assert template is not None
    plt.close(template)


def test_explicit_template_volume_remains_supported(tmp_path):
    cache_dir = tmp_path / "explicit-cache"

    template = plotting.make_pickle_file(
        "SBS96", return_plot_template=True, volume=str(cache_dir)
    )

    assert (cache_dir / "SBS96.pkl").is_file()
    plt.close(template)


def test_explicit_unwritable_template_volume_remains_strict(tmp_path):
    cache_blocker = tmp_path / "not-a-directory"
    cache_blocker.write_text("block cache directory creation")

    with pytest.raises(NotADirectoryError):
        plotting.make_pickle_file(
            "SBS96",
            return_plot_template=True,
            volume=str(cache_blocker / "cache"),
        )


def test_implicit_unwritable_template_cache_falls_back_to_memory(monkeypatch, tmp_path):
    cache_blocker = tmp_path / "not-a-directory"
    cache_blocker.write_text("block cache directory creation")
    monkeypatch.delenv("SIGPROFILERPLOTTING_VOLUME", raising=False)
    monkeypatch.setattr(
        plotting, "_default_template_cache_dir", lambda: str(cache_blocker / "cache")
    )

    template = plotting.make_pickle_file("SBS96", return_plot_template=True)

    assert template is not None
    plt.close(template)


def test_corrupt_implicit_template_cache_is_regenerated(monkeypatch, tmp_path):
    cache_dir = tmp_path / "implicit-cache"
    cache_dir.mkdir()
    cache_file = cache_dir / "SBS96.pkl"
    cache_file.write_bytes(b"not a pickle")
    monkeypatch.delenv("SIGPROFILERPLOTTING_VOLUME", raising=False)
    monkeypatch.setattr(
        plotting, "_default_template_cache_dir", lambda: str(cache_dir)
    )

    template = plotting.make_pickle_file("SBS96", return_plot_template=True)

    assert template is not None
    with cache_file.open("rb") as template_file:
        regenerated = plotting.pickle.load(template_file)
    assert regenerated is not None
    plt.close(regenerated)
    plt.close(template)


def test_wrong_type_in_implicit_template_cache_is_regenerated(monkeypatch, tmp_path):
    cache_dir = tmp_path / "implicit-cache"
    cache_dir.mkdir()
    cache_file = cache_dir / "SBS96.pkl"
    with cache_file.open("wb") as template_file:
        plotting.pickle.dump(None, template_file)
    monkeypatch.delenv("SIGPROFILERPLOTTING_VOLUME", raising=False)
    monkeypatch.setattr(
        plotting, "_default_template_cache_dir", lambda: str(cache_dir)
    )

    template = plotting.make_pickle_file("SBS96", return_plot_template=True)

    assert isinstance(template, plotting.Figure)
    plt.close(template)


def test_wrong_type_in_explicit_template_cache_remains_strict(tmp_path):
    cache_dir = tmp_path / "explicit-cache"
    cache_dir.mkdir()
    with (cache_dir / "SBS96.pkl").open("wb") as template_file:
        plotting.pickle.dump(None, template_file)

    with pytest.raises(TypeError, match="matplotlib.figure.Figure"):
        plotting.make_pickle_file(
            "SBS96", return_plot_template=True, volume=str(cache_dir)
        )


def test_empty_template_volume_environment_uses_implicit_fallback(
    monkeypatch, tmp_path
):
    cache_blocker = tmp_path / "not-a-directory"
    cache_blocker.write_text("block cache directory creation")
    monkeypatch.setenv("SIGPROFILERPLOTTING_VOLUME", "")
    monkeypatch.setattr(
        plotting, "_default_template_cache_dir", lambda: str(cache_blocker / "cache")
    )

    template = plotting.make_pickle_file("SBS96", return_plot_template=True)

    assert isinstance(template, plotting.Figure)
    plt.close(template)


@pytest.mark.parametrize(
    ("function", "plot_type"),
    [(sigplt.plotSBS, "384"), (sigplt.plotID, "28"), (sigplt.plotDBS, "186")],
)
def test_legacy_context_rejects_non_pdf_format(function, plot_type, tmp_path):
    with pytest.raises(
        ValueError, match=rf"context '{plot_type}'.*only savefig_format='pdf'"
    ):
        function(
            "matrix-is-not-read-before-format-validation",
            str(tmp_path),
            "project",
            plot_type,
            savefig_format="png",
        )


def test_import_does_not_globally_suppress_warnings():
    code = (
        "import warnings; "
        "warnings.resetwarnings(); "
        "import sigProfilerPlotting; "
        "warnings.warn('sigprofiler-warning-marker', UserWarning)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": REPOSITORY_ROOT},
    )

    assert "sigprofiler-warning-marker" in result.stderr


def test_portable_image_comparison_rejects_missing_content(tmp_path):
    reference_path = os.path.join(SPP_STANDARD_PNG, "SBS_96_plots_bars.png")
    altered_path = tmp_path / "missing-content.png"
    with Image.open(reference_path) as reference:
        altered = np.asarray(reference.convert("RGB")).copy()
    altered[altered.shape[0] // 3 : 2 * altered.shape[0] // 3, :, :] = 255
    Image.fromarray(altered).save(altered_path)

    structure, _ = image_similarities(altered_path, reference_path)
    assert structure < MINIMUM_STRUCTURAL_SIMILARITY


def test_portable_image_comparison_accepts_small_canvas_change(tmp_path):
    reference_path = os.path.join(SPP_STANDARD_PNG, "CNV_48_plots_Random.png")
    altered_path = tmp_path / "resized.png"
    with Image.open(reference_path) as reference:
        reference.resize(
            (round(reference.width * 1.03), round(reference.height * 1.03)),
            Image.Resampling.LANCZOS,
        ).save(altered_path)

    structure, chroma = image_similarities(altered_path, reference_path)
    assert structure >= MINIMUM_STRUCTURAL_SIMILARITY
    assert chroma >= MINIMUM_CHROMA_SIMILARITY


def test_portable_image_comparison_rejects_large_canvas_change(tmp_path):
    reference_path = os.path.join(SPP_STANDARD_PNG, "CNV_48_plots_Random.png")
    altered_path = tmp_path / "stretched.png"
    with Image.open(reference_path) as reference:
        reference.resize(
            (round(reference.width * 1.25), reference.height),
            Image.Resampling.LANCZOS,
        ).save(altered_path)

    assert image_similarities(altered_path, reference_path) == (0.0, 0.0)


@pytest.mark.parametrize(
    "reference_name", ["CNV_48_plots_bars.png", "SV_32_plots_bars.png"]
)
def test_portable_image_comparison_rejects_missing_bars(tmp_path, reference_name):
    reference_path = os.path.join(SPP_STANDARD_PNG, reference_name)
    altered_path = tmp_path / "missing-bars.png"
    with Image.open(reference_path) as reference:
        altered = np.asarray(reference.convert("RGB")).copy()
    altered[2 * altered.shape[0] // 3 :, :, :] = 255
    Image.fromarray(altered).save(altered_path)

    structure, chroma = image_similarities(altered_path, reference_path)
    assert (
        structure < MINIMUM_STRUCTURAL_SIMILARITY
        or chroma < MINIMUM_CHROMA_SIMILARITY
    )


@pytest.mark.parametrize(
    "reference_name",
    [
        "SBS_96_plots_bars.png",
        "CNV_48_plots_Random.png",
        "CNV_48_plots_bars.png",
        "SV_32_plots_Random.png",
        "SV_32_plots_bars.png",
    ],
)
def test_portable_image_comparison_rejects_color_loss(tmp_path, reference_name):
    reference_path = os.path.join(SPP_STANDARD_PNG, reference_name)
    altered_path = tmp_path / "grayscale.png"
    with Image.open(reference_path) as reference:
        reference.convert("L").convert("RGB").save(altered_path)

    _, chroma = image_similarities(altered_path, reference_path)
    assert chroma < MINIMUM_CHROMA_SIMILARITY

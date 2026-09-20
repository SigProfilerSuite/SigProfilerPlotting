#!/usr/bin/env python3

# Author: Erik Bergstrom

# Contact: ebergstr@eng.ucsd.edu

import argparse
import copy
import errno
import io
import itertools
import logging
import os
import pickle
import re
import string
import sys
from bdb import set_trace
from collections import OrderedDict

import matplotlib
import matplotlib.font_manager
import matplotlib.lines as lines
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import sklearn
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.ticker import LinearLocator
from PIL import Image
from sklearn.preprocessing import LabelEncoder

import sigProfilerPlotting as spplt
from .version import short_version as SPP_VERSION

matplotlib.use("Agg")

MUTTYPE = "MutationType"
INDEX_VALS = ["MutationType", "index", "Mutation Types", "classification"]
SPP_PATH = spplt.__path__[0]
SPP_TEMPLATES = os.path.join(SPP_PATH, "templates/")
SPP_FONTS = os.path.join(SPP_PATH, "fonts/")
SPP_REFERENCE = os.path.join(SPP_PATH, "reference_formats/")

_FONTS_LOADED = False

logging.getLogger("matplotlib.font_manager").disabled = False

type_dict = {
    "96": "SBS96.txt",
    "sbs": "SBS96.txt",
    "sbs96": "SBS96.txt",
    "288": "SBS288.txt",
    "sbs288": "SBS288.txt",
    "384": "SBS384.txt",
    "sbs384": "SBS384.txt",
    "sbs1536": "SBS1536.txt",
    "1536": "SBS1536.txt",
    "4608": "SBS4608.txt",
    "sbs4608": "SBS4608.txt",
    "sbs6144": "SBS6144.txt",
    "6144": "SBS6144.txt",
    "78": "DBS78.txt",
    "dbs": "DBS78.txt",
    "dbs78": "DBS78.txt",
    "dinuc": "DBS78.txt",
    "83": "ID83.txt",
    "id": "ID83.txt",
    "id83": "ID83.txt",
    "cnv48": "CNV48.txt",
    "48": "CNV48.txt",
    "sv32": "SV32.txt",
    "32": "SV32.txt",
}


# Loads fonts required for plotting
def load_custom_fonts():
    global _FONTS_LOADED
    if not _FONTS_LOADED:
        for font_file in os.listdir(SPP_FONTS):
            if font_file.endswith(".ttf"):
                try:
                    font_path = os.path.join(SPP_FONTS, font_file)
                    matplotlib.font_manager.fontManager.addfont(font_path)
                except:
                    print("ERROR loading font: " + font_file)
    _FONTS_LOADED = True


# Note that plt.close(), plt.clf(), and plt.cla() would not close memory
# Referenced the following post for the function below:
# https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
def clear_plotting_memory():
    usedbackend = matplotlib.get_backend()
    matplotlib.use(usedbackend)
    allfignums = matplotlib.pyplot.get_fignums()
    for i in allfignums:
        fig = matplotlib.pyplot.figure(i)
        fig.clear()
        matplotlib.pyplot.close(fig)


def _validate_savefig_format(savefig_format):
    if not isinstance(savefig_format, str) or savefig_format.lower() not in {
        "pdf",
        "png",
        "pil_image",
    }:
        raise ValueError("ERROR: savefig_format must be 'pdf', 'png', or 'PIL_Image'.")


def _require_pdf_for_legacy_context(
    savefig_format, plot_type, modern_contexts, function_name
):
    _validate_savefig_format(savefig_format)
    if plot_type not in modern_contexts and savefig_format.lower() != "pdf":
        raise ValueError(
            f"{function_name} context '{plot_type}' supports only savefig_format='pdf'."
        )


def _parse_mutation_value(value, percentage, sample=None, mutation_type=None):
    """Parse a matrix value without allowing a failed conversion to be ignored."""
    try:
        parsed = float(value) if percentage else int(value)
    except (TypeError, ValueError) as error:
        location = []
        if sample is not None:
            location.append(f"sample '{sample}'")
        if mutation_type is not None:
            location.append(f"mutation type '{mutation_type}'")
        location_text = " for " + ", ".join(location) if location else ""
        expected = "a numeric value" if percentage else "an integer mutation count"
        raise ValueError(
            f"Invalid matrix value {value!r}{location_text}; expected {expected}."
        ) from error

    if not np.isfinite(parsed):
        raise ValueError(f"Invalid non-finite matrix value {value!r}.")
    return parsed


def _open_canonical_legacy_matrix(matrix_path, plot_type):
    """Open a standard legacy SBS matrix in canonical reference order.

    PCAWG comma-separated inputs have a different multi-column context layout and
    retain their existing parser. Standard tab-separated matrices are validated
    and reindexed through ``process_input`` before their legacy plotter reads them.
    """
    with open(matrix_path, encoding="utf-8") as matrix_file:
        header = matrix_file.readline()
    if "," in header and "\t" not in header:
        return open(matrix_path, encoding="utf-8")

    data = process_input(matrix_path, plot_type)
    canonical_matrix = io.StringIO()
    data.to_csv(canonical_matrix, sep="\t", index_label=MUTTYPE)
    canonical_matrix.seek(0)
    return canonical_matrix


# Saves figures to files, unless savefig_format is "PIL_Image", in which case
# the figures are saved to a dictionary of buffers
def output_results(savefig_format, output_path, project, figs, context_type, dpi=100):
    _validate_savefig_format(savefig_format)
    if savefig_format.lower() == "pdf":
        file_path = os.path.join(output_path, f"{context_type}_plots_{project}.pdf")
        pp = PdfPages(file_path)
        for fig in figs:
            if context_type in ("CNV_48", "SV_32"):
                figs[fig].savefig(pp, format="pdf", bbox_inches="tight")
            else:
                figs[fig].savefig(pp, format="pdf")
        pp.close()
        clear_plotting_memory()
    elif savefig_format.lower() == "png":
        for fig in figs:
            file_path = os.path.join(
                output_path, f"{context_type}_plots_{fig}.png"
            )
            if context_type in ("CNV_48", "SV_32"):
                figs[fig].savefig(
                    file_path,
                    dpi=dpi,
                    bbox_inches="tight",
                )
            else:
                figs[fig].savefig(file_path, dpi=dpi)
        clear_plotting_memory()
    elif savefig_format.lower() == "pil_image":
        image_list = {}
        for fig in figs:
            tmp_buffer = io.BytesIO()
            if context_type in ("CNV_48", "SV_32"):
                figs[fig].savefig(
                    tmp_buffer, format="png", bbox_inches="tight", dpi=dpi
                )
            else:
                figs[fig].savefig(tmp_buffer, format="png", dpi=dpi)
            # convert tmp_buffer to a PIL and close buffer
            tmp_buffer.seek(0)
            tmp_image = Image.open(tmp_buffer)
            # add the image to the image list for return
            image_list[fig] = tmp_image
        clear_plotting_memory()
        return image_list
    return None


# Get corresponding reference index from our reference_format folder
def get_context_reference(plot_type):
    ref_index = []
    if plot_type.lower() in type_dict:
        SPP_TYPE = type_dict[plot_type.lower()]
    else:
        raise ValueError(
            "ERROR: SigProfilerPlotting is currently not supporting this input plot_type."
        )

    ref_index = pd.read_csv(SPP_REFERENCE + SPP_TYPE, sep="\t", header=None)
    ref_index = ref_index.iloc[:, 0].tolist()

    return ref_index


def process_input(matrix_path, plot_type):
    # input data is a DataFrame
    if isinstance(matrix_path, pd.DataFrame):
        # copy dataframe with deepcopy
        data = matrix_path.copy()
        # index is a non-standard value
        if data.index.name not in INDEX_VALS:
            if MUTTYPE in data.columns:
                data = data.set_index(MUTTYPE, drop=True)
            # the first column is non-MUTTYPE and non-integer
            elif not data.iloc[:, 0].apply(lambda x: isinstance(x, int)).all():
                data.rename(columns={data.columns[0]: MUTTYPE}, inplace=True)
                data = data.set_index(data.columns[0], drop=True)
            else:
                data = data.reset_index()
                data.rename(columns={data.columns[0]: MUTTYPE}, inplace=True)
                data = data.set_index(MUTTYPE, drop=True)
        else:
            # Note: set the index to MUTTYPE for consistency with the rest of the code
            data.index.name = MUTTYPE
    # input data is a file path
    elif isinstance(matrix_path, str):
        data = pd.read_csv(matrix_path, sep="\t", index_col=0)
        data = data.dropna(axis=1, how="all")
        data.index.name = MUTTYPE
    # input data is a numpy array
    elif isinstance(matrix_path, np.ndarray):
        # Note: ndarray does not have index or column names and is not recommended
        data = pd.DataFrame(matrix_path)
        # add index of mutation type to the dataframe
        if plot_type.lower() in type_dict:
            data.index = get_context_reference(plot_type)
    else:
        raise ValueError(
            "ERROR: matrix_path requires pd.DataFrame, path to file, or np.ndarray, not "
            + f"{type(matrix_path)}."
        )

    if data.isnull().values.any():
        raise ValueError("ERROR: matrix_path contains Nans.")

    def order_input_context(plot_type, input_data):
        if plot_type.lower() in type_dict:
            ref_format = get_context_reference(plot_type)
            if input_data.shape[0] != len(ref_format):
                raise ValueError(
                    "Input matrix file should have "
                    + str(len(ref_format))
                    + " rows"
                )

            duplicate_contexts = input_data.index[input_data.index.duplicated()].unique()
            if len(duplicate_contexts):
                raise ValueError(
                    "Input matrix contains duplicate mutation contexts: "
                    + ", ".join(map(str, duplicate_contexts))
                )

            actual_contexts = set(input_data.index)
            expected_contexts = set(ref_format)
            missing_contexts = expected_contexts - actual_contexts
            unexpected_contexts = actual_contexts - expected_contexts
            if missing_contexts or unexpected_contexts:
                details = []
                if missing_contexts:
                    details.append(
                        "missing: " + ", ".join(sorted(map(str, missing_contexts)))
                    )
                if unexpected_contexts:
                    details.append(
                        "unexpected: "
                        + ", ".join(sorted(map(str, unexpected_contexts)))
                    )
                raise ValueError(
                    "Input matrix mutation contexts do not match the reference ("
                    + "; ".join(details)
                    + ")."
                )

            reindexed_data = input_data.reindex(ref_format)
        else:
            # If a non-standard context is used, no sort is applied
            reindexed_data = input_data
        return reindexed_data

    return order_input_context(plot_type, data)


def get_default_96labels():
    first = ["A", "C", "G", "T"]
    inner_bracket = [[x] * 16 for x in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = [x for x in list(itertools.product(first, first))]
    result = [
        outter_bracket[f % 16][0]
        + "["
        + inner_bracket[f]
        + "]"
        + outter_bracket[f % 16][1]
        for f in range(0, 96)
    ]
    return result


def _default_template_cache_dir():
    """Return a writable, per-user cache location for generated templates."""
    if sys.platform == "darwin":
        cache_root = os.path.expanduser("~/Library/Caches")
    elif os.name == "nt":
        cache_root = os.getenv(
            "LOCALAPPDATA", os.path.expanduser(r"~\AppData\Local")
        )
    else:
        cache_root = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
    compatibility_key = (
        f"spp-{SPP_VERSION}_matplotlib-{matplotlib.__version__}_"
        f"py-{sys.version_info.major}.{sys.version_info.minor}"
    )
    return os.path.join(
        cache_root, "SigProfilerPlotting", "templates", compatibility_key
    )


def _finish_plot_template(plot, path, return_plot_template, strict_cache):
    if path is not None:
        try:
            with open(path, "wb") as template_file:
                pickle.dump(plot, template_file)
        except OSError:
            if strict_cache:
                raise
    if return_plot_template:
        return plot


def make_pickle_file(context="SBS96", return_plot_template=False, volume=None):

    # The environmental variable takes precedence over the volume argument
    # If the environmental variable is not set, the volume argument is used
    configured_volume = os.getenv("SIGPROFILERPLOTTING_VOLUME") or None
    requested_volume = volume or None
    strict_cache = configured_volume is not None or requested_volume is not None
    volume = configured_volume or requested_volume

    # Use the default volume when no environmental variable or volume argument is provided
    if volume is None:
        volume = _default_template_cache_dir()

    path = os.path.join(volume, context + ".pkl")

    # if the pickle file already exists, return the template
    if os.path.exists(path):
        try:
            with open(path, "rb") as template_file:
                cached_template = pickle.load(template_file)
            if not isinstance(cached_template, Figure):
                raise TypeError(
                    "Cached plot template must be a matplotlib.figure.Figure, not "
                    f"{type(cached_template).__name__}."
                )
            return cached_template
        except (
            OSError,
            EOFError,
            pickle.PickleError,
            AttributeError,
            ImportError,
            ValueError,
            TypeError,
        ):
            if strict_cache:
                raise
            # An implicit cache is disposable. Regenerate stale or incomplete
            # templates instead of making plotting fail permanently.
            try:
                os.remove(path)
            except OSError:
                pass

    # check if the template directory exists, create if not
    try:
        os.makedirs(volume, exist_ok=True)
    except OSError:
        if strict_cache:
            raise
        path = None

    if context == "SBS96":
        plot_custom_text = False
        sig_probs = False
        pcawg = False

        # total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
        # , extent=[-5, 80, -5, 30])
        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="lightgray")
        panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
        seq96 = [
            "A[C>A]A",
            "A[C>A]C",
            "A[C>A]G",
            "A[C>A]T",
            "C[C>A]A",
            "C[C>A]C",
            "C[C>A]G",
            "C[C>A]T",
            "G[C>A]A",
            "G[C>A]C",
            "G[C>A]G",
            "G[C>A]T",
            "T[C>A]A",
            "T[C>A]C",
            "T[C>A]G",
            "T[C>A]T",
            "A[C>G]A",
            "A[C>G]C",
            "A[C>G]G",
            "A[C>G]T",
            "C[C>G]A",
            "C[C>G]C",
            "C[C>G]G",
            "C[C>G]T",
            "G[C>G]A",
            "G[C>G]C",
            "G[C>G]G",
            "G[C>G]T",
            "T[C>G]A",
            "T[C>G]C",
            "T[C>G]G",
            "T[C>G]T",
            "A[C>T]A",
            "A[C>T]C",
            "A[C>T]G",
            "A[C>T]T",
            "C[C>T]A",
            "C[C>T]C",
            "C[C>T]G",
            "C[C>T]T",
            "G[C>T]A",
            "G[C>T]C",
            "G[C>T]G",
            "G[C>T]T",
            "T[C>T]A",
            "T[C>T]C",
            "T[C>T]G",
            "T[C>T]T",
            "A[T>A]A",
            "A[T>A]C",
            "A[T>A]G",
            "A[T>A]T",
            "C[T>A]A",
            "C[T>A]C",
            "C[T>A]G",
            "C[T>A]T",
            "G[T>A]A",
            "G[T>A]C",
            "G[T>A]G",
            "G[T>A]T",
            "T[T>A]A",
            "T[T>A]C",
            "T[T>A]G",
            "T[T>A]T",
            "A[T>C]A",
            "A[T>C]C",
            "A[T>C]G",
            "A[T>C]T",
            "C[T>C]A",
            "C[T>C]C",
            "C[T>C]G",
            "C[T>C]T",
            "G[T>C]A",
            "G[T>C]C",
            "G[T>C]G",
            "G[T>C]T",
            "T[T>C]A",
            "T[T>C]C",
            "T[T>C]G",
            "T[T>C]T",
            "A[T>G]A",
            "A[T>G]C",
            "A[T>G]G",
            "A[T>G]T",
            "C[T>G]A",
            "C[T>G]C",
            "C[T>G]G",
            "C[T>G]T",
            "G[T>G]A",
            "G[T>G]C",
            "G[T>G]G",
            "G[T>G]T",
            "T[T>G]A",
            "T[T>G]C",
            "T[T>G]G",
            "T[T>G]T",
        ]
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        xlabels = [seq[0] + seq[2] + seq[6] for seq in seq96]
        i = 0

        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 6, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y3),
                    0.15,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.159

        yText = y3 + 0.06
        plt.text(
            0.1,
            yText,
            "C>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "C>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.415,
            yText,
            "C>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.575,
            yText,
            "T>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.735,
            yText,
            "T>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.89,
            yText,
            "T>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        # ytick_offest = int(y/4)
        y = ymax / 1.025
        ytick_offest = float(y / 3)
        labs = np.arange(0.375, 96.375, 1)

        panel1.set_xlim([0, 96])
        # panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        # panel1.set_yticks(ylabs)
        count = 0
        m = 0
        for i in range(0, 96, 1):
            plt.text(
                i / 101 + 0.0415,
                0.02,
                xlabels[i][0],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.044,
                xlabels[i][1],
                fontsize=30,
                color=colors[m],
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                fontweight="bold",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.071,
                xlabels[i][2],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            count += 1
            if count == 16:
                count = 0
                m += 1

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]
        return _finish_plot_template(
            plot1, path, return_plot_template, strict_cache
        )
    elif context == "SBS288":
        plot_custom_text = False
        sig_probs = False
        pcawg = False

        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="lightgray")
        panel1 = plt.axes([0.04, 0.09, 0.7, 0.77])
        panel2 = plt.axes([0.77, 0.09, 0.21, 0.77])
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        i = 0
        result = get_default_96labels()
        xlabels = [seq[0] + seq[2] + seq[6] for seq in result]
        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 6, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y3),
                    0.11,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.117

        yText = y3 + 0.06

        plt.text(
            0.082,
            yText,
            "C>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.1975,
            yText,
            "C>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.315,
            yText,
            "C>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.43,
            yText,
            "T>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.55,
            yText,
            "T>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.665,
            yText,
            "T>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        y = ymax / 1.025
        ytick_offest = float(y / 3)
        font_label_size = 30
        labs = np.arange(0.375, 96.375, 1)

        panel1.set_xlim([0, 96])
        panel1.set_ylim([0, y])
        panel1.set_xticks(labs)

        count = 0
        m = 0
        for i in range(0, 96, 1):
            plt.text(
                i / 137 + 0.04,
                0.02,
                xlabels[i][0],
                fontsize=25,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 137 + 0.04,
                0.044,
                xlabels[i][1],
                fontsize=25,
                color=colors[m],
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                fontweight="bold",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 137 + 0.04,
                0.071,
                xlabels[i][2],
                fontsize=25,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            count += 1
            if count == 16:
                count = 0
                m += 1

        panel1.yaxis.grid(True)
        panel1.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in panel1.get_yticklabels()]

        yp2 = 28
        labels = []
        y2max = 0
        tsbColors = [
            [1 / 256, 70 / 256, 102 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            "green",
        ]

        y = int(y2max * 1.1)
        if y <= 4:
            y += 4
        while y % 4 != 0:
            y += 1
        ytick_offest = int(y / 4)

        panel2.spines["right"].set_visible(False)
        panel2.spines["top"].set_visible(False)
        labels.reverse()
        panel2.set_yticks([3, 7, 11, 15, 19, 23, 27])
        panel2.set_yticklabels(labels, fontsize=30, fontname="Arial", weight="bold")
        panel2.set_xticklabels(xlabels, fontsize=30)
        handles, labels = panel2.get_legend_handles_labels()
        panel2.legend(handles[:3], labels[:3], loc="best", prop={"size": 30})
        return _finish_plot_template(
            plot1, path, return_plot_template, strict_cache
        )
    elif context == "DBS78":
        plot_custom_text = False
        pcawg = False
        sig_probs = False
        plt.rcParams["axes.linewidth"] = 4
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="grey")
        panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [3 / 256, 102 / 256, 204 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [1 / 256, 102 / 256, 1 / 256],
            [255 / 256, 153 / 256, 153 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [255 / 256, 178 / 256, 102 / 256],
            [255 / 256, 128 / 256, 1 / 256],
            [204 / 256, 153 / 256, 255 / 256],
            [76 / 256, 1 / 256, 153 / 256],
        ]

        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        i = 0
        panel1.add_patch(
            plt.Rectangle(
                (0.043, y3),
                0.101,
                0.05,
                facecolor=colors[0],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.151, y3),
                0.067,
                0.05,
                facecolor=colors[1],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.225, y3),
                0.102,
                0.05,
                facecolor=colors[2],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.334, y3),
                0.067,
                0.05,
                facecolor=colors[3],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.408, y3),
                0.102,
                0.05,
                facecolor=colors[4],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.517, y3),
                0.067,
                0.05,
                facecolor=colors[5],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.591, y3),
                0.067,
                0.05,
                facecolor=colors[6],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.665, y3),
                0.102,
                0.05,
                facecolor=colors[7],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.774, y3),
                0.102,
                0.05,
                facecolor=colors[8],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.883, y3),
                0.102,
                0.05,
                facecolor=colors[9],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )

        yText = y3 + 0.06
        plt.text(
            0.07,
            yText,
            "AC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.163,
            yText,
            "AT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "CC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.345,
            yText,
            "CG>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.435,
            yText,
            "CT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.527,
            yText,
            "GC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.6,
            yText,
            "TA>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.69,
            yText,
            "TC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.8,
            yText,
            "TG>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.915,
            yText,
            "TT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        ytick_offest = int(y / 4)

        labs = np.arange(0.44, 78.44, 1)
        panel1.set_xlim([0, 78])
        panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        panel1.set_xticklabels(
            xlabels,
            rotation="vertical",
            fontsize=30,
            color="grey",
            fontname="Courier New",
            verticalalignment="top",
            fontweight="bold",
        )

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")
        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=True,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]
        [i.set_color("grey") for i in plt.gca().get_xticklabels()]
        return _finish_plot_template(
            plot1, path, return_plot_template, strict_cache
        )
    elif context == "ID83":
        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 12))
        plt.rc("axes", edgecolor="black")
        panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [253 / 256, 190 / 256, 111 / 256],
            [255 / 256, 128 / 256, 2 / 256],
            [176 / 256, 221 / 256, 139 / 256],
            [54 / 256, 161 / 256, 46 / 256],
            [253 / 256, 202 / 256, 181 / 256],
            [252 / 256, 138 / 256, 106 / 256],
            [241 / 256, 68 / 256, 50 / 256],
            [188 / 256, 25 / 256, 26 / 256],
            [208 / 256, 225 / 256, 242 / 256],
            [148 / 256, 196 / 256, 223 / 256],
            [74 / 256, 152 / 256, 201 / 256],
            [23 / 256, 100 / 256, 171 / 256],
            [226 / 256, 226 / 256, 239 / 256],
            [182 / 256, 182 / 256, 216 / 256],
            [134 / 256, 131 / 256, 189 / 256],
            [98 / 256, 64 / 256, 155 / 256],
        ]

        x = 0.0475
        y_top = 0.827
        y_bottom = 0.114
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 12, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y_top),
                    0.0595,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            panel1.add_patch(
                plt.Rectangle(
                    (x, y_bottom),
                    0.0595,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.0665

        panel1.add_patch(
            plt.Rectangle(
                (x - 0.001, y_top),
                0.006,
                0.05,
                facecolor=colors[12],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x - 0.001, y_bottom),
                0.006,
                0.05,
                facecolor=colors[12],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.011
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.0155,
                0.05,
                facecolor=colors[13],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.0155,
                0.05,
                facecolor=colors[13],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.022
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.027,
                0.05,
                facecolor=colors[14],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.027,
                0.05,
                facecolor=colors[14],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.0335
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.049,
                0.05,
                facecolor=colors[15],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.049,
                0.05,
                facecolor=colors[15],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )

        yText = y_top + 0.01
        plt.text(
            0.072,
            yText,
            "C",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.1385,
            yText,
            "T",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.205,
            yText,
            "C",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.2715,
            yText,
            "T",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.338,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.4045,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.471,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.5375,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.604,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.6705,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.737,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.8035,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.844,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.861,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.888,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.93,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )

        yText_labels_top = yText + 0.075
        yText_labels_bottom = y_bottom - 0.03
        yText_labels_bottom_sec = yText_labels_bottom - 0.045

        plt.text(
            0.08,
            yText_labels_top,
            "1bp Deletion",
            fontsize=40,
            fontname="Times New Roman",
            weight="bold",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.21,
            yText_labels_top,
            "1bp Insertion",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.375,
            yText_labels_top,
            ">1bp Deletion at Repeats\n      (Deletion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.64,
            yText_labels_top,
            ">1bp Insertion at Repeats\n       (Insertion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.85,
            yText_labels_top,
            " Microhomology\n(Deletion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        plt.text(
            0.058,
            yText_labels_bottom_sec,
            "Homopolymer Length",
            fontsize=35,
            fontname="Times New Roman",
            weight="bold",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.19,
            yText_labels_bottom_sec,
            "Homopolymer Length",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.39,
            yText_labels_bottom_sec,
            "Number of Repeat Units",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.65,
            yText_labels_bottom_sec,
            "Number of Repeat Units",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.85,
            yText_labels_bottom_sec,
            "Microhomology Length",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        x = 0.0477
        for i in range(0, 8, 1):
            if i != 2 and i != 3:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

            x += 0.0665

        for i in range(0, 4, 1):
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=32,
                fontweight="bold",
                fontname="Times New Roman",
                color="black",
                transform=plt.gcf().transFigure,
            )
            x += 0.0665

        plt.text(
            x,
            yText_labels_bottom,
            "1",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.011
        plt.text(
            x,
            yText_labels_bottom,
            "1  2",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.022
        plt.text(
            x,
            yText_labels_bottom,
            "1  2  3",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.0335
        plt.text(
            x,
            yText_labels_bottom,
            "1  2  3  4  5+",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        labs = np.arange(0.375, 83.375, 1)
        panel1.set_xlim([0, 83])
        panel1.set_ylim([0, y])
        panel1.set_xticks(labs)

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=False,
            labelleft=True,
            right=False,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="gray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]

        return _finish_plot_template(
            plot1, path, return_plot_template, strict_cache
        )


def getylabels(ylabels):
    if max(ylabels) >= 10**9:
        ylabels = ["{:.2e}".format(x) for x in ylabels]
        ylabels[0] = "0.00"
    else:
        if max(ylabels) <= 1000:
            ylabels = ["{:,.0f}".format(x) for x in ylabels]
            ylabels[0] = "0"
        elif max(ylabels) < 10**5 and max(ylabels) > 1000:
            ylabels = ["{:,.0f}".format(x / 1000) + "k" for x in ylabels]
            ylabels[0] = "0"
        else:  # if max(ylabels)>= 10**5:
            ylabels = ["{:,.0f}".format(x / (10**6)) + "m" for x in ylabels]
            ylabels[0] = "0"
    return ylabels


def getxlabels(xlabels):
    if max(xlabels) >= 10**10:
        xlabels = ["{:.2e}".format(x) for x in xlabels]
        xlabels[0] = "0.00"
    else:
        if max(xlabels) <= 1000:
            xlabels = ["{:,.1f}".format(x) for x in xlabels]
            xlabels[0] = "0"
        elif max(xlabels) < 10**6 and max(xlabels) > 1000:
            xlabels = ["{:,.2f}".format(x / 1000) + "k" for x in xlabels]
            xlabels[0] = "0.00"
        else:  # max(xlabels)>= 10**6:
            xlabels = ["{:,.2f}".format(x / (10**6)) + "m" for x in xlabels]
            xlabels[0] = "0.00"
    return xlabels


def reindex_sbs96(data_f):
    first = ["A", "C", "G", "T"]
    inner_bracket = [[x] * 16 for x in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = [x for x in list(itertools.product(first, first))]
    result = [
        outter_bracket[f % 16][0]
        + "["
        + inner_bracket[f]
        + "]"
        + outter_bracket[f % 16][1]
        for f in range(0, 96)
    ]
    data_f = data_f.reindex(result)
    return data_f


def reindex_sbs288(data_f):
    result = get_default_96labels()
    mutations_df = pd.DataFrame(index=result, columns=data_f.columns)
    T_mutations_df = pd.DataFrame(index=result, columns=data_f.columns)
    U_mutations_df = pd.DataFrame(index=result, columns=data_f.columns)
    N_mutations_df = pd.DataFrame(index=result, columns=data_f.columns)

    for row in result:
        T_mutations_df.loc[row] = data_f.loc["T:" + row]
        U_mutations_df.loc[row] = data_f.loc["U:" + row]
        N_mutations_df.loc[row] = data_f.loc["N:" + row]
        mutations_df.loc[row] = (
            data_f.loc["T:" + row] + data_f.loc["U:" + row] + data_f.loc["N:" + row]
        )

    mutations_TSB_df_T = pd.DataFrame(
        T_mutations_df.values.reshape((-1, 16, len(data_f.columns))).sum(axis=1),
        index=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
        columns=data_f.columns,
    )
    mutations_TSB_df_U = pd.DataFrame(
        U_mutations_df.values.reshape((-1, 16, len(data_f.columns))).sum(axis=1),
        index=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
        columns=data_f.columns,
    )
    mutations_TSB_df_N = pd.DataFrame(
        N_mutations_df.values.reshape((-1, 16, len(data_f.columns))).sum(axis=1),
        index=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
        columns=data_f.columns,
    )
    # Concatenate the sum row to the DataFrame
    mutations_TSB_df_T = pd.concat(
        [mutations_TSB_df_T, mutations_TSB_df_T.sum().rename("All").to_frame().T]
    )
    mutations_TSB_df_U = pd.concat(
        [mutations_TSB_df_U, mutations_TSB_df_U.sum().rename("All").to_frame().T]
    )
    mutations_TSB_df_N = pd.concat(
        [mutations_TSB_df_N, mutations_TSB_df_N.sum().rename("All").to_frame().T]
    )

    # Reindex the DataFrame to roll the index
    mutations_TSB_df_T = mutations_TSB_df_T.reindex(
        np.roll(mutations_TSB_df_T.index, shift=1)
    )
    mutations_TSB_df_U = mutations_TSB_df_U.reindex(
        np.roll(mutations_TSB_df_U.index, shift=1)
    )
    mutations_TSB_df_N = mutations_TSB_df_N.reindex(
        np.roll(mutations_TSB_df_N.index, shift=1)
    )
    return (
        mutations_df,
        {
            "T": mutations_TSB_df_T,
            "U": mutations_TSB_df_U,
            "N": mutations_TSB_df_N,
        },
    )


def plotSV(
    matrix_path,
    output_path,
    project,
    percentage=False,
    aggregate=False,
    savefig_format="pdf",
    dpi=100,
):
    """Outputs a pdf containing Rearrangement signature plots

    :param matrix_path: path to matrix with 32 channels as rows and samples as columns
    :param output_path: path to output pdf file containing plots
    :param project: name of project
    :param plot_type: output type of plot (default:pdf)
    :param percentage: True if y-axis is displayed as percentage of CNV events, False if displayed as counts (default:False)
    :param aggregate: True if output is a single pdf of counts aggregated across samples(e.g for a given cancer type, y-axis will be counts per sample), False if output is a multi-page pdf of counts for each sample

    # >>> plotSV()

    """

    # inner function to construct plot
    def plot(counts, labels, sample, project, percentage, aggregate=False):
        color_mapping = {
            "del": {
                ">10Mb": "deeppink",
                "1Mb-10Mb": "hotpink",
                "10-100Kb": "lightpink",
                "100Kb-1Mb": "palevioletred",
                "1-10Kb": "lavenderblush",
            },
            "tds": {
                ">10Mb": "saddlebrown",
                "1Mb-10Mb": "sienna",
                "10-100Kb": "sandybrown",
                "100Kb-1Mb": "peru",
                "1-10Kb": "linen",
            },
            "inv": {
                ">10Mb": "rebeccapurple",
                "1Mb-10Mb": "blueviolet",
                "10-100Kb": "plum",
                "100Kb-1Mb": "mediumorchid",
                "1-10Kb": "thistle",
            },
        }

        alpha_dict = dict(enumerate(string.ascii_lowercase))
        x_labels = ["1-10kb", "10-100kb", "100kb-1Mb", "1Mb-10Mb", ">10Mb"]
        super_class = ["clustered", "non-clustered"]
        sub_class = ["del", "tds", "inv", "trans"]
        N = 32
        ticks = [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
        ]
        width = 0.27
        xticks = []
        i = -1  # used to distinguish first bar from the rest
        fig, ax = plt.subplots(figsize=(16, 8))

        # Custom Formatting
        plt.style.use("ggplot")
        plt.rcParams["axes.facecolor"] = "white"
        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_locator(ticker.LinearLocator(5))
        ax.spines["bottom"].set_color("black")
        ax.spines["top"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.spines["left"].set_color("black")
        plt.xlim(xmin=-0.5, xmax=len(labels) - 0.5)
        tmp_max = max(counts)
        plt.ylim(ymax=1.25 * tmp_max)
        # Add light gray horizontal lines at y-ticks
        ax.grid(linestyle="-", linewidth=1, color="#EDEDED", axis="y")

        for count, label in zip(counts, labels):
            categories = label.split("_")
            if len(categories) > 2:
                rearrangement_class = categories[1]
                size_class = categories[2]
            i += 1  # position of bar
            # print (categories)

            if (
                len(categories) == 2
            ):  # clustered translocation or non-clustered translocation
                ax.bar(
                    ticks[i], count, color="dimgray", edgecolor="black"
                )  # translocation only has one color
            else:
                ax.bar(
                    ticks[i],
                    count,
                    color=color_mapping[rearrangement_class][size_class],
                    edgecolor="black",
                )

            xticks.append(ticks[i])
        ax.set_xticks(xticks)
        ax.set_xticklabels(
            x_labels * 3 + [" "] + x_labels * 3 + [" "],
            rotation=90,
            weight="bold",
            fontsize=16,
            fontname="Arial",
            color="black",
        )
        ax.tick_params(labelleft=True, left=False, bottom=False)
        ax.tick_params(axis="y", which="major", pad=0, labelsize=30)

        # Set patch dimensions
        patch_height = 0.05
        patch_width = 2.8
        loh_width = 2.5
        loh_len = 4.8

        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

        #### CLUSTERED PATCHES ####
        ax.add_patch(
            plt.Rectangle(
                (-0.5, 1.095),
                15.9,
                patch_height * 1.5,
                clip_on=False,
                facecolor="gray",
                transform=trans,
            )
        )
        plt.text(
            6,
            1.1125,
            "Clustered",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (-0.5, 1.01),
                loh_len + 0.1,
                patch_height * 1.5,
                clip_on=False,
                facecolor="maroon",
                transform=trans,
            )
        )
        plt.text(
            1.3,
            1.03,
            "Del",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (4.6, 1.01),
                loh_len,
                patch_height * 1.5,
                clip_on=False,
                facecolor="darkorange",
                transform=trans,
            )
        )
        plt.text(
            6.27,
            1.03,
            "Tds",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (9.6, 1.01),
                loh_len,
                patch_height * 1.5,
                clip_on=False,
                facecolor="slateblue",
                transform=trans,
            )
        )
        plt.text(
            11.35,
            1.03,
            "Inv",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (14.6, 1.01),
                0.8,
                patch_height * 1.5,
                clip_on=False,
                facecolor="dimgray",
                transform=trans,
            )
        )
        plt.text(
            14.75,
            1.03,
            "T",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )

        # add vertical black lines
        ax.axvline(x=15.5, color="black", linewidth=1)

        #### NON-CLUSTERED PATCHES ####
        ax.add_patch(
            plt.Rectangle(
                (15.6, 1.095),
                15.9,
                patch_height * 1.5,
                clip_on=False,
                facecolor="black",
                transform=trans,
            )
        )
        plt.text(
            21,
            1.1125,
            "Non-Clustered",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (15.6, 1.01),
                loh_len,
                patch_height * 1.5,
                clip_on=False,
                facecolor="maroon",
                transform=trans,
            )
        )
        plt.text(
            17.35,
            1.03,
            "Del",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (20.6, 1.01),
                loh_len,
                patch_height * 1.5,
                clip_on=False,
                facecolor="darkorange",
                transform=trans,
            )
        )
        plt.text(
            22.25,
            1.03,
            "Tds",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (25.6, 1.01),
                loh_len,
                patch_height * 1.5,
                clip_on=False,
                facecolor="slateblue",
                transform=trans,
            )
        )
        plt.text(
            27.37,
            1.03,
            "Inv",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (30.6, 1.01),
                0.9,
                patch_height * 1.5,
                clip_on=False,
                facecolor="dimgray",
                transform=trans,
            )
        )
        plt.text(
            30.82,
            1.03,
            "T",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )

        # format the set_yticklabels labels
        if percentage:
            tmp_y_labels = [
                "{0:0.1f}%".format(round(x, 1)) for x in ax.get_yticks().tolist()
            ]
        else:
            tmp_y_labels = [round(x, 1) for x in ax.get_yticks().tolist()]
        # ax.yaxis.labelpad = 300

        # set the y-axis ticks and labels
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(
            tmp_y_labels, fontname="Arial", weight="bold", fontsize=16, color="black"
        )

        # y-axis titles
        if aggregate and not percentage:
            ax.set_ylabel(
                "Number of events per sample",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif aggregate and percentage:
            ax.set_ylabel(
                "Percentage of SV's",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif not aggregate and not percentage:
            ax.set_ylabel(
                "Number of events",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif not aggregate and percentage:
            ax.set_ylabel(
                "Percentage(%)",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )

        # TITLE
        if not aggregate:
            plt.text(
                0,
                0.90,
                sample,
                fontsize=20,
                fontname="Arial",
                fontweight="bold",
                color="black",
                transform=trans,
            )
        else:
            plt.text(
                0,
                0.90,
                project,
                fontsize=20,
                fontname="Arial",
                fontweight="bold",
                color="black",
                transform=trans,
            )

        return fig

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    # To reindex the input data
    df = process_input(matrix_path, "32")
    df.reset_index(inplace=True)
    label = df.columns[0]
    labels = df[label]

    figs = {}
    if aggregate:
        num_samples = len(df.columns) - 1
        df["total_count"] = df.select_dtypes(include=[np.number]).sum(axis=1) / num_samples  # NORMALIZE BY # of SAMPLES
        counts = list(df["total_count"])
        if percentage and sum(counts) != 0:
            counts = [(x / sum(counts)) * 100 for x in counts]
        sample = ""
        figs[sample] = plot(counts, labels, sample, project, percentage, aggregate=True)
    else:
        # each column vector in dataframe contains counts for a specific sample
        samples = list(df)[1:]
        for i, (col, sample) in enumerate(zip(df.columns[1:], samples)):
            counts = list(df[col])
            if percentage and sum(counts) != 0:
                counts = [(x / sum(counts)) * 100 for x in counts]
            assert (len(counts)) == 32
            assert (len(labels)) == 32
            figs[sample] = plot(counts, labels, sample, project, percentage)

    return output_results(savefig_format, output_path, project, figs, "SV_32", dpi=dpi)


def plotCNV(
    matrix_path,
    output_path,
    project,
    percentage=False,
    aggregate=False,
    read_from_file=True,
    savefig_format="pdf",
    dpi=100,
):
    """Outputs a pdf containing CNV signature plots

    :param matrix_path: path to matrix generated by CNVMatrixGenerator
    :param output_path: path to output pdf file containing plots
    :param project: name of project
    :param percentage: True if y-axis is displayed as percentage of CNV events, False if displayed as counts (default:False)
    :param aggregate: True if output is a single pdf of counts aggregated across samples(e.g for a given cancer type, y-axis will be counts per sample), False if output is a multi-page pdf of counts for each sample
    :param read_from_file: Retained for API compatibility; input type is detected from matrix_path.
    >>> plotCNV()

    """

    # inner function to construct plot
    def plot(counts, labels, sample, project, percentage, aggregate=False):
        counts_ordered = list()
        labels_ordered = list()
        labels_updated = list()

        # index order will be: homdel, LOH, then het
        for i in range(0, len(labels)):
            labels_ordered.append(str(counts[i]) + ":" + labels[i])
        labels_ordered = pd.Series(labels_ordered)
        l2 = pd.Series(sorted(labels_ordered.values, key=lambda x: x.split(":", 3)[2]))
        lab_homdel = l2[l2.str.contains(":homdel:")]
        lab_homdel.index = [homdel for homdel in range(0, len(lab_homdel))]
        lab_LOH = l2[l2.str.contains(":LOH:")]
        lab_LOH.index = [
            loh for loh in range(len(lab_homdel), len(lab_homdel) + len(lab_LOH))
        ]
        lab_het = l2[l2.str.contains(":het:")]
        lab_het.index = [
            het
            for het in range(
                len(lab_homdel) + len(lab_LOH),
                len(lab_homdel) + len(lab_LOH) + len(lab_het),
            )
        ]
        labels_ordered = pd.concat([lab_homdel, lab_LOH, lab_het])
        for i in labels_ordered:
            tmp_count = i.split(":", 1)
            counts_ordered.append(float(tmp_count[0]))
            labels_updated.append(tmp_count[1])

        labels = pd.Series(labels_updated)
        counts = counts_ordered

        super_class = ["Het", "LOH", "Hom del"]
        hom_del_class = ["0 - 100kb", "100kb - 1Mb", ">1Mb"]
        loh_subclass = ["1", "2", "3-4", "5-8", "9+"]
        het_sub_class = ["2", "3-4", "5-8", "9+"]
        x_labels = ["0 - 100kb", "100kb - 1Mb", "1Mb - 10Mb", "10Mb - 40Mb", ">40Mb"]
        color_mapping = {
            "0:0-100kb": "#F0F8FF",
            "0:100kb-1Mb": "#787CE6",
            "0:>1Mb": "#0000CD",
            "1:0-100kb": "#EBEBEB",
            "1:100kb-1Mb": "#C5C5C5",
            "1:1Mb-10Mb": "#9F9F9F",
            "1:10Mb-40Mb": "#797979",
            "1:>40Mb": "#545454",
            "2:0-100kb": "#F5FFFA",
            "2:100kb-1Mb": "#C0E2C3",
            "2:1Mb-10Mb": "#8BC48E",
            "2:10Mb-40Mb": "#56A858",
            "2:>40Mb": "#228B22",
            "3-4:0-100kb": "#FFF0F5",
            "3-4:100kb-1Mb": "#DEBDEB",
            "3-4:1Mb-10Mb": "#BE8BE1",
            "3-4:10Mb-40Mb": "#9D58D7",
            "3-4:>40Mb": "#7D26CD",
            "5-8:0-100kb": "#FFFAF0",
            "5-8:100kb-1Mb": "#F2DCB3",
            "5-8:1Mb-10Mb": "#E6BF78",
            "5-8:10Mb-40Mb": "#D9A23C",
            "5-8:>40Mb": "#CD8500",
            "9+:0-100kb": "#FFE4E1",
            "9+:100kb-1Mb": "#E2ADBC",
            "9+:1Mb-10Mb": "#C47798",
            "9+:10Mb-40Mb": "#A84074",
            "9+:>40Mb": "#8B0A50",
        }
        colors = ["#0000CD", "#545454", "#228B22", "#7D26CD", "#CD8500", "#8B0A50"]

        N = 48
        ticks = np.arange(N)
        width = 0.27
        xticks = []
        i = -1  # used to distinguish first bar from the rest

        plt.style.use("ggplot")
        plt.rcParams["axes.facecolor"] = "white"

        # create the subplot
        fig, ax = plt.subplots(figsize=(16, 10))
        # Create the plot layout (axis and grid lines)
        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        plt.rcParams["axes.linewidth"] = 1
        ax.yaxis.set_major_locator(ticker.LinearLocator(5))
        ax.spines["bottom"].set_color("black")
        ax.spines["top"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.spines["left"].set_color("black")
        # Add light gray horizontal lines at y-ticks
        ax.grid(linestyle="-", linewidth=1, color="#EDEDED", axis="y")
        plt.xlim(xmin=-0.5, xmax=len(labels) - 0.5)
        tmp_max = max(counts)
        # give buffer of space for title at top left of plot
        plt.ylim(ymax=1.25 * tmp_max)

        for count, label in zip(counts, labels):
            tmp = label.split(":", 2)[1]
            categories = label.split(":")
            cnv_class = categories[0]
            size_class = categories[2]
            # hom del has different color scheme and size classification
            hom_del = False
            if categories[1] == "homdel":
                hom_del = True
            i += 1  # position of bar

            ax.bar(
                ticks[i],
                count,
                color=color_mapping[cnv_class + ":" + size_class],
                edgecolor="black",
                align="center",
            )

            xticks.append(ticks[i])

        # ADD PATCHES AND TEXT
        patch_height = 0.05
        patch_width = 2.8
        loh_width = 2.5
        loh_len = 4.85

        # add vertical black lines
        ax.axvline(x=2.5, color="black", linewidth=1)
        ax.axvline(x=patch_width + loh_len * 5.09, color="black", linewidth=1)

        categories = het_sub_class + loh_subclass + ["Hom" + "\n" + "Del"]
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

        patch_locs = np.arange(0, 45, 5)  # position of patches in data coordinates
        line_locs = (
            []
        )  # for recording positions of top evel patches and separation lines

        # homdel patch
        ax.add_patch(
            plt.Rectangle(
                (-0.5, 1.065),
                2.925,
                patch_height,
                clip_on=False,
                facecolor="#a9a9a9",
                transform=trans,
            )
        )
        ax.add_patch(
            plt.Rectangle(
                (-0.5, 1.01),
                2.925,
                patch_height,
                clip_on=False,
                facecolor=colors[0],
                transform=trans,
            )
        )
        plt.text(
            0.65,
            1.02,
            "0",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        plt.text(
            0.15,
            1.075,
            "HD",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )

        # LOH Patches
        ax.add_patch(
            plt.Rectangle(
                (2.575, 1.065),
                24.825,
                patch_height,
                clip_on=False,
                facecolor="#a9a9a9",
                transform=trans,
            )
        )
        plt.text(
            patch_width + loh_len * 2.25,
            1.075,
            "LOH",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (2.575, 1.01),
                loh_len - 0.025,
                patch_height,
                clip_on=False,
                facecolor=colors[1],
                transform=trans,
            )
        )
        plt.text(
            4.725,
            1.02,
            "1",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (7.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[2],
                transform=trans,
            )
        )
        plt.text(
            9.725,
            1.02,
            "2",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (12.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[3],
                transform=trans,
            )
        )
        plt.text(
            14,
            1.02,
            "3-4",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (17.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[4],
                transform=trans,
            )
        )
        plt.text(
            19.025,
            1.02,
            "5-8",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (22.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[5],
                transform=trans,
            )
        )
        plt.text(
            24.375,
            1.02,
            "9+",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )

        # Heterozygous patches
        ax.add_patch(
            plt.Rectangle(
                (27.55, 1.065),
                loh_len * 4.0938,
                patch_height,
                clip_on=False,
                facecolor="#a9a9a9",
                transform=trans,
            )
        )
        plt.text(
            33.25,
            1.075,
            "Heterozygous",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (27.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[2],
                transform=trans,
            )
        )
        plt.text(
            29.75,
            1.02,
            "2",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (32.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[3],
                transform=trans,
            )
        )
        plt.text(
            34.25,
            1.02,
            "3-4",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (37.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[4],
                transform=trans,
            )
        )
        plt.text(
            39.175,
            1.02,
            "5-8",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )
        ax.add_patch(
            plt.Rectangle(
                (42.55, 1.01),
                loh_len,
                patch_height,
                clip_on=False,
                facecolor=colors[5],
                transform=trans,
            )
        )
        plt.text(
            44.525,
            1.02,
            "9+",
            fontsize=23,
            fontname="Arial",
            fontweight="bold",
            color="white",
            transform=trans,
        )

        # This is the x-axis
        ax.set_xticks(xticks)
        ax.set_xticklabels(
            hom_del_class + x_labels * 9,
            rotation=90,
            weight="bold",
            fontsize=16,
            fontname="Arial",
            color="black",
        )
        ax.tick_params(labelleft=True, left=False, bottom=False)
        ax.tick_params(axis="y", which="major", pad=0, labelsize=60)

        # format the y-axis labels
        if percentage:
            tmp_y_labels = [
                "{0:0.1f}%".format(round(x, 1)) for x in ax.get_yticks().tolist()
            ]
        else:
            tmp_y_labels = [round(x, 1) for x in ax.get_yticks().tolist()]

        # set the y-axis ticks and labels
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(
            tmp_y_labels, fontname="Arial", weight="bold", fontsize=16, color="black"
        )

        # y-axis title
        if aggregate and not percentage:
            ax.set_ylabel(
                "Number of Events Per Sample",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif not aggregate and percentage:
            ax.set_ylabel(
                "Percentage of Copy Number Segments",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif not aggregate and not percentage:
            ax.set_ylabel(
                "Number of Events",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )
        elif aggregate and percentage:
            ax.set_ylabel(
                "Percentage of Copy Number Segments",
                fontsize=24,
                fontname="Arial",
                weight="bold",
                labelpad=15,
                color="black",
            )

        # Add the sample name
        plt.text(
            3,
            0.90,
            sample,
            fontsize=20,
            fontname="Arial",
            fontweight="bold",
            color="black",
            transform=trans,
        )

        return fig

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    # ``process_input`` detects file and in-memory inputs directly. Keep the
    # read_from_file argument for API compatibility, but do not let its default
    # reject a DataFrame before input processing can inspect it.
    if isinstance(matrix_path, str) and not os.path.exists(matrix_path):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), matrix_path
        )

    # To reindex the input data
    df = process_input(matrix_path, "48")
    df.reset_index(inplace=True)
    label = df.columns[0]
    labels = df[label]
    figs = {}
    if aggregate:
        sample_data = df.iloc[:, 1:]
        if sample_data.shape[1] == 0:
            raise ValueError("CNV aggregation requires at least one sample column.")
        try:
            sample_data = sample_data.apply(pd.to_numeric)
        except (TypeError, ValueError) as error:
            raise ValueError("CNV sample columns must contain numeric values.") from error
        num_samples = sample_data.shape[1]
        # Normalize the total event count in each channel by the number of samples.
        df["total_count"] = sample_data.sum(axis=1) / num_samples
        counts = list(df["total_count"])
        if percentage and sum(counts) != 0:
            counts = [(x / sum(counts)) * 100 for x in counts]
        sample = ""
        figs[sample] = plot(
            counts,
            labels,
            sample,
            project,
            percentage,
            aggregate=True,
        )
    else:
        # each column vector in dataframe contains counts for a specific sample
        samples = list(df)[1:]
        for i, (col, sample) in enumerate(zip(df.columns[1:], samples)):
            counts = list(df[col])
            if percentage and sum(counts) != 0:
                counts = [(x / sum(counts)) * 100 for x in counts]
            assert len(counts) == 48
            assert len(labels) == 48
            figs[sample] = plot(
                counts, labels, sample, project, percentage, aggregate=False
            )

    return output_results(savefig_format, output_path, project, figs, "CNV_48", dpi=dpi)


def plotSBS(
    matrix_path,
    output_path,
    project,
    plot_type,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,
):
    """Use an input matrix to create a SBS plot.

    Args:
            matrix_path: The path to a text file or a pandas DataFrame.
            output_path: Path to a directory for saving the output.
            project: Name of unique sample set
            plot_type: Context of the mutational matrix (ie. 96, 288, 384, 1536)
            savefig_format: Format of the output plot. SBS96 and SBS288 support
                pdf, png, and PIL_Image; legacy contexts support pdf only.
            volume: Path to the .pkl file containing the plot template. For Docker.
    Returns:
            Plot of the given input matrix.
    """
    _require_pdf_for_legacy_context(
        savefig_format, plot_type, {"96", "288"}, "plotSBS"
    )
    plot_custom_text = False
    sig_probs = False
    pcawg = False

    # load custom fonts for plotting
    load_custom_fonts()

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    if plot_type == "96":
        data = process_input(matrix_path, plot_type)
        data = reindex_sbs96(data)
        sample_count = 0

        buf = io.BytesIO()
        fig_orig = make_pickle_file(
            context="SBS96", return_plot_template=True, volume=volume
        )
        pickle.dump(fig_orig, buf)

        figs = {}
        buff_list = {}
        ctx = data.index  # [seq[0]+seq[2]+seq[6] for seq in data.index]
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        colorsall = [
            [colors[j] for i in range(int(len(ctx) / 6))] for j in range(6)
        ]
        colors_flat_list = [item for sublist in colorsall for item in sublist]

        for sample in data.columns:
            buf.seek(0)
            figs[sample] = pickle.load(buf)
            panel1 = figs[sample].axes[0]

            total_count = np.sum(data[sample].values)
            x = 0.4
            ymax = 0
            i = 0
            muts = data[sample].values
            if percentage:
                if total_count > 0:
                    plt.bar(
                        np.arange(len(ctx)) + x,
                        muts / total_count * 100,
                        width=0.4,
                        color=colors_flat_list,
                        align="center",
                        zorder=1000,
                    )
                    ymax = np.max(muts / total_count * 100)
                sig_probs = True
            else:
                plt.bar(
                    np.arange(len(ctx)) + x,
                    muts,
                    width=0.4,
                    color=colors_flat_list,
                    align="center",
                    zorder=1000,
                )
                ymax = np.max(muts)

            x = 0.043
            y3 = 0.87
            y = int(ymax * 1.25)
            y2 = y + 2

            if y <= 4:
                y += 4

            while y % 4 != 0:
                y += 1

            y = ymax / 1.025 if ymax > 0 else 4
            ytick_offest = float(y / 3)

            if percentage:
                ylabs = [
                    0,
                    round(ytick_offest, 1),
                    round(ytick_offest * 2, 1),
                    round(ytick_offest * 3, 1),
                    round(ytick_offest * 4, 1),
                ]
                ylabels = [
                    str(0),
                    str(round(ytick_offest, 1)) + "%",
                    str(round(ytick_offest * 2, 1)) + "%",
                    str(round(ytick_offest * 3, 1)) + "%",
                    str(round(ytick_offest * 4, 1)) + "%",
                ]
            else:
                ylabs = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]
                ylabels = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]

            labs = np.arange(0.375, 96.375, 1)

            font_label_size = 30
            if not percentage:
                if int(ylabels[3]) >= 1000:
                    font_label_size = 20

            if percentage:
                if len(ylabels) > 2:
                    font_label_size = 20

            if not percentage:
                ylabels = getylabels(ylabels)

            panel1.set_xlim([0, 96])
            panel1.set_ylim([0, y])
            panel1.set_yticks(ylabs)
            if sig_probs:
                plt.text(
                    0.045,
                    0.75,
                    sample,
                    fontsize=60,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    0.045,
                    0.75,
                    sample + ": " + "{:,}".format(int(total_count)) + " subs",
                    fontsize=60,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

            panel1.set_yticklabels(ylabels, fontsize=font_label_size)
            plt.gca().yaxis.grid(True)
            plt.gca().grid(
                which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
            )
            panel1.set_xlabel("")
            panel1.set_ylabel("")

            custom_text_upper_plot = ""
            try:
                custom_text_upper[sample_count]
            except:
                custom_text_upper = False
            try:
                custom_text_middle[sample_count]
            except:
                custom_text_middle = False
            try:
                custom_text_bottom[sample_count]
            except:
                custom_text_bottom = False

            if custom_text_upper:
                plot_custom_text = True
                if len(custom_text_upper[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False
            if custom_text_middle:
                if len(custom_text_middle[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False

            if plot_custom_text:
                x_pos_custom = 0.98
                if custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = (
                        custom_text_upper[sample_count]
                        + "\n"
                        + custom_text_middle[sample_count]
                    )
                    if custom_text_bottom:
                        custom_text_upper_plot += (
                            "\n" + custom_text_bottom[sample_count]
                        )

                if custom_text_upper and not custom_text_middle:
                    custom_text_upper_plot = custom_text_upper[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                    )

                elif custom_text_upper and custom_text_middle:
                    if not custom_text_bottom:
                        panel1.text(
                            x_pos_custom,
                            0.72,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )
                    else:
                        panel1.text(
                            x_pos_custom,
                            0.68,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                elif not custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = custom_text_middle[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                    )

            if percentage:
                plt.ylabel(
                    "Percentage of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )
            else:
                plt.ylabel(
                    "Number of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )

            panel1.tick_params(
                axis="both",
                which="both",
                bottom=False,
                labelbottom=False,
                left=True,
                labelleft=True,
                right=True,
                labelright=False,
                top=False,
                labeltop=False,
                direction="in",
                length=25,
                colors="lightgray",
                width=2,
            )

            [i.set_color("black") for i in plt.gca().get_yticklabels()]
            sample_count += 1

        return output_results(
            savefig_format, output_path, project, figs, "SBS_96", dpi=dpi
        )

    elif plot_type == "192" or plot_type == "96SB" or plot_type == "384":
        with _open_canonical_legacy_matrix(matrix_path, "384") as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][6] == ">" or first_line[0][3] == ">":
                pcawg = True
            if (
                first_line[0][7] != "]"
                and first_line[0][6] != ">"
                and first_line[0][3] != ">"
            ):
                sys.exit(
                    "The matrix does not match the correct SBS192 format. Please check you formatting and rerun this plotting function."
                )
        file_path = os.path.join(output_path, f"SBS_384_plots_{project}.pdf")
        pp = PdfPages(file_path)
        mutations = OrderedDict()
        try:
            with _open_canonical_legacy_matrix(matrix_path, "384") as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[3:]
                    samples = [x.replace('"', "") for x in samples]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]

                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["C>A"] = OrderedDict()
                    mutations[sample]["C>G"] = OrderedDict()
                    mutations[sample]["C>T"] = OrderedDict()
                    mutations[sample]["T>A"] = OrderedDict()
                    mutations[sample]["T>C"] = OrderedDict()
                    mutations[sample]["T>G"] = OrderedDict()

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        line = [x.replace('"', "") for x in line]
                        nuc = line[2][0] + "[" + line[1] + "]" + line[2][2]
                        bias = line[0][0]
                    else:
                        line = lines.strip().split()
                        nuc = line[0][2:]
                        bias = line[0][0]
                    if bias == "N" or bias == "B":
                        continue
                    else:
                        if pcawg:
                            mut_type = line[1]
                            sample_index = 3
                        else:
                            mut_type = line[0][4:7]
                            sample_index = 1

                        for sample in samples:
                            mutCount = _parse_mutation_value(
                                line[sample_index], percentage, sample, line[0]
                            )
                            if percentage and 0 < mutCount < 1:
                                sig_probs = True
                            if nuc not in mutations[sample][mut_type].keys():
                                mutations[sample][mut_type][nuc] = [0, 0]
                            if bias == "T":
                                mutations[sample][mut_type][nuc][0] = mutCount
                            else:
                                mutations[sample][mut_type][nuc][1] = mutCount
                            sample_index += 1

            sample_count = 0
            for sample in mutations.keys():
                total_count = sum(
                    sum(sum(tsb) for tsb in nuc.values())
                    for nuc in mutations[sample].values()
                )
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(43.93, 9.92))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
                xlabels = []

                x = 0.7
                ymax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0
                for key in mutations[sample]:
                    for seq in mutations[sample][key]:
                        xlabels.append(seq[0] + seq[2] + seq[6])
                        if percentage:
                            if total_count > 0:
                                trans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][0] / total_count * 100,
                                    width=0.75,
                                    color=[1 / 256, 70 / 256, 102 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-transcribed Strand",
                                )
                                x += 0.75
                                untrans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][1] / total_count * 100,
                                    width=0.75,
                                    color=[228 / 256, 41 / 256, 38 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-untranscribed Strand",
                                )
                                x += 0.2475
                                if (
                                    mutations[sample][key][seq][0] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][0]
                                        / total_count
                                        * 100
                                    )
                                if (
                                    mutations[sample][key][seq][1] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][1]
                                        / total_count
                                        * 100
                                    )

                        else:
                            trans = plt.bar(
                                x,
                                mutations[sample][key][seq][0],
                                width=0.75,
                                color=[1 / 256, 70 / 256, 102 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-transcribed Strand",
                            )
                            x += 0.75
                            untrans = plt.bar(
                                x,
                                mutations[sample][key][seq][1],
                                width=0.75,
                                color=[228 / 256, 41 / 256, 38 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-untranscribed Strand",
                            )
                            x += 0.2475
                            if mutations[sample][key][seq][0] > ymax:
                                ymax = mutations[sample][key][seq][0]
                            if mutations[sample][key][seq][1] > ymax:
                                ymax = mutations[sample][key][seq][1]
                        x += 1
                    i += 1

                x = 0.0415
                y3 = 0.87
                y = int(ymax * 1.25)
                x_plot = 0

                yText = y3 + 0.06
                plt.text(
                    0.1,
                    yText,
                    "C>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.255,
                    yText,
                    "C>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.415,
                    yText,
                    "C>T",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.575,
                    yText,
                    "T>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.735,
                    yText,
                    "T>C",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.89,
                    yText,
                    "T>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1

                y = ymax / 1.025

                ytick_offest = float(y / 3)
                for i in range(0, 6, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y3),
                            0.155,
                            0.05,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    panel1.add_patch(
                        plt.Rectangle(
                            (x_plot, 0),
                            32,
                            round(ytick_offest * 4, 1),
                            facecolor=colors[i],
                            zorder=0,
                            alpha=0.25,
                            edgecolor="grey",
                        )
                    )
                    x += 0.1585
                    x_plot += 32

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.750, 192.750, 1)

                font_label_size = 30
                if not percentage:
                    if int(ylabels[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels) > 2:
                        font_label_size = 20

                if not percentage:
                    ylabels = getylabels(ylabels)

                panel1.set_xlim([0, 96])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                count = 0
                m = 0
                for i in range(0, 96, 1):
                    plt.text(
                        i / 101 + 0.0415,
                        0.02,
                        xlabels[i][0],
                        fontsize=30,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 101 + 0.0415,
                        0.044,
                        xlabels[i][1],
                        fontsize=30,
                        color=colors[m],
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 101 + 0.0415,
                        0.071,
                        xlabels[i][2],
                        fontsize=30,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    if count == 16:
                        count = 0
                        m += 1

                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count))
                        + " transcribed subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_bottom:
                    if len(custom_text_bottom[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.84
                    if custom_text_upper and custom_text_bottom:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_bottom[sample_count]
                        )

                    if custom_text_upper and not custom_text_bottom:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_bottom:
                        panel1.text(
                            x_pos_custom,
                            0.72,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif not custom_text_upper and custom_text_bottom:
                        custom_text_upper_plot = custom_text_bottom[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=font_label_size)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                plt.legend(handles=[trans, untrans], prop={"size": 30})
                if percentage:
                    plt.ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors=[0.6, 0.6, 0.6],
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]

                pp.savefig(plot1)
                plt.close()
                sample_count += 1
            pp.close()

        except Exception:
            pdf_path = output_path + "SBS_384_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)
            raise

    elif (
        plot_type == "192_extended"
        or plot_type == "96SB_extended"
        or plot_type == "384_extended"
    ):
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][6] == ">" or first_line[0][3] == ">":
                pcawg = True
            if (
                first_line[0][7] != "]"
                and first_line[0][6] != ">"
                and first_line[0][3] != ">"
            ):
                sys.exit(
                    "The matrix does not match the correct SBS288 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"SBS_384_extended_plots_{project}.pdf")
        pp = PdfPages(file_path)
        mutations = OrderedDict()
        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[3:]
                    samples = [x.replace('"', "") for x in samples]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]

                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["C>A"] = OrderedDict()
                    mutations[sample]["C>G"] = OrderedDict()
                    mutations[sample]["C>T"] = OrderedDict()
                    mutations[sample]["T>A"] = OrderedDict()
                    mutations[sample]["T>C"] = OrderedDict()
                    mutations[sample]["T>G"] = OrderedDict()

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        line = [x.replace('"', "") for x in line]
                        nuc = line[2][0] + "[" + line[1] + "]" + line[2][2]
                        bias = line[0][0]
                    else:
                        line = lines.strip().split()
                        nuc = line[0][2:]
                        bias = line[0][0]

                    if pcawg:
                        mut_type = line[1]
                        sample_index = 3
                    else:
                        mut_type = line[0][4:7]
                        sample_index = 1

                    for sample in samples:
                        if percentage:
                            mutCount = float(line[sample_index])
                            if mutCount < 1 and mutCount > 0:
                                sig_probs = True
                        else:
                            try:
                                mutCount = int(line[sample_index])
                            except:
                                print(
                                    "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                    end="",
                                )

                            # mutCount = int(line[sample_index])
                        if nuc not in mutations[sample][mut_type].keys():
                            mutations[sample][mut_type][nuc] = [0, 0, 0]
                        if bias == "T":
                            mutations[sample][mut_type][nuc][0] = mutCount
                        elif bias == "U":
                            mutations[sample][mut_type][nuc][1] = mutCount
                        else:
                            mutations[sample][mut_type][nuc][2] += mutCount
                        sample_index += 1

            sample_count = 0
            for sample in mutations.keys():
                total_count = sum(
                    sum(sum(tsb) for tsb in nuc.values())
                    for nuc in mutations[sample].values()
                )
                total_trans = sum(
                    sum(sum(tsb[:-1]) for tsb in nuc.values())
                    for nuc in mutations[sample].values()
                )
                total_nontrans = sum(
                    sum(tsb[-1] for tsb in nuc.values())
                    for nuc in mutations[sample].values()
                )
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(43.93, 9.92))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
                xlabels = []

                x = 0.7
                ymax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0
                for key in mutations[sample]:
                    for seq in mutations[sample][key]:
                        xlabels.append(seq[0] + seq[2] + seq[6])
                        if percentage:
                            if total_count > 0:
                                trans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][0] / total_count * 100,
                                    width=0.75,
                                    color=[1 / 256, 70 / 256, 102 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-transcribed Strand",
                                )
                                x += 0.75
                                untrans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][1] / total_count * 100,
                                    width=0.75,
                                    color=[228 / 256, 41 / 256, 38 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-untranscribed Strand",
                                )
                                x += 0.75
                                nontrans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][2] / total_count * 100,
                                    width=0.75,
                                    color="green",
                                    align="center",
                                    zorder=1000,
                                    label="Intergenic",
                                )
                                x += 0.2475
                                if (
                                    mutations[sample][key][seq][0] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][0]
                                        / total_count
                                        * 100
                                    )
                                if (
                                    mutations[sample][key][seq][1] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][1]
                                        / total_count
                                        * 100
                                    )
                                if (
                                    mutations[sample][key][seq][2] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][2]
                                        / total_count
                                        * 100
                                    )

                        else:
                            trans = plt.bar(
                                x,
                                mutations[sample][key][seq][0],
                                width=0.75,
                                color=[1 / 256, 70 / 256, 102 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-transcribed Strand",
                            )
                            x += 0.75
                            untrans = plt.bar(
                                x,
                                mutations[sample][key][seq][1],
                                width=0.75,
                                color=[228 / 256, 41 / 256, 38 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-untranscribed Strand",
                            )
                            x += 0.75
                            nontrans = plt.bar(
                                x,
                                mutations[sample][key][seq][2],
                                width=0.75,
                                color="green",
                                align="center",
                                zorder=1000,
                                label="Intergenic",
                            )
                            x += 0.2475
                            if mutations[sample][key][seq][0] > ymax:
                                ymax = mutations[sample][key][seq][0]
                            if mutations[sample][key][seq][1] > ymax:
                                ymax = mutations[sample][key][seq][1]
                            if mutations[sample][key][seq][2] > ymax:
                                ymax = mutations[sample][key][seq][2]
                        x += 1
                    i += 1

                x = 0.0415
                y3 = 0.87
                y = int(ymax * 1.25)
                x_plot = 0

                yText = y3 + 0.06
                plt.text(
                    0.1,
                    yText,
                    "C>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.255,
                    yText,
                    "C>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.415,
                    yText,
                    "C>T",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.575,
                    yText,
                    "T>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.735,
                    yText,
                    "T>C",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.89,
                    yText,
                    "T>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1

                y = ymax / 1.025

                ytick_offest = float(y / 3)
                for i in range(0, 6, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y3),
                            0.155,
                            0.05,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    panel1.add_patch(
                        plt.Rectangle(
                            (x_plot, 0),
                            44.1,
                            round(ytick_offest * 4, 1),
                            facecolor=colors[i],
                            zorder=0,
                            alpha=0.25,
                            edgecolor="grey",
                        )
                    )
                    x += 0.1585
                    x_plot += 44.1

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.750, 192.750, 1)

                font_label_size = 30
                if not percentage:
                    if int(ylabels[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels) > 2:
                        font_label_size = 20

                if not percentage:
                    ylabels = getylabels(ylabels)

                panel1.set_xlim([0, 264])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                count = 0
                m = 0
                for i in range(0, 96, 1):
                    plt.text(
                        i / 101 + 0.0415,
                        0.02,
                        xlabels[i][0],
                        fontsize=30,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 101 + 0.0415,
                        0.044,
                        xlabels[i][1],
                        fontsize=30,
                        color=colors[m],
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 101 + 0.0415,
                        0.071,
                        xlabels[i][2],
                        fontsize=30,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    if count == 16:
                        count = 0
                        m += 1

                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample
                        + ": "
                        + "{:,}".format(int(total_nontrans))
                        + " non-trans & "
                        + "{:,}".format(int(total_trans))
                        + " trans subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_bottom:
                    if len(custom_text_bottom[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.84
                    if custom_text_upper and custom_text_bottom:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_bottom[sample_count]
                        )

                    if custom_text_upper and not custom_text_bottom:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_bottom:
                        panel1.text(
                            x_pos_custom,
                            0.72,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif not custom_text_upper and custom_text_bottom:
                        custom_text_upper_plot = custom_text_bottom[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=font_label_size)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                plt.legend(
                    handles=[trans, untrans, nontrans],
                    prop={"size": 20},
                    loc="upper right",
                )
                if percentage:
                    plt.ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors=[0.6, 0.6, 0.6],
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]

                pp.savefig(plot1)
                plt.close()
                sample_count += 1
            pp.close()

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "SBS_384_extended_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif plot_type == "6":
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if len(first_line[0]) > 3:
                sys.exit(
                    "The matrix does not match the correct SBS6 format. Please check you formatting and rerun this plotting function."
                )
        file_path = os.path.join(output_path, f"SBS_6_plots_{project}.pdf")
        pp = PdfPages(file_path)

        mutations = OrderedDict()
        total_count = []
        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                samples = first_line.strip().split("\t")
                samples = samples[1:]
                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["C>A"] = 0
                    mutations[sample]["C>G"] = 0
                    mutations[sample]["C>T"] = 0
                    mutations[sample]["T>A"] = 0
                    mutations[sample]["T>C"] = 0
                    mutations[sample]["T>G"] = 0

                for lines in f:
                    line = lines.strip().split()
                    nuc = line[0]
                    mut_type = line[0]
                    sample_index = 1

                    for sample in samples:
                        if percentage:
                            mutCount = float(line[sample_index])
                            if mutCount < 1 and mutCount > 0:
                                sig_probs = True
                        else:
                            try:
                                mutCount = int(line[sample_index])
                            except:
                                print(
                                    "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                    end="",
                                )

                            # mutCount = int(line[sample_index])
                        mutations[sample][mut_type] = mutCount
                        sample_index += 1

            for sample in mutations:
                total_count = sum(mutations[sample].values())
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(15, 10))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])
                xlabels = []

                y = -0.5
                xmax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0
                for key in mutations[sample]:
                    xlabels.append(key)
                    if percentage:
                        if total_count > 0:
                            plt.barh(
                                y,
                                mutations[sample][key] / total_count * 100,
                                height=0.7,
                                color=colors[i],
                                align="center",
                                zorder=1000,
                            )
                            if mutations[sample][key] / total_count * 100 > xmax:
                                xmax = mutations[sample][key] / total_count * 100
                    else:
                        plt.barh(
                            y,
                            mutations[sample][key],
                            height=0.7,
                            color=colors[i],
                            align="center",
                            zorder=1000,
                        )
                        if mutations[sample][key] > xmax:
                            xmax = mutations[sample][key]
                    y -= 1
                    i += 1

                y = 0.043
                x3 = 0.87
                x = int(xmax * 1.1)

                while x % 4 != 0:
                    x += 1
                xtick_offest = int(x / 4)

                if percentage:
                    xlabs = [
                        0,
                        round(xtick_offest, 1),
                        round(xtick_offest * 2, 1),
                        round(xtick_offest * 3, 1),
                        round(xtick_offest * 4, 1),
                    ]
                    xlabels = [
                        str(0),
                        str(round(xtick_offest, 1)) + "%",
                        str(round(xtick_offest * 2, 1)) + "%",
                        str(round(xtick_offest * 3, 1)) + "%",
                        str(round(xtick_offest * 4, 1)) + "%",
                    ]
                else:
                    xlabs = [
                        0,
                        xtick_offest,
                        xtick_offest * 2,
                        xtick_offest * 3,
                        xtick_offest * 4,
                    ]
                    xlabels = [
                        0,
                        xtick_offest,
                        xtick_offest * 2,
                        xtick_offest * 3,
                        xtick_offest * 4,
                    ]

                # if not percentage:
                #   xlabels = ['{:,}'.format(int(x)) for x in xlabels]

                if not percentage:
                    xlabels = ["{:,}".format(int(x)) for x in xlabels]
                    if len(xlabels[-1]) > 3:
                        xlabels_temp = []
                        if len(xlabels[-1]) > 7:
                            for label in xlabels:
                                if len(label) > 7:
                                    xlabels_temp.append(label[0:-8] + "m")
                                elif len(label) > 3:
                                    xlabels_temp.append(label[0:-4] + "k")
                                else:
                                    xlabels_temp.append(label)

                        else:
                            for label in xlabels:
                                if len(label) > 3:
                                    xlabels_temp.append(label[0:-4] + "k")
                                else:
                                    xlabels_temp.append(label)
                        xlabels = xlabels_temp

                ylabs = np.arange(-5.5, 0.5, 1)
                ylabels = ["T>G", "T>C", "T>A", "C>T", "C>G", "C>A"]
                panel1.set_xlim([0, x])
                panel1.set_ylim([-6, 0])
                panel1.set_xticks(xlabs)
                panel1.set_yticks(ylabs)
                panel1.set_xticklabels(xlabels, fontsize=30)
                panel1.set_yticklabels(ylabels, fontsize=30)
                panel1.spines["right"].set_visible(False)
                panel1.spines["top"].set_visible(False)

                if sig_probs:
                    plt.text(
                        0.125,
                        0.9,
                        sample,
                        fontsize=40,
                        fontweight="bold",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.125,
                        0.9,
                        sample + ": " + "{:,}".format(int(total_count)) + " subs",
                        fontsize=40,
                        fontweight="bold",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                if percentage:
                    plt.xlabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.xlabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.set_ylabel("")

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=True,
                    labelbottom=True,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    width=2,
                )

                pp.savefig(plot1)
                plt.close()
            pp.close()

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "SBS_6_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif plot_type == "12" or plot_type == "6SB" or plot_type == "24":
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][1] != ":" or len(first_line[0]) != 5:
                sys.exit(
                    "The matrix does not match the correct SBS192 format. Please check you formatting and rerun this plotting function."
                )
        file_path = os.path.join(output_path, f"SBS_24_plots_{project}.pdf")
        pp = PdfPages(file_path)
        mutations = OrderedDict()

        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                samples = first_line.strip().split("\t")
                samples = samples[1:]
                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["C>A"] = [0, 0]
                    mutations[sample]["C>G"] = [0, 0]
                    mutations[sample]["C>T"] = [0, 0]
                    mutations[sample]["T>A"] = [0, 0]
                    mutations[sample]["T>C"] = [0, 0]
                    mutations[sample]["T>G"] = [0, 0]

                for lines in f:
                    line = lines.strip().split()
                    nuc = line[0][2:]
                    bias = line[0][0]
                    if bias == "N" or bias == "B":
                        continue
                    else:
                        sample_index = 1
                        for sample in samples:
                            if percentage:
                                mutCount = float(line[sample_index])
                                if mutCount < 1 and mutCount > 0:
                                    sig_probs = True
                            else:
                                try:
                                    mutCount = int(line[sample_index])
                                except:
                                    print(
                                        "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                        end="",
                                    )

                                # mutCount = int(line[sample_index])
                            if bias == "T":
                                mutations[sample][nuc][0] = mutCount
                            else:
                                mutations[sample][nuc][1] = mutCount
                            sample_index += 1
            for sample in mutations:
                total_count = sum(sum(tsb) for tsb in mutations[sample].values())
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(15, 10))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])

                y = 12.485
                xmax = 0
                for key in mutations[sample]:
                    if percentage:
                        if total_count > 0:
                            trans = plt.barh(
                                y,
                                mutations[sample][key][0] / total_count * 100,
                                height=0.75,
                                color=[1 / 256, 70 / 256, 102 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-transcribed Strand",
                            )
                            y -= 0.75
                            untrans = plt.barh(
                                y,
                                mutations[sample][key][1] / total_count * 100,
                                height=0.75,
                                color=[228 / 256, 41 / 256, 38 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-untranscribed Strand",
                            )
                            y -= 0.2475
                            if mutations[sample][key][0] / total_count * 100 > xmax:
                                xmax = mutations[sample][key][0] / total_count * 100
                            if mutations[sample][key][1] / total_count * 100 > xmax:
                                xmax = mutations[sample][key][1] / total_count * 100

                    else:
                        trans = plt.barh(
                            y,
                            mutations[sample][key][0],
                            height=0.75,
                            color=[1 / 256, 70 / 256, 102 / 256],
                            align="center",
                            zorder=1000,
                            label="Genic-transcribed Strand",
                        )
                        y -= 0.75
                        untrans = plt.barh(
                            y,
                            mutations[sample][key][1],
                            height=0.75,
                            color=[228 / 256, 41 / 256, 38 / 256],
                            align="center",
                            zorder=1000,
                            label="Genic-untranscribed Strand",
                        )
                        y -= 0.2475
                        if mutations[sample][key][0] > xmax:
                            xmax = mutations[sample][key][0]
                        if mutations[sample][key][1] > xmax:
                            xmax = mutations[sample][key][1]
                    y -= 1

                x = int(xmax * 1.1)

                while x % 4 != 0:
                    x += 1

                xtick_offest = int(x / 4)

                if percentage:
                    xlabs = [
                        0,
                        round(xtick_offest, 1),
                        round(xtick_offest * 2, 1),
                        round(xtick_offest * 3, 1),
                        round(xtick_offest * 4, 1),
                    ]
                    xlabels = [
                        str(0),
                        str(round(xtick_offest, 1)) + "%",
                        str(round(xtick_offest * 2, 1)) + "%",
                        str(round(xtick_offest * 3, 1)) + "%",
                        str(round(xtick_offest * 4, 1)) + "%",
                    ]
                else:
                    xlabs = [
                        0,
                        xtick_offest,
                        xtick_offest * 2,
                        xtick_offest * 3,
                        xtick_offest * 4,
                    ]
                    xlabels = [
                        0,
                        xtick_offest,
                        xtick_offest * 2,
                        xtick_offest * 3,
                        xtick_offest * 4,
                    ]
                if not percentage:
                    xlabels = ["{:,}".format(int(x)) for x in xlabels]
                    if len(xlabels[-1]) > 3:
                        xlabels_temp = []
                        if len(xlabels[-1]) > 7:
                            for label in xlabels:
                                if len(label) > 7:
                                    xlabels_temp.append(label[0:-8] + "m")
                                elif len(label) > 3:
                                    xlabels_temp.append(label[0:-4] + "k")
                                else:
                                    xlabels_temp.append(label)

                        else:
                            for label in xlabels:
                                if len(label) > 3:
                                    xlabels_temp.append(label[0:-4] + "k")
                                else:
                                    xlabels_temp.append(label)
                        xlabels = xlabels_temp

                ylabs = np.arange(2.15, 13, 2)
                ylabels = ["T>G", "T>C", "T>A", "C>T", "C>G", "C>A"]
                panel1.set_xlim([0, x])
                panel1.set_ylim([1.2524, 13.235])
                panel1.set_yticks(ylabs)
                panel1.set_xticks(xlabs)
                panel1.set_xticklabels(xlabels, fontsize=30)
                panel1.set_yticklabels(ylabels, fontsize=30)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                panel1.spines["right"].set_visible(False)
                panel1.spines["top"].set_visible(False)

                if sig_probs:
                    plt.text(
                        0.125,
                        0.9,
                        sample,
                        fontsize=40,
                        fontweight="bold",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.125,
                        0.9,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count))
                        + " transcribed subs",
                        fontsize=40,
                        fontweight="bold",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                if percentage:
                    plt.xlabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.xlabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=True,
                    labelbottom=True,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    width=2,
                )

                plt.legend(handles=[trans, untrans], prop={"size": 25})

                pp.savefig(plot1)
                plt.close()
            pp.close()

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "SBS_24_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif plot_type == "1536":
        with _open_canonical_legacy_matrix(matrix_path, "1536") as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][1] == ">":
                pcawg = True
            if first_line[0][6] != "]" and first_line[0][1] != ">":
                sys.exit(
                    "The matrix does not match the correct SBS1536 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"SBS_1536_plots_{project}.pdf")
        pp = PdfPages(file_path)

        mutations_96 = OrderedDict()
        path_list = matrix_path.split("/")
        extension = path_list[-1].split(".")
        extension = extension[-1]
        matrix_path_96 = (
            "/".join([x for x in path_list[:-1]])
            + "/"
            + project
            + ".SBS96."
            + extension
        )
        mutations = OrderedDict()
        mutations_5 = OrderedDict()
        mutations_3 = OrderedDict()
        max_count = {}
        max_all = {}
        max_5 = {}
        max_3 = {}
        total_count = []
        total_counts = {
            "TT": 0,
            "TG": 0,
            "TC": 0,
            "TA": 0,
            "GT": 0,
            "GG": 0,
            "GC": 0,
            "GA": 0,
            "CT": 0,
            "CG": 0,
            "CC": 0,
            "CA": 0,
            "AT": 0,
            "AG": 0,
            "AC": 0,
            "AA": 0,
        }
        total_counts_5 = {"T": 0, "G": 0, "C": 0, "A": 0}
        total_counts_3 = {"T": 0, "G": 0, "C": 0, "A": 0}

        try:
            with _open_canonical_legacy_matrix(matrix_path, "1536") as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[2:]
                    samples = [x.replace('"', "") for x in samples]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]

                for sample in samples:
                    max_all[sample] = 0
                    max_5[sample] = 0
                    max_3[sample] = 0
                    total_counts[sample] = {
                        "TT": 0,
                        "TG": 0,
                        "TC": 0,
                        "TA": 0,
                        "GT": 0,
                        "GG": 0,
                        "GC": 0,
                        "GA": 0,
                        "CT": 0,
                        "CG": 0,
                        "CC": 0,
                        "CA": 0,
                        "AT": 0,
                        "AG": 0,
                        "AC": 0,
                        "AA": 0,
                    }
                    total_counts_5[sample] = {"T": 0, "G": 0, "C": 0, "A": 0}
                    total_counts_3[sample] = {"T": 0, "G": 0, "C": 0, "A": 0}

                    mutations_96[sample] = OrderedDict()
                    mutations_96[sample]["C>A"] = OrderedDict()
                    mutations_96[sample]["C>G"] = OrderedDict()
                    mutations_96[sample]["C>T"] = OrderedDict()
                    mutations_96[sample]["T>A"] = OrderedDict()
                    mutations_96[sample]["T>C"] = OrderedDict()
                    mutations_96[sample]["T>G"] = OrderedDict()

                    max_count[sample] = 0
                    mutations[sample] = OrderedDict()
                    mutations_5[sample] = OrderedDict()
                    mutations_3[sample] = OrderedDict()

                    mutations[sample]["C>A"] = OrderedDict()
                    mutations_5[sample]["C>A"] = OrderedDict()
                    mutations_3[sample]["C>A"] = OrderedDict()
                    mutations[sample]["C>A"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["C>G"] = OrderedDict()
                    mutations_5[sample]["C>G"] = OrderedDict()
                    mutations_3[sample]["C>G"] = OrderedDict()
                    mutations[sample]["C>G"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["C>T"] = OrderedDict()
                    mutations_5[sample]["C>T"] = OrderedDict()
                    mutations_3[sample]["C>T"] = OrderedDict()
                    mutations[sample]["C>T"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>T"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>T"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>A"] = OrderedDict()
                    mutations_5[sample]["T>A"] = OrderedDict()
                    mutations_3[sample]["T>A"] = OrderedDict()
                    mutations[sample]["T>A"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>C"] = OrderedDict()
                    mutations_5[sample]["T>C"] = OrderedDict()
                    mutations_3[sample]["T>C"] = OrderedDict()
                    mutations[sample]["T>C"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>C"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>C"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>G"] = OrderedDict()
                    mutations_5[sample]["T>G"] = OrderedDict()
                    mutations_3[sample]["T>G"] = OrderedDict()
                    mutations[sample]["T>G"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        line = [x.replace('"', "") for x in line]
                        nuc = line[1][0:2] + "[" + line[0] + "]" + line[1][3:]
                        mut_type = line[0]
                        penta_key = line[1][0] + line[1][-1]
                        tri_key = line[1][1] + line[1][-2]
                        sample_index = 2
                    else:
                        line = lines.strip().split()
                        nuc = line[0]
                        mut_type = line[0][3:6]
                        penta_key = line[0][0] + line[0][-1]
                        tri_key = line[0][1] + line[0][-2]
                        sample_index = 1

                        tri = line[0][1:8]

                    for sample in samples:
                        if tri not in mutations_96[sample][mut_type]:
                            mutations_96[sample][mut_type][tri] = 0
                        mutCount = _parse_mutation_value(
                            line[sample_index], percentage, sample, line[0]
                        )
                        if percentage and 0 < mutCount < 1:
                            sig_probs = True

                        if pcawg:
                            sample_ref = sample_index - 2
                        else:
                            sample_ref = sample_index - 1
                        if mutCount > max_count[samples[sample_ref]]:
                            max_count[samples[sample_ref]] = mutCount

                        if mutCount > max_all[sample]:
                            max_all[sample] = mutCount

                        mutations[sample][mut_type][penta_key][tri_key] = mutCount
                        total_counts[sample][penta_key] += mutCount
                        total_counts_5[sample][penta_key[0]] += mutCount
                        total_counts_3[sample][penta_key[1]] += mutCount
                        penta_key_short = penta_key[0]
                        mutations_5[sample][mut_type][penta_key_short][tri_key] = 0
                        mutations_3[sample][mut_type][penta_key_short][tri_key] = 0
                        mutations_96[sample][mut_type][tri] += mutCount
                        sample_index += 1

            sample_count = 0
            for sample in mutations.keys():
                total_count_sample = sum(
                    sum(nuc.values()) for nuc in mutations_96[sample].values()
                )
                if total_count_sample == 0:
                    continue
                total_count = max_all[sample] * 1.1
                ratio = total_count / total_count_sample
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(43.93, 22.5))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.03, 0.0677, 0.92, 0.267])  # 1536 context panel
                panel2 = plt.axes([0.03, 0.67, 0.92, 0.247])  # 96 context panel
                panel3 = plt.axes([0.03, 0.35, 0.92, 0.1335])  # 3' context
                panel4 = plt.axes([0.03, 0.5, 0.92, 0.1335])  # 5' context
                xlabels = []
                ylabels = []
                ylabels_5 = []
                ylabels_3 = []

                # Set up all of the color maps
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                colors_heat = [
                    np.linspace(56 / 255, 255 / 255, 5),
                    np.linspace(66 / 255, 225 / 255, 5),
                    np.linspace(157 / 255, 40 / 255, 5),
                ]
                colors_heat_compact = [
                    np.linspace(56 / 255, 255 / 255, 5),
                    np.linspace(66 / 255, 225 / 255, 5),
                    np.linspace(157 / 255, 40 / 255, 5),
                ]

                # Plot the 1536 matrix and collect the relevant info for the 96, 5' and 3' plots
                i = 0
                x_pos = 0
                x_inter = 0
                for key in mutations[sample]:
                    y_pos = 15
                    for penta in mutations[sample][key]:
                        key_5 = penta[0]
                        key_3 = penta[1]
                        ylabels.append(penta[0] + "---" + penta[1])
                        for tri in mutations[sample][key][penta]:
                            tri_nuc = tri[0] + "[" + key + "]" + tri[1]
                            normalized = mutations_96[sample][key][tri_nuc]
                            try:
                                mut_count = int(
                                    int(
                                        20
                                        * round(
                                            float(
                                                mutations[sample][key][penta][tri]
                                                / total_count_sample
                                                / ratio
                                                * 100
                                            )
                                        )
                                        / 20
                                    )
                                    / 20
                                )
                                mutations_5[sample][key][key_5][tri] += float(
                                    mutations[sample][key][penta][tri]
                                )
                                mutations_3[sample][key][key_3][tri] += float(
                                    mutations[sample][key][penta][tri]
                                )
                                if mutations_5[sample][key][key_5][tri] > max_5[sample]:
                                    max_5[sample] = mutations_5[sample][key][key_5][tri]
                                if mutations_3[sample][key][key_3][tri] > max_3[sample]:
                                    max_3[sample] = mutations_3[sample][key][key_3][tri]
                            except:
                                mut_count = 0
                            xlabels.append(tri[0] + "-" + tri[1])
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat[0][mut_count],
                                    colors_heat[1][mut_count],
                                    colors_heat[2][mut_count],
                                ),
                            )
                            panel1.add_patch(rectangle)
                            x_pos += 1
                        y_pos -= 1
                        x_pos = x_inter

                    x_inter += 17
                    x_pos = x_inter
                    i += 1

                # Plot 5' and 3' context matrices
                x_pos = 0
                x_inter = 0
                total_count_5 = max_5[sample] * 1.1
                total_count_3 = max_3[sample] * 1.1
                ratio_5 = total_count_5 / total_count_sample
                ratio_3 = total_count_3 / total_count_sample
                ratio_total = max(ratio_5, ratio_3)
                for key in mutations_5[sample]:
                    y_pos = 3
                    for penta in mutations_5[sample][key]:
                        ylabels_5.append(penta + "---N")
                        ylabels_3.append("N---" + penta)
                        for tri in mutations_5[sample][key][penta]:
                            mut_count = int(
                                int(
                                    20
                                    * round(
                                        float(mutations_5[sample][key][penta][tri])
                                        / total_count_sample
                                        / ratio_total
                                        * 100
                                    )
                                    / 20
                                )
                                / 20
                            )
                            mut_count_3 = int(
                                int(
                                    20
                                    * round(
                                        float(mutations_3[sample][key][penta][tri])
                                        / total_count_sample
                                        / ratio_total
                                        * 100
                                    )
                                    / 20
                                )
                                / 20
                            )
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat_compact[0][mut_count],
                                    colors_heat_compact[1][mut_count],
                                    colors_heat_compact[2][mut_count],
                                ),
                            )

                            panel4.add_patch(rectangle)
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat[0][mut_count_3],
                                    colors_heat[1][mut_count_3],
                                    colors_heat[2][mut_count_3],
                                ),
                            )
                            panel3.add_patch(rectangle)
                            x_pos += 1
                        y_pos -= 1
                        x_pos = x_inter
                    x_inter += 17
                    x_pos = x_inter
                    i += 1

                # Plot the 96 bar plot
                x = 0.5
                ymax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0

                for key in mutations_96[sample]:
                    for seq in mutations_96[sample][key]:
                        xlabels.append(seq[0] + seq[2] + seq[6])
                        if percentage:
                            if total_count_sample > 0:
                                panel2.bar(
                                    x,
                                    mutations_96[sample][key][seq]
                                    / total_count_sample
                                    * 100,
                                    width=0.5,
                                    color=colors[i],
                                    align="center",
                                    zorder=1000,
                                )
                                if (
                                    mutations_96[sample][key][seq]
                                    / total_count_sample
                                    * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations_96[sample][key][seq]
                                        / total_count_sample
                                        * 100
                                    )
                        else:
                            panel2.bar(
                                x,
                                mutations_96[sample][key][seq],
                                width=0.5,
                                color=colors[i],
                                align="center",
                                zorder=1000,
                            )
                            if mutations_96[sample][key][seq] > ymax:
                                ymax = mutations_96[sample][key][seq]
                        x += 1
                    x += 1
                    i += 1

                # scale bar for bottom 1536 plot
                y_grad = 0.267 / len(colors_heat[0])
                y_start = 0.0677
                for l in range(0, len(colors_heat[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.96, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat[0][l],
                            colors_heat[1][l],
                            colors_heat[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat[0][l],
                            colors_heat[1][l],
                            colors_heat[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad

                # scale bar for top 1536 plot
                y_grad = 0.1335 / len(colors_heat_compact[0])
                y_start = 0.5
                for l in range(0, len(colors_heat_compact[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.96, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad

                # scale bar for middle 1536 plot
                y_grad = 0.1335 / len(colors_heat_compact[0])
                y_start = 0.35
                for l in range(0, len(colors_heat_compact[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.96, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad
                y_tick_grad = max_count[sample] / 2

                # scale numbers for bottom 1536 plot
                plt.text(
                    0.9825,
                    0.0677,
                    "0",
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.2012,
                    str(ratio / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.325,
                    str(ratio)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )

                # scale numbers for top 1536 plot
                plt.text(
                    0.9825,
                    0.5,
                    "0",
                    fontsize=20,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.56675,
                    str(ratio_total / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.625,
                    str(ratio_total)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )

                # scale numbers for middle 1536 plot
                plt.text(
                    0.9825,
                    0.35,
                    "0",
                    fontsize=20,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.41675,
                    str(ratio_total / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.9825,
                    0.475,
                    str(ratio_total)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                y = int(ymax * 1.25)

                x = 0.033
                y3 = 0.92
                for i in range(0, 6, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y3),
                            0.143,
                            0.03,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    x += 0.154

                # Plot the labels for above the SBS96 bar plot
                y3 = 0.9
                yText = y3 + 0.06
                plt.text(
                    0.085,
                    yText,
                    "C>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.24,
                    yText,
                    "C>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.395,
                    yText,
                    "C>T",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.552,
                    yText,
                    "T>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.705,
                    yText,
                    "T>C",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.86,
                    yText,
                    "T>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                # Set up the parameters for the yscale ticks and labels
                y = ymax / 1.025
                ytick_offest = float(y / 3)

                ylabels_96 = []
                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels_96 = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels_96 = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                # Adjust fontsize if required due to large numbers on y-axis
                font_label_size = 25
                if not percentage:
                    if int(ylabels_96[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels_96) > 2:
                        font_label_size = 20

                if not percentage:
                    ylabels_96 = getylabels(ylabels_96)

                # Set all panel limits
                panel1.set_xlim([0, 101])
                panel1.set_ylim([0, 16])
                panel2.set_xlim([0, 101])
                panel2.set_yticks(ylabs)
                panel1.set_yticks([])
                panel1.set_xticks([])
                panel2.set_xticks([])
                panel4.set_xlim([0, 101])
                panel4.set_ylim([0, 4])
                panel3.set_xlim([0, 101])
                panel3.set_ylim([0, 4])
                panel4.set_yticks([])
                panel3.set_yticks([])
                panel3.set_xticks([])
                panel4.set_xticks([])

                m = 0
                count = 0
                x_letter = 0
                for i in range(0, 96, 1):
                    # Bottom labels
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.04,
                        xlabels[i][0],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.05,
                        xlabels[i][1],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.06,
                        xlabels[i][2],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )

                    # Top labels
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.64,
                        xlabels[i][0],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.65,
                        xlabels[i][1],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.032,
                        0.66,
                        xlabels[i][2],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    x_letter += 0.92
                    if (i + 1) % 16 == 0 and i != 0:
                        x_letter += 0.92
                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 bottom plot
                m = 0
                count = 0
                y_letter = 5.2
                for i in range(0, 16, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.2675

                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 top matrix plot
                m = 0
                count = 0
                y_letter = 9.85
                for i in range(0, 4, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels_5[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels_5[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels_5[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.53

                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 middle matrix plot
                m = 0
                count = 0
                y_letter = 7.45
                for i in range(0, 4, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels_3[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels_3[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels_3[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.53

                    if count == 16:
                        count = 0
                        m += 1

                # Set up all parameters for custom text in upper righthand corner
                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.94
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel2.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel2.text(
                                x_pos_custom,
                                0.86,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel2.text(
                                x_pos_custom,
                                0.835,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel2.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                # Plot the sample name in upper left corner of plot
                if sig_probs:
                    panel2.text(
                        0.04,
                        0.875,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    panel2.text(
                        0.04,
                        0.875,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count_sample))
                        + " subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                # Finish setting all axis grids and tick parameters
                panel2.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                panel1.set_yticklabels([])
                panel1.set_xticklabels([])
                panel2.set_xticklabels([])

                if percentage:
                    panel2.set_ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=28,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    panel2.set_ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=28,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.axis("off")
                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel2.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel3.axis("off")
                panel3.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel4.axis("off")
                panel4.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )

                # Change axis depending if percentage parameter is selected
                if percentage:
                    panel2.set_yticklabels(
                        ylabels_96, fontsize=font_label_size - 1, color="black"
                    )
                else:
                    panel2.set_yticklabels(
                        ylabels_96, fontsize=font_label_size, color="black"
                    )

                pp.savefig(plot1)
                plt.close()
                sample_count += 1
            pp.close()
        except Exception:
            pdf_path = output_path + "SBS_1536_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)
            raise

    elif plot_type == "4608":
        with _open_canonical_legacy_matrix(matrix_path, "4608") as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][1] == ">":
                pcawg = True
            if first_line[0][8] != "]" and first_line[0][6] != ">":
                sys.exit(
                    "The matrix does not match the correct SBS4608 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"SBS_4608_plots_{project}.pdf")
        pp = PdfPages(file_path)

        path_list = matrix_path.split("/")
        extension = path_list[-1].split(".")
        extension = extension[-1]
        mutations = OrderedDict()
        mutations_5 = OrderedDict()
        mutations_3 = OrderedDict()
        mutations_TSB = OrderedDict()
        mutations_96 = OrderedDict()
        max_count = {}
        max_all = {}
        max_5 = {}
        max_3 = {}
        total_count = []
        total_counts = {
            "TT": 0,
            "TG": 0,
            "TC": 0,
            "TA": 0,
            "GT": 0,
            "GG": 0,
            "GC": 0,
            "GA": 0,
            "CT": 0,
            "CG": 0,
            "CC": 0,
            "CA": 0,
            "AT": 0,
            "AG": 0,
            "AC": 0,
            "AA": 0,
        }
        total_counts_5 = {"T": 0, "G": 0, "C": 0, "A": 0}
        total_counts_3 = {"T": 0, "G": 0, "C": 0, "A": 0}

        # try:
        if True:
            with _open_canonical_legacy_matrix(matrix_path, "4608") as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[2:]
                    samples = [x.replace('"', "") for x in samples]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]

                for sample in samples:
                    max_all[sample] = 0
                    max_5[sample] = 0
                    max_3[sample] = 0
                    total_counts[sample] = {
                        "TT": 0,
                        "TG": 0,
                        "TC": 0,
                        "TA": 0,
                        "GT": 0,
                        "GG": 0,
                        "GC": 0,
                        "GA": 0,
                        "CT": 0,
                        "CG": 0,
                        "CC": 0,
                        "CA": 0,
                        "AT": 0,
                        "AG": 0,
                        "AC": 0,
                        "AA": 0,
                    }
                    total_counts_5[sample] = {"T": 0, "G": 0, "C": 0, "A": 0}
                    total_counts_3[sample] = {"T": 0, "G": 0, "C": 0, "A": 0}

                    mutations_96[sample] = OrderedDict()
                    mutations_96[sample]["C>A"] = OrderedDict()
                    mutations_96[sample]["C>G"] = OrderedDict()
                    mutations_96[sample]["C>T"] = OrderedDict()
                    mutations_96[sample]["T>A"] = OrderedDict()
                    mutations_96[sample]["T>C"] = OrderedDict()
                    mutations_96[sample]["T>G"] = OrderedDict()

                    max_count[sample] = 0
                    mutations[sample] = OrderedDict()
                    mutations_5[sample] = OrderedDict()
                    mutations_3[sample] = OrderedDict()

                    mutations[sample]["C>A"] = OrderedDict()
                    mutations_5[sample]["C>A"] = OrderedDict()
                    mutations_3[sample]["C>A"] = OrderedDict()
                    mutations[sample]["C>A"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["C>G"] = OrderedDict()
                    mutations_5[sample]["C>G"] = OrderedDict()
                    mutations_3[sample]["C>G"] = OrderedDict()
                    mutations[sample]["C>G"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["C>T"] = OrderedDict()
                    mutations_5[sample]["C>T"] = OrderedDict()
                    mutations_3[sample]["C>T"] = OrderedDict()
                    mutations[sample]["C>T"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["C>T"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["C>T"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>A"] = OrderedDict()
                    mutations_5[sample]["T>A"] = OrderedDict()
                    mutations_3[sample]["T>A"] = OrderedDict()
                    mutations[sample]["T>A"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>A"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>C"] = OrderedDict()
                    mutations_5[sample]["T>C"] = OrderedDict()
                    mutations_3[sample]["T>C"] = OrderedDict()
                    mutations[sample]["T>C"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>C"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>C"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations[sample]["T>G"] = OrderedDict()
                    mutations_5[sample]["T>G"] = OrderedDict()
                    mutations_3[sample]["T>G"] = OrderedDict()
                    mutations[sample]["T>G"] = {
                        "TT": OrderedDict(),
                        "TG": OrderedDict(),
                        "TC": OrderedDict(),
                        "TA": OrderedDict(),
                        "GT": OrderedDict(),
                        "GG": OrderedDict(),
                        "GC": OrderedDict(),
                        "GA": OrderedDict(),
                        "CT": OrderedDict(),
                        "CG": OrderedDict(),
                        "CC": OrderedDict(),
                        "CA": OrderedDict(),
                        "AT": OrderedDict(),
                        "AG": OrderedDict(),
                        "AC": OrderedDict(),
                        "AA": OrderedDict(),
                    }
                    mutations_5[sample]["T>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }
                    mutations_3[sample]["T>G"] = {
                        "T": OrderedDict(),
                        "G": OrderedDict(),
                        "C": OrderedDict(),
                        "A": OrderedDict(),
                    }

                    mutations_TSB[sample] = OrderedDict()
                    mutations_TSB[sample]["All"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>A"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>G"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>T"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>A"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>C"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>G"] = OrderedDict({"T": 0, "U": 0, "N": 0})

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        line = [x.replace('"', "") for x in line]
                        nuc = line[2][0:2] + "[" + line[0] + "]" + line[1][3:]
                        mut_type = line[0]
                        penta_key = line[1][0] + line[1][-1]
                        tri_key = line[1][1] + line[1][-2]
                        sample_index = 2
                    else:
                        line = lines.strip().split()
                        nuc = line[0]
                        mut_type = line[0][5:8]
                        penta_key = line[0][2] + line[0][-1]
                        tri_key = line[0][3] + line[0][-2]
                        sample_index = 1
                        tsb = nuc[0]
                        tri = line[0][3:10]

                    for sample in samples:
                        if tri not in mutations_96[sample][mut_type]:
                            mutations_96[sample][mut_type][tri] = 0
                        mutCount = _parse_mutation_value(
                            line[sample_index], percentage, sample, line[0]
                        )
                        if percentage and 0 < mutCount < 1:
                            sig_probs = True

                        if pcawg:
                            sample_ref = sample_index - 2
                        else:
                            sample_ref = sample_index - 1
                        if mutCount > max_count[samples[sample_ref]]:
                            max_count[samples[sample_ref]] = mutCount

                        if mutCount > max_all[sample]:
                            max_all[sample] = mutCount

                        mutations[sample][mut_type][penta_key][tri_key] = mutCount
                        total_counts[sample][penta_key] += mutCount
                        total_counts_5[sample][penta_key[0]] += mutCount
                        total_counts_3[sample][penta_key[1]] += mutCount
                        penta_key_short = penta_key[0]
                        mutations_5[sample][mut_type][penta_key_short][tri_key] = 0
                        mutations_3[sample][mut_type][penta_key_short][tri_key] = 0
                        mutations_96[sample][mut_type][tri] += mutCount
                        mutations_TSB[sample][mut_type][tsb] += mutCount
                        mutations_TSB[sample]["All"][tsb] += mutCount
                        sample_index += 1

            sample_count = 0
            for sample in mutations.keys():
                total_count_sample = sum(
                    sum(nuc.values()) for nuc in mutations_96[sample].values()
                )
                if total_count_sample == 0:
                    continue
                total_count = max_all[sample] * 1.1
                ratio = total_count / total_count_sample
                plt.rcParams["axes.linewidth"] = 2
                # plot1 = plt.figure(figsize=(43.93,22.5))
                plot1 = plt.figure(figsize=(43.93, 18))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.03, 0.0677, 0.76, 0.267])  # 1536 context panel
                panel2 = plt.axes([0.03, 0.67, 0.76, 0.247])  # 96 context panel
                panel3 = plt.axes([0.03, 0.35, 0.76, 0.1335])  # 3' context
                panel4 = plt.axes([0.03, 0.5, 0.76, 0.1335])  # 5' context
                panel5 = plt.axes([0.82, 0.67, 0.17, 0.247])  # SBS288 (TSB) context
                xlabels = []
                ylabels = []
                ylabels_5 = []
                ylabels_3 = []

                # Set up all of the color maps
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                colors_heat = [
                    np.linspace(56 / 255, 255 / 255, 5),
                    np.linspace(66 / 255, 225 / 255, 5),
                    np.linspace(157 / 255, 40 / 255, 5),
                ]
                colors_heat_compact = [
                    np.linspace(56 / 255, 255 / 255, 5),
                    np.linspace(66 / 255, 225 / 255, 5),
                    np.linspace(157 / 255, 40 / 255, 5),
                ]

                # Plot the 1536 matrix and collect the relevant info for the 96, 5' and 3' plots
                i = 0
                x_pos = 0
                x_inter = 0
                for key in mutations[sample]:
                    y_pos = 15
                    for penta in mutations[sample][key]:
                        key_5 = penta[0]
                        key_3 = penta[1]
                        ylabels.append(penta[0] + "---" + penta[1])
                        for tri in mutations[sample][key][penta]:
                            tri_nuc = tri[0] + "[" + key + "]" + tri[1]
                            normalized = mutations_96[sample][key][tri_nuc]
                            try:
                                mut_count = int(
                                    int(
                                        20
                                        * round(
                                            float(
                                                mutations[sample][key][penta][tri]
                                                / total_count_sample
                                                / ratio
                                                * 100
                                            )
                                        )
                                        / 20
                                    )
                                    / 20
                                )
                                mutations_5[sample][key][key_5][tri] += float(
                                    mutations[sample][key][penta][tri]
                                )
                                mutations_3[sample][key][key_3][tri] += float(
                                    mutations[sample][key][penta][tri]
                                )
                                if mutations_5[sample][key][key_5][tri] > max_5[sample]:
                                    max_5[sample] = mutations_5[sample][key][key_5][tri]
                                if mutations_3[sample][key][key_3][tri] > max_3[sample]:
                                    max_3[sample] = mutations_3[sample][key][key_3][tri]
                            except:
                                mut_count = 0
                            xlabels.append(tri[0] + "-" + tri[1])
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat[0][mut_count],
                                    colors_heat[1][mut_count],
                                    colors_heat[2][mut_count],
                                ),
                            )
                            panel1.add_patch(rectangle)
                            x_pos += 1
                        y_pos -= 1
                        x_pos = x_inter

                    x_inter += 17
                    x_pos = x_inter
                    i += 1

                # Plot 5' and 3' context matrices
                x_pos = 0
                x_inter = 0
                total_count_5 = max_5[sample] * 1.1
                total_count_3 = max_3[sample] * 1.1
                ratio_5 = total_count_5 / total_count_sample
                ratio_3 = total_count_3 / total_count_sample
                ratio_total = max(ratio_5, ratio_3)
                for key in mutations_5[sample]:
                    y_pos = 3
                    for penta in mutations_5[sample][key]:
                        ylabels_5.append(penta + "---N")
                        ylabels_3.append("N---" + penta)
                        for tri in mutations_5[sample][key][penta]:
                            mut_count = int(
                                int(
                                    20
                                    * round(
                                        float(mutations_5[sample][key][penta][tri])
                                        / total_count_sample
                                        / ratio_total
                                        * 100
                                    )
                                    / 20
                                )
                                / 20
                            )
                            mut_count_3 = int(
                                int(
                                    20
                                    * round(
                                        float(mutations_3[sample][key][penta][tri])
                                        / total_count_sample
                                        / ratio_total
                                        * 100
                                    )
                                    / 20
                                )
                                / 20
                            )
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat_compact[0][mut_count],
                                    colors_heat_compact[1][mut_count],
                                    colors_heat_compact[2][mut_count],
                                ),
                            )

                            panel4.add_patch(rectangle)
                            rectangle = mplpatches.Rectangle(
                                (x_pos, y_pos),
                                1,
                                1,
                                linewidth=1,
                                facecolor=(
                                    colors_heat[0][mut_count_3],
                                    colors_heat[1][mut_count_3],
                                    colors_heat[2][mut_count_3],
                                ),
                            )
                            panel3.add_patch(rectangle)
                            x_pos += 1
                        y_pos -= 1
                        x_pos = x_inter
                    x_inter += 17
                    x_pos = x_inter
                    i += 1

                # Plot the 96 bar plot
                x = 0.5
                ymax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0

                for key in mutations_96[sample]:
                    for seq in mutations_96[sample][key]:
                        xlabels.append(seq[0] + seq[2] + seq[6])
                        if percentage:
                            if total_count_sample > 0:
                                panel2.bar(
                                    x,
                                    mutations_96[sample][key][seq]
                                    / total_count_sample
                                    * 100,
                                    width=0.5,
                                    color=colors[i],
                                    align="center",
                                    zorder=1000,
                                )
                                if (
                                    mutations_96[sample][key][seq]
                                    / total_count_sample
                                    * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations_96[sample][key][seq]
                                        / total_count_sample
                                        * 100
                                    )
                        else:
                            panel2.bar(
                                x,
                                mutations_96[sample][key][seq],
                                width=0.5,
                                color=colors[i],
                                align="center",
                                zorder=1000,
                            )
                            if mutations_96[sample][key][seq] > ymax:
                                ymax = mutations_96[sample][key][seq]
                        x += 1
                    x += 1
                    i += 1

                # scale bar for bottom 1536 plot
                y_grad = 0.267 / len(colors_heat[0])
                y_start = 0.0677
                for l in range(0, len(colors_heat[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.797, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat[0][l],
                            colors_heat[1][l],
                            colors_heat[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat[0][l],
                            colors_heat[1][l],
                            colors_heat[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad

                # scale bar for top 1536 plot
                y_grad = 0.1335 / len(colors_heat_compact[0])
                y_start = 0.5
                for l in range(0, len(colors_heat_compact[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.797, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad

                # scale bar for middle 1536 plot
                y_grad = 0.1335 / len(colors_heat_compact[0])
                y_start = 0.35
                for l in range(0, len(colors_heat_compact[0]), 1):
                    rectangle = mplpatches.Rectangle(
                        (0.797, y_start),
                        0.02,
                        y_grad,
                        linewidth=1,
                        facecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                        transform=plt.gcf().transFigure,
                        clip_on=False,
                        edgecolor=(
                            colors_heat_compact[0][l],
                            colors_heat_compact[1][l],
                            colors_heat_compact[2][l],
                        ),
                    )
                    panel1.add_patch(rectangle)
                    y_start += y_grad
                y_tick_grad = max_count[sample] / 2

                # scale numbers for bottom 1536 plot
                plt.text(
                    0.82,
                    0.0677,
                    "0",
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.2012,
                    str(ratio / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.325,
                    str(ratio)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )

                # scale numbers for top 1536 plot
                plt.text(
                    0.82,
                    0.5,
                    "0",
                    fontsize=20,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.56675,
                    str(ratio_total / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.625,
                    str(ratio_total)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )

                # scale numbers for middle 1536 plot
                plt.text(
                    0.82,
                    0.35,
                    "0",
                    fontsize=20,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.41675,
                    str(ratio_total / 2)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.82,
                    0.475,
                    str(ratio_total)[:5],
                    fontsize=15,
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                y = int(ymax * 1.25)

                x = 0.033
                y3 = 0.92
                for i in range(0, 6, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y3),
                            0.117,
                            0.03,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    x += 0.1275

                # Plot the labels for above the SBS96 bar plot
                y3 = 0.9
                yText = y3 + 0.06
                plt.text(
                    0.075,
                    yText,
                    "C>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.2025,
                    yText,
                    "C>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.33,
                    yText,
                    "C>T",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.455,
                    yText,
                    "T>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.585,
                    yText,
                    "T>C",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.71,
                    yText,
                    "T>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                # Set up the parameters for the yscale ticks and labels
                y = ymax / 1.025
                ytick_offest = float(y / 3)

                ylabels_96 = []
                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels_96 = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels_96 = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                # Adjust fontsize if required due to large numbers on y-axis
                font_label_size = 25
                if not percentage:
                    if int(ylabels_96[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels_96) > 2:
                        font_label_size = 20

                if not percentage:
                    ylabels_96 = getylabels(ylabels_96)

                # Set all panel limits
                panel1.set_xlim([0, 101])
                panel1.set_ylim([0, 16])
                panel2.set_xlim([0, 101])
                panel2.set_yticks(ylabs)
                panel1.set_yticks([])
                panel1.set_xticks([])
                panel2.set_xticks([])
                panel4.set_xlim([0, 101])
                panel4.set_ylim([0, 4])
                panel3.set_xlim([0, 101])
                panel3.set_ylim([0, 4])
                panel4.set_yticks([])
                panel3.set_yticks([])
                panel3.set_xticks([])
                panel4.set_xticks([])

                # x-axis 1536 bottom plot
                m = 0
                count = 0
                x_letter = 0
                for i in range(0, 96, 1):
                    # Bottom labels
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.04,
                        xlabels[i][0],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.05,
                        xlabels[i][1],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.06,
                        xlabels[i][2],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )

                    # Top labels
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.64,
                        xlabels[i][0],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.65,
                        xlabels[i][1],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        x_letter / 101 + 0.031,
                        0.66,
                        xlabels[i][2],
                        fontsize=25,
                        color="black",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    x_letter += 0.76

                    if (i + 1) % 16 == 0 and i != 0:
                        x_letter += 0.76
                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 bottom plot
                m = 0
                count = 0
                y_letter = 5.2
                for i in range(0, 16, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.2675

                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 top matrix plot
                m = 0
                count = 0
                y_letter = 9.85
                for i in range(0, 4, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels_5[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels_5[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels_5[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.53

                    if count == 16:
                        count = 0
                        m += 1

                # y-axis 1536 middle matrix plot
                m = 0
                count = 0
                y_letter = 7.45
                for i in range(0, 4, 1):
                    plt.text(
                        0.003,
                        y_letter / 16,
                        ylabels_3[i][0],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.008,
                        y_letter / 16 + 0,
                        ylabels_3[i][1:4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        0.022,
                        y_letter / 16 + 0,
                        ylabels_3[i][4],
                        fontsize=25,
                        color="black",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    y_letter -= 0.53

                    if count == 16:
                        count = 0
                        m += 1

                # Set up all parameters for custom text in upper righthand corner
                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                x_pos_custom = 0.85
                y_pos_upper = 0.125
                y_pos_middle = 0.10
                y_pos_bottom = 0.075
                zorder_value = 10

                fontsize_custom = 25
                fontweight_custom = "bold"
                fontcolor_custom = "black"
                fontname_custom = "Arial"
                ha_custom = "left"

                if custom_text_upper:
                    plt.text(
                        x_pos_custom,
                        y_pos_upper,
                        custom_text_upper[sample_count],
                        fontsize=fontsize_custom,
                        weight=fontweight_custom,
                        color=fontcolor_custom,
                        fontname=fontname_custom,
                        transform=plt.gcf().transFigure,
                        ha=ha_custom,
                        zorder=zorder_value,
                    )

                if custom_text_middle:
                    plt.text(
                        x_pos_custom,
                        y_pos_middle,
                        custom_text_middle[sample_count],
                        fontsize=fontsize_custom,
                        weight=fontweight_custom,
                        color=fontcolor_custom,
                        fontname=fontname_custom,
                        transform=plt.gcf().transFigure,
                        ha=ha_custom,
                        zorder=zorder_value,
                    )

                if custom_text_bottom:
                    plt.text(
                        x_pos_custom,
                        y_pos_bottom,
                        custom_text_bottom[sample_count],
                        fontsize=fontsize_custom,
                        weight=fontweight_custom,
                        color=fontcolor_custom,
                        fontname=fontname_custom,
                        transform=plt.gcf().transFigure,
                        ha=ha_custom,
                        zorder=zorder_value,
                    )

                # Plot the sample name in upper left corner of plot
                if sig_probs:
                    panel2.text(
                        0.04,
                        0.875,
                        sample,
                        fontsize=50,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    panel2.text(
                        0.04,
                        0.875,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count_sample))
                        + " subs",
                        fontsize=50,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                # Finish setting all axis grids and tick parameters
                panel2.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                panel1.set_yticklabels([])
                panel1.set_xticklabels([])
                panel2.set_xticklabels([])

                if percentage:
                    panel2.set_ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=26,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    panel2.set_ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=26,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.axis("off")
                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel2.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel3.axis("off")
                panel3.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )
                panel4.axis("off")
                panel4.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=False,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="white",
                    width=2,
                )

                # Change axis depending if percentage parameter is selected
                if percentage:
                    panel2.set_yticklabels(
                        ylabels_96, fontsize=font_label_size - 1, color="black"
                    )
                else:
                    panel2.set_yticklabels(
                        ylabels_96, fontsize=font_label_size, color="black"
                    )

                yp2 = 28
                labels = []
                y2max = 0
                tsbColors = [
                    [1 / 256, 70 / 256, 102 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    "green",
                ]
                for mut in mutations_TSB[sample]:
                    labels.append(mut)
                    i = 0
                    for tsb in mutations_TSB[sample][mut]:
                        if tsb == "T":
                            label = "Genic-transcribed"
                        elif tsb == "U":
                            label = "Genic-untranscribed"
                        else:
                            label = "Intergenic"
                        if percentage:
                            if total_count_sample > 0:
                                panel5.barh(
                                    yp2,
                                    mutations_TSB[sample][mut][tsb]
                                    / total_count_sample
                                    * 100,
                                    color=tsbColors[i],
                                    label=label,
                                )
                                if (
                                    mutations_TSB[sample][mut][tsb]
                                    / total_count_sample
                                    * 100
                                    > y2max
                                ):
                                    y2max = (
                                        mutations_TSB[sample][mut][tsb]
                                        / total_count_sample
                                        * 100
                                    )
                        else:
                            if mutations_TSB[sample][mut][tsb] > y2max:
                                y2max = mutations_TSB[sample][mut][tsb]

                            panel5.barh(
                                yp2,
                                mutations_TSB[sample][mut][tsb],
                                color=tsbColors[i],
                                label=label,
                            )
                        yp2 -= 1
                        i += 1
                    yp2 -= 1
                y = int(y2max * 1.1)
                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)
                if percentage:
                    xlabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    xlabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    xlabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    xlabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                if not percentage:
                    xlabels = getxlabels(xlabels)

                panel5.spines["right"].set_visible(False)
                panel5.spines["top"].set_visible(False)
                labels.reverse()
                panel5.set_yticks([3, 7, 11, 15, 19, 23, 27])
                panel5.set_yticklabels(
                    labels, fontsize=30, fontname="Arial", weight="bold"
                )
                panel5.set_xticklabels(xlabels, fontsize=30)
                panel5.set_xticks(xlabs)
                handles, labels = panel5.get_legend_handles_labels()
                panel5.legend(handles[:3], labels[:3], loc="best", prop={"size": 20})

                pp.savefig(plot1)
                plt.close()
                sample_count += 1
            pp.close()

    elif plot_type == "288":
        plot_custom_text = False
        sig_probs = False
        pcawg = False

        data = process_input(matrix_path, plot_type)

        sample_count = 0

        data, tsb_mats = reindex_sbs288(data)
        buf = io.BytesIO()
        fig_orig = make_pickle_file(
            context="SBS288", return_plot_template=True, volume=volume
        )
        pickle.dump(fig_orig, buf)

        figs = {}
        buff_list = {}
        ctx = data.index
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        colorsall = [
            [colors[j] for i in range(int(len(ctx) / 6))] for j in range(6)
        ]
        colors_flat_list = [item for sublist in colorsall for item in sublist]

        for sample in data.columns:
            buf.seek(0)
            figs[sample] = pickle.load(buf)
            panel1 = figs[sample].axes[0]
            panel2 = figs[sample].axes[1]

            total_count = np.sum(data[sample].values)
            muts = data[sample].values
            x = 0.4
            ymax = 0
            if percentage:
                if total_count > 0:
                    panel1.bar(
                        np.arange(len(ctx)) + x,
                        muts / total_count * 100,
                        width=0.4,
                        color=colors_flat_list,
                        align="center",
                        zorder=1000,
                    )
                    ymax = np.max(muts / total_count * 100)
                sig_probs = True
            else:
                panel1.bar(
                    np.arange(len(ctx)) + x,
                    muts,
                    width=0.4,
                    color=colors_flat_list,
                    align="center",
                    zorder=1000,
                )
                ymax = np.max(muts)

            y = int(ymax * 1.25)
            if y <= 4:
                y += 4

            while y % 4 != 0:
                y += 1
            y = ymax / 1.025 if ymax > 0 else 4
            ytick_offest = float(y / 3)

            if percentage:
                ylabs = [
                    0,
                    round(ytick_offest, 1),
                    round(ytick_offest * 2, 1),
                    round(ytick_offest * 3, 1),
                    round(ytick_offest * 4, 1),
                ]
                ylabels = [
                    str(0),
                    str(round(ytick_offest, 1)) + "%",
                    str(round(ytick_offest * 2, 1)) + "%",
                    str(round(ytick_offest * 3, 1)) + "%",
                    str(round(ytick_offest * 4, 1)) + "%",
                ]

            else:
                ylabs = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]
                ylabels = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]

            font_label_size = 30
            if not percentage:
                if int(ylabels[3]) >= 1000:
                    font_label_size = 20

            if percentage:
                if len(ylabels) > 2:
                    font_label_size = 20

            labs = np.arange(0.375, 96.375, 1)

            if not percentage:
                ylabels = spplt.getylabels(ylabels)

            panel1.set_xlim([0, 96])
            panel1.set_ylim([0, y])

            panel1.set_yticks(ylabs)

            if sig_probs:
                plt.text(
                    0.045,
                    0.75,
                    sample,
                    fontsize=60,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    0.045,
                    0.75,
                    sample + ": " + "{:,}".format(int(total_count)) + " subs",
                    fontsize=60,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

            custom_text_upper_plot = ""
            try:
                custom_text_upper[sample_count]
            except:
                custom_text_upper = False
            try:
                custom_text_middle[sample_count]
            except:
                custom_text_middle = False
            try:
                custom_text_bottom[sample_count]
            except:
                custom_text_bottom = False

            if custom_text_upper:
                plot_custom_text = True
                if len(custom_text_upper[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False
            if custom_text_middle:
                if len(custom_text_middle[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False

            if plot_custom_text:
                x_pos_custom = 0.73
                if custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = (
                        custom_text_upper[sample_count]
                        + "\n"
                        + custom_text_middle[sample_count]
                    )
                    if custom_text_bottom:
                        custom_text_upper_plot += (
                            "\n" + custom_text_bottom[sample_count]
                        )

                if custom_text_upper and not custom_text_middle:
                    custom_text_upper_plot = custom_text_upper[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                        zorder=1,
                    )

                elif custom_text_upper and custom_text_middle:
                    if not custom_text_bottom:
                        panel1.text(
                            x_pos_custom,
                            0.72,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )
                    else:
                        panel1.text(
                            x_pos_custom,
                            0.68,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                elif not custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = custom_text_middle[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                    )

            panel1.set_yticklabels(ylabels, fontsize=font_label_size)
            panel1.yaxis.grid(True)
            panel1.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
            panel1.set_xlabel("")
            panel1.set_ylabel("")

            if percentage:
                panel1.set_ylabel(
                    "Percentage of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )
            else:
                panel1.set_ylabel(
                    "Number of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )

            panel1.tick_params(
                axis="both",
                which="both",
                bottom=False,
                labelbottom=False,
                left=True,
                labelleft=True,
                right=True,
                labelright=False,
                top=False,
                labeltop=False,
                direction="in",
                length=25,
                colors="lightgray",
                width=2,
            )

            [i.set_color("black") for i in panel1.get_yticklabels()]

            yp2 = 28
            labels = []
            y2max = 0
            tsbColors = [
                [1 / 256, 70 / 256, 102 / 256],
                [228 / 256, 41 / 256, 38 / 256],
                "green",
            ]

            if percentage and total_count > 0:
                y2max = (
                    np.max(
                        [
                            tsb_mats["T"][sample]["All"],
                            tsb_mats["U"][sample]["All"],
                            tsb_mats["N"][sample]["All"],
                        ]
                    )
                    / total_count
                    * 100
                )

                panel2.barh(
                    range(28, 1, -4),
                    tsb_mats["T"][sample].values / total_count * 100,
                    color=tsbColors[0],
                    label="Genic-transcribed",
                )
                panel2.barh(
                    range(27, 1, -4),
                    tsb_mats["U"][sample].values / total_count * 100,
                    color=tsbColors[1],
                    label="Genic-untranscribed",
                )
                panel2.barh(
                    range(26, 1, -4),
                    tsb_mats["N"][sample].values / total_count * 100,
                    color=tsbColors[2],
                    label="Intergenic",
                )
            elif percentage:
                y2max = 0
                zeros = np.zeros(len(tsb_mats["T"][sample]))
                panel2.barh(
                    range(28, 1, -4),
                    zeros,
                    color=tsbColors[0],
                    label="Genic-transcribed",
                )
                panel2.barh(
                    range(27, 1, -4),
                    zeros,
                    color=tsbColors[1],
                    label="Genic-untranscribed",
                )
                panel2.barh(
                    range(26, 1, -4),
                    zeros,
                    color=tsbColors[2],
                    label="Intergenic",
                )
            else:
                y2max = np.max(
                    [
                        tsb_mats["T"][sample]["All"],
                        tsb_mats["U"][sample]["All"],
                        tsb_mats["N"][sample]["All"],
                    ]
                )
                panel2.barh(
                    range(28, 1, -4),
                    tsb_mats["T"][sample].values,
                    color=tsbColors[0],
                    label="Genic-transcribed",
                )
                panel2.barh(
                    range(27, 1, -4),
                    tsb_mats["U"][sample].values,
                    color=tsbColors[1],
                    label="Genic-untranscribed",
                )
                panel2.barh(
                    range(26, 1, -4),
                    tsb_mats["N"][sample].values,
                    color=tsbColors[2],
                    label="Intergenic",
                )

            labels = list(tsb_mats["T"][sample].index)

            y = int(y2max * 1.1)
            if y <= 4:
                y += 4
            while y % 4 != 0:
                y += 1
            ytick_offest = int(y / 4)

            if percentage:
                xlabs = [
                    0,
                    round(ytick_offest, 1),
                    round(ytick_offest * 2, 1),
                    round(ytick_offest * 3, 1),
                    round(ytick_offest * 4, 1),
                ]
                xlabels = [
                    str(0),
                    str(round(ytick_offest, 1)) + "%",
                    str(round(ytick_offest * 2, 1)) + "%",
                    str(round(ytick_offest * 3, 1)) + "%",
                    str(round(ytick_offest * 4, 1)) + "%",
                ]
            else:
                xlabs = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]
                xlabels = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]

            if not percentage:
                xlabels = spplt.getxlabels(xlabels)

            panel2.spines["right"].set_visible(False)
            panel2.spines["top"].set_visible(False)
            labels.reverse()
            panel2.set_yticks([3, 7, 11, 15, 19, 23, 27])
            panel2.set_yticklabels(
                labels, fontsize=30, fontname="Arial", weight="bold"
            )
            panel2.set_xticks(xlabs)
            panel2.set_xticklabels(xlabels, fontsize=30)
            handles, labels = panel2.get_legend_handles_labels()
            panel2.legend(handles[:3], labels[:3], loc="best", prop={"size": 30})
            sample_count += 1

        return output_results(
            savefig_format, output_path, project, figs, "SBS_288", dpi=dpi
        )

    elif plot_type == "288_Normalized":
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            if first_line[0][6] == ">" or first_line[0][3] == ">":
                pcawg = True
            if (
                first_line[0][7] != "]"
                and first_line[0][6] != ">"
                and first_line[0][3] != ">"
            ):
                sys.exit(
                    "The matrix does not match the correct SBS288 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"SBS_288_Normalized_plots_{project}.pdf")
        pp = PdfPages(file_path)

        mutations = OrderedDict()
        mutations_TSB = OrderedDict()
        total_count = []
        if True:
            # try:
            with open(matrix_path) as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[2:]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]

                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["C>A"] = OrderedDict()
                    mutations[sample]["C>G"] = OrderedDict()
                    mutations[sample]["C>T"] = OrderedDict()
                    mutations[sample]["T>A"] = OrderedDict()
                    mutations[sample]["T>C"] = OrderedDict()
                    mutations[sample]["T>G"] = OrderedDict()

                    mutations_TSB[sample] = OrderedDict()
                    mutations_TSB[sample]["All"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>A"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>G"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["C>T"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>A"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>C"] = OrderedDict({"T": 0, "U": 0, "N": 0})
                    mutations_TSB[sample]["T>G"] = OrderedDict({"T": 0, "U": 0, "N": 0})

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        mut_type = line[0]
                        nuc = line[1][0] + "[" + mut_type + "]" + line[1][2]
                        sample_index = 2
                    else:
                        line = lines.strip().split()
                        nuc = line[0]
                        mut_type = line[0][4:7]
                        sample_index = 1
                        tsb = nuc[0]

                    for sample in samples:
                        if percentage:
                            mutCount = float(line[sample_index])
                            if mutCount < 1 and mutCount > 0:
                                sig_probs = True
                        else:
                            try:
                                mutCount = int(line[sample_index])
                            except:
                                print(
                                    "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                    end="",
                                )
                        if nuc[2:] not in mutations[sample][mut_type]:
                            mutations[sample][mut_type][nuc[2:]] = mutCount
                        else:
                            mutations[sample][mut_type][nuc[2:]] += mutCount
                        mutations_TSB[sample][mut_type][tsb] += mutCount
                        mutations_TSB[sample]["All"][tsb] += mutCount
                        sample_index += 1

            sample_count = 0

            for sample in mutations.keys():
                total_count = sum(
                    sum(nuc.values()) for nuc in mutations[sample].values()
                )
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(43.93, 9.92))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.04, 0.09, 0.7, 0.77])
                panel2 = plt.axes([0.77, 0.09, 0.21, 0.77])
                xlabels = []

                x = 0.4
                ymax = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [1 / 256, 1 / 256, 1 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [203 / 256, 202 / 256, 202 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [236 / 256, 199 / 256, 197 / 256],
                ]
                i = 0
                for key in mutations[sample]:
                    for seq in mutations[sample][key]:
                        xlabels.append(seq[0] + seq[2] + seq[6])
                        if percentage:
                            if total_count > 0:
                                panel1.bar(
                                    x,
                                    mutations[sample][key][seq] / total_count * 100,
                                    width=0.4,
                                    color=colors[i],
                                    align="center",
                                    zorder=1000,
                                )
                                if (
                                    mutations[sample][key][seq] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq] / total_count * 100
                                    )
                        else:
                            panel1.bar(
                                x,
                                mutations[sample][key][seq],
                                width=0.4,
                                color=colors[i],
                                align="center",
                                zorder=1000,
                            )
                            if mutations[sample][key][seq] > ymax:
                                ymax = mutations[sample][key][seq]
                        x += 1
                    i += 1

                x = 0.043
                y3 = 0.87
                y = int(ymax * 1.25)
                y2 = y + 2
                for i in range(0, 6, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y3),
                            0.11,
                            0.05,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    x += 0.117

                yText = y3 + 0.06

                plt.text(
                    0.082,
                    yText,
                    "C>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.1975,
                    yText,
                    "C>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.315,
                    yText,
                    "C>T",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.43,
                    yText,
                    "T>A",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.55,
                    yText,
                    "T>C",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.665,
                    yText,
                    "T>G",
                    fontsize=55,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                # ytick_offest = int(y/4)
                y = ymax / 1.025
                ytick_offest = float(y / 3)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                font_label_size = 30
                if not percentage:
                    if int(ylabels[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels) > 2:
                        font_label_size = 20

                labs = np.arange(0.375, 96.375, 1)

                if not percentage:
                    ylabels = getylabels(ylabels)

                panel1.set_xlim([0, 96])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                count = 0
                m = 0
                for i in range(0, 96, 1):
                    plt.text(
                        i / 137 + 0.04,
                        0.02,
                        xlabels[i][0],
                        fontsize=25,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 137 + 0.04,
                        0.044,
                        xlabels[i][1],
                        fontsize=25,
                        color=colors[m],
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        fontweight="bold",
                        transform=plt.gcf().transFigure,
                    )
                    plt.text(
                        i / 137 + 0.04,
                        0.071,
                        xlabels[i][2],
                        fontsize=25,
                        color="gray",
                        rotation="vertical",
                        verticalalignment="center",
                        fontname="Courier New",
                        transform=plt.gcf().transFigure,
                    )
                    count += 1
                    if count == 16:
                        count = 0
                        m += 1

                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample + ": " + "{:,}".format(int(total_count)) + " subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.73
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                            zorder=1,
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.72,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.68,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=font_label_size)
                panel1.yaxis.grid(True)
                panel1.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                if percentage:
                    panel1.set_ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    panel1.set_ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="lightgray",
                    width=2,
                )

                [i.set_color("black") for i in panel1.get_yticklabels()]

                yp2 = 25.5
                labels = []
                y2max = 0
                tsbColors = [
                    [1 / 256, 70 / 256, 102 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    "green",
                ]
                for mut in mutations_TSB[sample]:
                    labels.append(mut)
                    if True:
                        totalMuts = sum(
                            [
                                mutations_TSB[sample][mut][x]
                                for x in mutations_TSB[sample][mut]
                            ]
                        )
                        panel2.barh(
                            yp2,
                            mutations_TSB[sample][mut]["N"] / totalMuts * 100,
                            1.25,
                            color=tsbColors[2],
                            label="Intergenic",
                            zorder=0,
                        )
                        panel2.barh(
                            yp2,
                            (
                                mutations_TSB[sample][mut]["N"]
                                + mutations_TSB[sample][mut]["T"]
                                + mutations_TSB[sample][mut]["U"]
                            )
                            / totalMuts
                            * 100,
                            1.25,
                            color="purple",
                            label="Genic-transcribed Regions",
                            zorder=-1,
                        )

                        yp2 -= 1.25

                        transcribedTotal = (
                            mutations_TSB[sample][mut]["T"]
                            + mutations_TSB[sample][mut]["U"]
                        )
                        panel2.barh(
                            yp2,
                            mutations_TSB[sample][mut]["T"] / transcribedTotal * 100,
                            1.25,
                            color=tsbColors[0],
                            label="Transcribed strand",
                            zorder=0,
                        )
                        panel2.barh(
                            yp2,
                            (
                                mutations_TSB[sample][mut]["T"]
                                + mutations_TSB[sample][mut]["U"]
                            )
                            / transcribedTotal
                            * 100,
                            1.25,
                            color=tsbColors[1],
                            label="Genic-untranscribed strand",
                            zorder=-1,
                        )

                        yp2 -= 1.25
                    yp2 -= 1.25
                panel2.grid(
                    which="major",
                    axis="x",
                    ls="--",
                    color=[0.95, 0.95, 0.95],
                    zorder=10000,
                    lw=2,
                )
                xlabels = ["0%", "25%", "50%", "75", "100%"]
                xlabs = [0, 25, 50, 75, 100]
                panel2.spines["right"].set_visible(False)
                panel2.spines["top"].set_visible(False)
                labels.reverse()
                panel2.set_yticks([2.5, 6.25, 10, 13.75, 17.75, 21.25, 25])
                panel2.set_yticklabels(
                    labels, fontsize=30, fontname="Arial", weight="bold"
                )
                panel2.set_xticks(xlabs)
                panel2.set_xticklabels(xlabels, fontsize=30)
                handles, labels = panel2.get_legend_handles_labels()
                panel2.legend(
                    handles[:4],
                    labels[:4],
                    loc="upper right",
                    prop={"size": 15},
                    bbox_to_anchor=(0.95, 1.15),
                )
                pp.savefig(plot1)
                plt.close()
                sample_count += 1

            pp.close()

    else:
        print(
            "Error: The function plotSBS does not support plot_type",
            plot_type,
            "so no plot has been generated."
        )


def plotID(
    matrix_path,
    output_path,
    project,
    plot_type,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,
):
    _require_pdf_for_legacy_context(
        savefig_format, plot_type, {"94", "ID94", "94ID", "83"}, "plotID"
    )
    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    # load custom fonts for plotting
    load_custom_fonts()

    plot_custom_text = False
    sig_probs = False
    pcawg = False
    if (
        plot_type == "94"
        or plot_type == "ID94"
        or plot_type == "94ID"
        or plot_type == "83"
    ):
        data = process_input(matrix_path, plot_type)

        try:
            sample_count = 0
            buf = io.BytesIO()
            fig_orig = make_pickle_file(
                context="ID83", return_plot_template=True, volume=volume
            )
            pickle.dump(fig_orig, buf)
            figs = {}
            colors = [
                [253 / 256, 190 / 256, 111 / 256],
                [255 / 256, 128 / 256, 2 / 256],
                [176 / 256, 221 / 256, 139 / 256],
                [54 / 256, 161 / 256, 46 / 256],
                [253 / 256, 202 / 256, 181 / 256],
                [252 / 256, 138 / 256, 106 / 256],
                [241 / 256, 68 / 256, 50 / 256],
                [188 / 256, 25 / 256, 26 / 256],
                [208 / 256, 225 / 256, 242 / 256],
                [148 / 256, 196 / 256, 223 / 256],
                [74 / 256, 152 / 256, 201 / 256],
                [23 / 256, 100 / 256, 171 / 256],
                [226 / 256, 226 / 256, 239 / 256],
                [182 / 256, 182 / 256, 216 / 256],
                [134 / 256, 131 / 256, 189 / 256],
                [98 / 256, 64 / 256, 155 / 256],
            ]
            ctx = data.index
            xlabels = [
                i.split(":")[0] + i.split(":")[1] + i.split(":")[2]
                for i in ctx.to_list()
            ]
            xlables_set = [
                "1DelC",
                "1DelT",
                "1InsC",
                "1InsT",
                "2DelR",
                "3DelR",
                "4DelR",
                "5DelR",
                "2InsR",
                "3InsR",
                "4InsR",
                "5InsR",
                "2DelM",
                "3DelM",
                "4DelM",
                "5DelM",
            ]
            colors_idx = copy.deepcopy(xlabels)
            for ii in range(0, len(xlables_set)):
                colors_idx = [ii if x == xlables_set[ii] else x for x in colors_idx]

            colors_flat_list = [colors[i] for i in colors_idx]

            for sample in data.columns:  # mutations.keys():
                buf.seek(0)
                figs[sample] = pickle.load(buf)
                panel1 = figs[sample].axes[0]
                muts = data[sample].values
                total_count = np.sum(muts)
                x = 0.4
                ymax = 0

                if percentage:
                    if total_count > 0:
                        plt.bar(
                            np.arange(len(ctx)) + x,
                            muts / total_count * 100,
                            width=0.4,
                            color=colors_flat_list,
                            align="center",
                            zorder=1000,
                        )
                        ymax = np.max(muts / total_count * 100)
                    sig_probs = True
                else:
                    plt.bar(
                        np.arange(len(ctx)) + x,
                        muts,
                        width=0.4,
                        color=colors_flat_list,
                        align="center",
                        zorder=1000,
                    )
                    ymax = np.max(muts)

                x = 0.0475
                y_top = 0.827
                y_bottom = 0.114
                y = int(ymax * 1.25)
                y2 = y + 2

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.375, 83.375, 1)

                if not percentage:
                    ylabels = spplt.getylabels(ylabels)

                panel1.set_xlim([0, 83])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                if sig_probs:
                    plt.text(
                        0.0475,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.0475,
                        0.75,
                        sample + ": " + "{:,}".format(int(total_count)) + " indels",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.95
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.72,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.68,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=30)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                if percentage:
                    plt.ylabel(
                        "Percentage of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="gray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                sample_count += 1

            return output_results(
                savefig_format, output_path, project, figs, "ID_83", dpi=dpi
            )
        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "ID_83_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif (
        plot_type == "INDEL_simple"
        or plot_type == "simple_INDEL"
        or plot_type == "ID_simple"
        or plot_type == "simple_ID"
        or plot_type == "28"
    ):
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            mutation_type = first_line[0]
            mutation_type_list = mutation_type.split(":")
            if len(mutation_type_list) != 4:
                sys.exit(
                    "The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function."
                )
        file_path = os.path.join(output_path, f"ID_simple_plots_{project}.pdf")
        pp = PdfPages(file_path)

        indel_types = [
            "1:Del:C:1",
            "1:Del:C:2",
            "1:Del:C:3",
            "1:Del:C:4",
            "1:Del:C:5",
            "1:Del:C:6" "1:Del:T:1",
            "1:Del:T:2",
            "1:Del:T:3",
            "1:Del:T:4",
            "1:Del:T:5",
            "1:Del:T:6" "1:Ins:C:0",
            "1:Ins:C:1",
            "1:Ins:C:2",
            "1:Ins:C:3",
            "1:Ins:C:4",
            "1:Ins:C:5",
            "1:Ins:T:0",
            "1:Ins:T:1",
            "1:Ins:T:2",
            "1:Ins:T:3",
            "1:Ins:T:4",
            "1:Ins:T:5",
            "long_Del",
            "long_Ins",
            "MH",
            "complex",
        ]

        mutations = OrderedDict()

        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                samples = first_line.strip().split("\t")
                samples = samples[1:]
                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["1DelC"] = [0, 0, 0, 0, 0, 0]
                    mutations[sample]["1DelT"] = [0, 0, 0, 0, 0, 0]
                    mutations[sample]["1InsC"] = [0, 0, 0, 0, 0, 0]
                    mutations[sample]["1InsT"] = [0, 0, 0, 0, 0, 0]
                    mutations[sample]["long_Del"] = [0]
                    mutations[sample]["long_Ins"] = [0]
                    mutations[sample]["MH"] = [0]
                    mutations[sample]["complex"] = [0]

                for lines in f:
                    line = lines.strip().split()
                    categories = line[0].split(":")
                    if len(categories) < 2:
                        mut_type = categories[0]
                        repeat_size = 0
                    else:
                        mut_type = categories[0] + categories[1] + categories[2]
                        repeat_size = int(categories[3])
                    sample_index = 1

                    for sample in samples:
                        # if mut_type in mutations[sample].keys():
                        if percentage:
                            mutCount = float(line[sample_index])
                            if mutCount < 1 and mutCount > 0:
                                sig_probs = True
                        else:
                            try:
                                mutCount = int(line[sample_index])
                            except:
                                print(
                                    "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                    end="",
                                )

                            # mutCount = int(line[sample_index])
                        mutations[sample][mut_type][repeat_size] = mutCount

                        # else:
                        #   if percentage:
                        #       mutCount = float(line[sample_index])
                        #   else:
                        #       mutCount = int(line[sample_index])
                        #   if int(mut_type[0]) > 1:
                        #       repeat_size = 0
                        #       if categories[2] == 'M':
                        #           mut_type = 'MH'
                        #           mutations[sample][mut_type][repeat_size] += mutCount
                        #       else:
                        #           if categories[1] == 'Del':
                        #               mut_type = 'Del'
                        #               mutations[sample][mut_type][repeat_size] += mutCount
                        #           else:
                        #               mut_type = 'Ins'
                        #               mutations[sample][mut_type][repeat_size] += mutCount

                        # continue
                        sample_index += 1

            for sample in mutations:
                total_count = sum(sum(nuc) for nuc in mutations[sample].values())
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(15, 13))
                plt.rc("axes", edgecolor="black")
                panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])
                xlabels = []

                x = 0.4
                ymax = 0
                colors = [
                    [253 / 256, 190 / 256, 111 / 256],
                    [255 / 256, 128 / 256, 2 / 256],
                    [176 / 256, 221 / 256, 139 / 256],
                    [54 / 256, 161 / 256, 46 / 256],
                    # [188/256,25/256,26/256],
                    [23 / 256, 100 / 256, 171 / 256],
                    [98 / 256, 64 / 256, 155 / 256],
                    [98 / 256, 64 / 256, 155 / 256],
                ]

                i = 0
                for key in mutations[sample]:
                    l = 1
                    for seq in mutations[sample][key]:
                        xlabels.append(l)
                        if percentage:
                            if total_count > 0:
                                plt.bar(
                                    x,
                                    seq / total_count * 100,
                                    width=0.4,
                                    color=colors[i],
                                    align="center",
                                    zorder=1000,
                                )
                                if seq / total_count * 100 > ymax:
                                    ymax = seq / total_count * 100
                        else:
                            plt.bar(
                                x,
                                seq,
                                width=0.4,
                                color=colors[i],
                                align="center",
                                zorder=1000,
                            )
                            if seq > ymax:
                                ymax = seq
                        x += 1
                        l += 1
                    if i < 4:
                        i += 1
                x = 0.126
                y_top = 0.9
                y_bottom = 0.075
                y = int(ymax * 1.25)
                y2 = y + 2

                for i in range(0, 4, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y_top),
                            0.154,
                            0.037,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y_bottom),
                            0.154,
                            0.037,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    x += 0.1715

                x -= 0.001
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_top),
                        0.098,
                        0.037,
                        facecolor=colors[i + 1],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_bottom),
                        0.098,
                        0.037,
                        facecolor=colors[i + 1],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )

                yText = y_top + 0.0055
                plt.text(
                    0.185,
                    yText,
                    "C",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.36,
                    yText,
                    "T",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.53,
                    yText,
                    "C",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.705,
                    yText,
                    "T",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    color="white",
                    transform=plt.gcf().transFigure,
                )

                yText_labels_top = yText + 0.045
                yText_labels_bottom = y_bottom - 0.03
                yText_labels_bottom_sec = yText_labels_bottom - 0.025

                plt.text(
                    0.2,
                    yText_labels_top,
                    "1bp Deletion",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.54,
                    yText_labels_top,
                    "1bp Insertion",
                    fontsize=35,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.155,
                    yText_labels_bottom_sec,
                    "Homopolymer Length",
                    fontsize=30,
                    fontname="Times New Roman",
                    weight="bold",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.505,
                    yText_labels_bottom_sec,
                    "Homopolymer Length",
                    fontsize=30,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.827,
                    yText_labels_top,
                    ">1bp",
                    fontsize=30,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.83,
                    yText_labels_bottom_sec,
                    "Type",
                    fontsize=30,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

                x = 0.127
                yText_labels_bottom = y_bottom - 0.025

                for l in range(0, 4, 1):
                    if l < 2:
                        for i in range(1, 6, 1):
                            plt.text(
                                x,
                                yText_labels_bottom,
                                str(i),
                                fontsize=25,
                                fontweight="bold",
                                fontname="Times New Roman",
                                color="black",
                                transform=plt.gcf().transFigure,
                            )
                            x += 0.028
                        x -= 0.005
                        plt.text(
                            x,
                            yText_labels_bottom,
                            "6+",
                            fontsize=25,
                            fontweight="bold",
                            fontname="Times New Roman",
                            color="black",
                            transform=plt.gcf().transFigure,
                        )
                        x += 0.037
                    else:
                        if l == 2:
                            x += 0
                        for i in range(0, 5, 1):
                            plt.text(
                                x,
                                yText_labels_bottom,
                                str(i),
                                fontsize=25,
                                fontweight="bold",
                                fontname="Times New Roman",
                                color="black",
                                transform=plt.gcf().transFigure,
                            )
                            x += 0.028
                        x -= 0.005
                        plt.text(
                            x,
                            yText_labels_bottom,
                            "5+",
                            fontsize=25,
                            fontweight="bold",
                            fontname="Times New Roman",
                            color="black",
                            transform=plt.gcf().transFigure,
                        )
                        x += 0.037

                yText_labels_bottom += 0.0
                plt.text(
                    x,
                    yText_labels_bottom,
                    "Del",
                    fontsize=17,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                    rotation="vertical",
                )
                x += 0.026
                plt.text(
                    x,
                    yText_labels_bottom,
                    "Ins",
                    fontsize=17,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                    rotation="vertical",
                )
                x += 0.0295
                yText_labels_bottom -= 0.001
                plt.text(
                    x,
                    yText_labels_bottom,
                    "MH",
                    fontsize=17,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                    rotation="vertical",
                )
                x += 0.0295
                yText_labels_bottom -= 0.002
                plt.text(
                    x,
                    yText_labels_bottom,
                    "COMP",
                    fontsize=10,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                    rotation="vertical",
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                if not percentage:
                    ylabels = getylabels(ylabels)

                panel1.set_xlim([0, 28])
                panel1.set_ylim([0, y])
                panel1.set_yticks(ylabs)

                if sig_probs:
                    plt.text(
                        0.13,
                        0.85,
                        sample,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.13,
                        0.85,
                        sample + ": " + "{:,}".format(int(total_count)) + " indels",
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                panel1.set_yticklabels(ylabels, fontsize=30)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                if percentage:
                    plt.ylabel(
                        "Percentage of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="gray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]

                pp.savefig(plot1)
                plt.close()
            pp.close()

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "ID_simple_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif (
        plot_type == "IDSB"
        or plot_type == "415"
    ):
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split("\t")
            mutation_type = first_line[0]
            mutation_type_list = mutation_type.split(":")
            if len(mutation_type_list) != 5:
                print(mutation_type_list)
                sys.exit(
                    "The matrix does not match the correct ID-96 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"ID_TSB_plots_{project}.pdf")
        pp = PdfPages(file_path)

        indel_types_tsb = []
        tsb_I = ["T", "U", "N", "B", "Q"]
        indel_types = [
            "1:Del:C:0",
            "1:Del:C:1",
            "1:Del:C:2",
            "1:Del:C:3",
            "1:Del:C:4",
            "1:Del:C:5",
            "1:Del:T:0",
            "1:Del:T:1",
            "1:Del:T:2",
            "1:Del:T:3",
            "1:Del:T:4",
            "1:Del:T:5",
            "1:Ins:C:0",
            "1:Ins:C:1",
            "1:Ins:C:2",
            "1:Ins:C:3",
            "1:Ins:C:4",
            "1:Ins:C:5",
            "1:Ins:T:0",
            "1:Ins:T:1",
            "1:Ins:T:2",
            "1:Ins:T:3",
            "1:Ins:T:4",
            "1:Ins:T:5",
            # >1bp INDELS
            "2:Del:R:0",
            "2:Del:R:1",
            "2:Del:R:2",
            "2:Del:R:3",
            "2:Del:R:4",
            "2:Del:R:5",
            "3:Del:R:0",
            "3:Del:R:1",
            "3:Del:R:2",
            "3:Del:R:3",
            "3:Del:R:4",
            "3:Del:R:5",
            "4:Del:R:0",
            "4:Del:R:1",
            "4:Del:R:2",
            "4:Del:R:3",
            "4:Del:R:4",
            "4:Del:R:5",
            "5:Del:R:0",
            "5:Del:R:1",
            "5:Del:R:2",
            "5:Del:R:3",
            "5:Del:R:4",
            "5:Del:R:5",
            "2:Ins:R:0",
            "2:Ins:R:1",
            "2:Ins:R:2",
            "2:Ins:R:3",
            "2:Ins:R:4",
            "2:Ins:R:5",
            "3:Ins:R:0",
            "3:Ins:R:1",
            "3:Ins:R:2",
            "3:Ins:R:3",
            "3:Ins:R:4",
            "3:Ins:R:5",
            "4:Ins:R:0",
            "4:Ins:R:1",
            "4:Ins:R:2",
            "4:Ins:R:3",
            "4:Ins:R:4",
            "4:Ins:R:5",
            "5:Ins:R:0",
            "5:Ins:R:1",
            "5:Ins:R:2",
            "5:Ins:R:3",
            "5:Ins:R:4",
            "5:Ins:R:5",
            # MicroHomology INDELS
            "2:Del:M:1",
            "3:Del:M:1",
            "3:Del:M:2",
            "4:Del:M:1",
            "4:Del:M:2",
            "4:Del:M:3",
            "5:Del:M:1",
            "5:Del:M:2",
            "5:Del:M:3",
            "5:Del:M:4",
            "5:Del:M:5",
        ]

        for indels in indel_types:
            for tsbs in tsb_I:
                indel_types_tsb.append(tsbs + ":" + indels)

        sig_probs = False
        mutations = OrderedDict()
        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                if pcawg:
                    samples = first_line.strip().split(",")
                    samples = samples[4:]
                    samples = [x.replace('"', "") for x in samples]
                else:
                    samples = first_line.strip().split("\t")
                    samples = samples[1:]
                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["1DelC"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["1DelT"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["1InsC"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["1InsT"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["2DelR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["3DelR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["4DelR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["5DelR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["2InsR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["3InsR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["3InsR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["4InsR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["5InsR"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]
                    mutations[sample]["2DelM"] = [[0, 0]]
                    mutations[sample]["3DelM"] = [[0, 0], [0, 0]]
                    mutations[sample]["4DelM"] = [[0, 0], [0, 0], [0, 0]]
                    mutations[sample]["5DelM"] = [
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                        [0, 0],
                    ]

                for lines in f:
                    if pcawg:
                        line = lines.strip().split(",")
                        line = [x.replace('"', "") for x in line]
                        if line[1] == "repeats":
                            mut_type = (
                                line[2][0]
                                + line[0][0]
                                + line[0][1].lower()
                                + line[0][2].lower()
                                + "R"
                            )
                        else:
                            mut_type = (
                                line[2][0]
                                + line[0][0]
                                + line[0][1].lower()
                                + line[0][2].lower()
                                + line[1][0]
                            )
                        try:
                            repeat_size = int(line[3])
                        except:
                            repeat_size = int(line[3][0])
                        if line[1] == "MH":
                            repeat_size -= 1
                        sample_index = 4
                    else:
                        line = lines.strip().split()
                        if line[0] not in indel_types_tsb:
                            continue
                        categories = line[0].split(":")
                        bias = categories[0]
                        if bias == "B" or bias == "N" or bias == "Q":
                            continue
                        mut_type = categories[1] + categories[2] + categories[3]

                        repeat_size = int(categories[4])
                        if categories[3] == "M":
                            repeat_size -= 1
                        sample_index = 1

                    for sample in samples:
                        if mut_type in mutations[sample].keys():
                            if percentage:
                                mutCount = float(line[sample_index])
                                if mutCount < 1 and mutCount > 0:
                                    sig_probs = True
                            else:
                                try:
                                    mutCount = int(line[sample_index])
                                except:
                                    print(
                                        "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                        end="",
                                    )

                                # mutCount = int(line[sample_index])
                            if bias == "T":
                                mutations[sample][mut_type][repeat_size][0] = mutCount
                            else:
                                mutations[sample][mut_type][repeat_size][1] = mutCount
                        else:
                            continue
                        sample_index += 1

            sample_count = 0
            for sample in mutations.keys():
                total_count = sum(
                    sum(sum(tsb) for tsb in nuc) for nuc in mutations[sample].values()
                )
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(43.93, 12))
                plt.rc("axes", edgecolor="black")
                panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
                xlabels = []

                x = 0.4
                ymax = 0
                colors = [
                    [253 / 256, 190 / 256, 111 / 256],
                    [255 / 256, 128 / 256, 2 / 256],
                    [176 / 256, 221 / 256, 139 / 256],
                    [54 / 256, 161 / 256, 46 / 256],
                    [253 / 256, 202 / 256, 181 / 256],
                    [252 / 256, 138 / 256, 106 / 256],
                    [241 / 256, 68 / 256, 50 / 256],
                    [188 / 256, 25 / 256, 26 / 256],
                    [208 / 256, 225 / 256, 242 / 256],
                    [148 / 256, 196 / 256, 223 / 256],
                    [74 / 256, 152 / 256, 201 / 256],
                    [23 / 256, 100 / 256, 171 / 256],
                    [226 / 256, 226 / 256, 239 / 256],
                    [182 / 256, 182 / 256, 216 / 256],
                    [134 / 256, 131 / 256, 189 / 256],
                    [98 / 256, 64 / 256, 155 / 256],
                ]

                i = 0
                for key in mutations[sample]:
                    l = 1
                    for seq in mutations[sample][key]:
                        xlabels.append(l)
                        if percentage:
                            if total_count > 0:
                                trans = plt.bar(
                                    x,
                                    seq[0] / total_count * 100,
                                    width=0.2,
                                    color=[1 / 256, 70 / 256, 102 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-transcribed Strand",
                                )
                                x += 0.2
                                untrans = plt.bar(
                                    x,
                                    seq[1] / total_count * 100,
                                    width=0.2,
                                    color=[228 / 256, 41 / 256, 38 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-untranscribed Strand",
                                )
                                if seq[0] / total_count * 100 > ymax:
                                    ymax = seq[0] / total_count * 100
                                if seq[1] / total_count * 100 > ymax:
                                    ymax = seq[1] / total_count * 100

                        else:
                            trans = plt.bar(
                                x,
                                seq[0],
                                width=0.2,
                                color=[1 / 256, 70 / 256, 102 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-transcribed Strand",
                            )
                            x += 0.2
                            untrans = plt.bar(
                                x,
                                seq[1],
                                width=0.2,
                                color=[228 / 256, 41 / 256, 38 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-untranscribed Strand",
                            )
                            if seq[0] > ymax:
                                ymax = seq[0]
                            if seq[1] > ymax:
                                ymax = seq[1]

                        x += 0.799
                        l += 1
                    i += 1

                x = 0.0475
                y_top = 0.827
                y_bottom = 0.114
                y = int(ymax * 1.25)
                y2 = y + 2
                for i in range(0, 12, 1):
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y_top),
                            0.0595,
                            0.05,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    panel1.add_patch(
                        plt.Rectangle(
                            (x, y_bottom),
                            0.0595,
                            0.05,
                            facecolor=colors[i],
                            clip_on=False,
                            transform=plt.gcf().transFigure,
                        )
                    )
                    x += 0.0665

                panel1.add_patch(
                    plt.Rectangle(
                        (x - 0.001, y_top),
                        0.006,
                        0.05,
                        facecolor=colors[12],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (x - 0.001, y_bottom),
                        0.006,
                        0.05,
                        facecolor=colors[12],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                x += 0.011
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_top),
                        0.0155,
                        0.05,
                        facecolor=colors[13],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_bottom),
                        0.0155,
                        0.05,
                        facecolor=colors[13],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                x += 0.022
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_top),
                        0.027,
                        0.05,
                        facecolor=colors[14],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_bottom),
                        0.027,
                        0.05,
                        facecolor=colors[14],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                x += 0.0335
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_top),
                        0.049,
                        0.05,
                        facecolor=colors[15],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (x, y_bottom),
                        0.049,
                        0.05,
                        facecolor=colors[15],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )

                yText = y_top + 0.01
                plt.text(
                    0.072,
                    yText,
                    "C",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.1385,
                    yText,
                    "T",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.205,
                    yText,
                    "C",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.2715,
                    yText,
                    "T",
                    fontsize=40,
                    fontname="Times New Roman",
                    fontweight="bold",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.338,
                    yText,
                    "2",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.4045,
                    yText,
                    "3",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.471,
                    yText,
                    "4",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.5375,
                    yText,
                    "5+",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.604,
                    yText,
                    "2",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.6705,
                    yText,
                    "3",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.737,
                    yText,
                    "4",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.8035,
                    yText,
                    "5+",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.844,
                    yText,
                    "2",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.861,
                    yText,
                    "3",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.888,
                    yText,
                    "4",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.93,
                    yText,
                    "5+",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )

                yText_labels_top = yText + 0.075
                yText_labels_bottom = y_bottom - 0.03
                yText_labels_bottom_sec = yText_labels_bottom - 0.045

                plt.text(
                    0.08,
                    yText_labels_top,
                    "1bp Deletion",
                    fontsize=40,
                    fontname="Times New Roman",
                    weight="bold",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.21,
                    yText_labels_top,
                    "1bp Insertion",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.375,
                    yText_labels_top,
                    ">1bp Deletion at Repeats\n      (Deletion Length)",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.64,
                    yText_labels_top,
                    ">1bp Insertion at Repeats\n       (Insertion Length)",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.85,
                    yText_labels_top,
                    " Microhomology\n(Deletion Length)",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

                plt.text(
                    0.058,
                    yText_labels_bottom_sec,
                    "Homopolymer Length",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.19,
                    yText_labels_bottom_sec,
                    "Homopolymer Length",
                    fontsize=35,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.39,
                    yText_labels_bottom_sec,
                    "Number of Repeat Units",
                    fontsize=35,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.65,
                    yText_labels_bottom_sec,
                    "Number of Repeat Units",
                    fontsize=35,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.85,
                    yText_labels_bottom_sec,
                    "Microhomology Length",
                    fontsize=35,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

                x = 0.0477
                for i in range(0, 8, 1):
                    if i != 2 and i != 3:
                        plt.text(
                            x,
                            yText_labels_bottom,
                            "1  2  3  4  5  6+",
                            fontsize=32,
                            fontweight="bold",
                            fontname="Times New Roman",
                            color="black",
                            transform=plt.gcf().transFigure,
                        )
                    else:
                        plt.text(
                            x,
                            yText_labels_bottom,
                            "0  1  2  3  4  5+",
                            fontsize=32,
                            fontweight="bold",
                            fontname="Times New Roman",
                            color="black",
                            transform=plt.gcf().transFigure,
                        )

                    x += 0.0665

                for i in range(0, 4, 1):
                    plt.text(
                        x,
                        yText_labels_bottom,
                        "0  1  2  3  4  5+",
                        fontsize=32,
                        fontweight="bold",
                        fontname="Times New Roman",
                        color="black",
                        transform=plt.gcf().transFigure,
                    )
                    x += 0.0665

                plt.text(
                    x,
                    yText_labels_bottom,
                    "1",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                x += 0.011
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                x += 0.022
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
                x += 0.0335
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5+",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.375, 83.375, 1)

                if not percentage:
                    ylabels = getylabels(ylabels)

                panel1.set_xlim([0, 83])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)

                if sig_probs:
                    plt.text(
                        0.0475,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.0475,
                        0.75,
                        sample + ": " + "{:,}".format(int(total_count)) + " indels",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.95
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.72,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.68,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=30)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                panel1.legend(handles=[trans, untrans], prop={"size": 30})

                if percentage:
                    plt.ylabel(
                        "Percentage of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="gray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]

                pp.savefig(plot1)
                plt.close()
                sample_count += 1
            pp.close()

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "ID_TSB_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    else:
        print(
            "Error: The function plotID does not support plot_type",
            plot_type,
            "so no plot has been generated."
        )


def plotDBS(
    matrix_path,
    output_path,
    project,
    plot_type,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,
):
    _require_pdf_for_legacy_context(
        savefig_format, plot_type, {"78", "78DBS", "DBS78"}, "plotDBS"
    )
    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    # load custom fonts for plotting
    load_custom_fonts()

    plot_custom_text = False
    pcawg = False
    sig_probs = False
    if plot_type == "78" or plot_type == "78DBS" or plot_type == "DBS78":
        data = process_input(matrix_path, plot_type)

        dinucs = [
            "TT>GG",
            "TT>CG",
            "TT>AG",
            "TT>GC",
            "TT>CC",
            "TT>AC",
            "TT>GA",
            "TT>CA",
            "TT>AA",
            "AC>CA",
            "AC>CG",
            "AC>CT",
            "AC>GA",
            "AC>GG",
            "AC>GT",
            "AC>TA",
            "AC>TG",
            "AC>TT",
            "CT>AA",
            "CT>AC",
            "CT>AG",
            "CT>GA",
            "CT>GC",
            "CT>GG",
            "CT>TG",
            "CT>TC",
            "CT>TA",
            "AT>CA",
            "AT>CC",
            "AT>CG",
            "AT>GA",
            "AT>GC",
            "AT>TA",
            "TG>GT",
            "TG>CT",
            "TG>AT",
            "TG>GC",
            "TG>CC",
            "TG>AC",
            "TG>GA",
            "TG>CA",
            "TG>AA",
            "CC>AA",
            "CC>AG",
            "CC>AT",
            "CC>GA",
            "CC>GG",
            "CC>GT",
            "CC>TA",
            "CC>TG",
            "CC>TT",
            "CG>AT",
            "CG>GC",
            "CG>GT",
            "CG>TC",
            "CG>TA",
            "CG>TT",
            "TC>GT",
            "TC>CT",
            "TC>AT",
            "TC>GG",
            "TC>CG",
            "TC>AG",
            "TC>GA",
            "TC>CA",
            "TC>AA",
            "GC>AA",
            "GC>AG",
            "GC>AT",
            "GC>CA",
            "GC>CG",
            "GC>TA",
            "TA>GT",
            "TA>CT",
            "TA>AT",
            "TA>GG",
            "TA>CG",
            "TA>GC",
        ]

        revcompl = lambda x: "".join(
            [{"A": "T", "C": "G", "G": "C", "T": "A"}[B] for B in x][::-1]
        )
        mutations = OrderedDict()

        try:
            vals = set(list(data.index)) - set(dinucs)
            vals = sorted(list(vals))
            for ech in vals:
                ech_mod = ech.split(">")[0] + ">" + revcompl(ech.split(">")[1])
                data.rename(index={ech: ech_mod}, inplace=True)

            data = data.sort_index()
            ctx = data.index
            xlabels = [dn.split(">")[1] for dn in ctx]
            colors = [
                [3 / 256, 189 / 256, 239 / 256],
                [3 / 256, 102 / 256, 204 / 256],
                [162 / 256, 207 / 256, 99 / 256],
                [1 / 256, 102 / 256, 1 / 256],
                [255 / 256, 153 / 256, 153 / 256],
                [228 / 256, 41 / 256, 38 / 256],
                [255 / 256, 178 / 256, 102 / 256],
                [255 / 256, 128 / 256, 1 / 256],
                [204 / 256, 153 / 256, 255 / 256],
                [76 / 256, 1 / 256, 153 / 256],
            ]
            mainlist = [dn.split(">")[0] for dn in ctx]
            le = LabelEncoder()
            colors_idxs = le.fit_transform(mainlist)
            colors_flat_list = [colors[i] for i in colors_idxs]
            sample_count = 0

            buf = io.BytesIO()
            fig_orig = make_pickle_file(
                context="DBS78", return_plot_template=True, volume=volume
            )
            pickle.dump(fig_orig, buf)

            figs = {}

            for sample in data.columns:
                buf.seek(0)
                figs[sample] = pickle.load(buf)
                panel1 = figs[sample].axes[0]
                total_count = np.sum(
                    data[sample].values
                )  # sum(sum(nuc.values()) for nuc in mutations[sample].values())

                x = 0.4
                muts = data[sample].values
                ymax = 0
                if percentage:
                    if total_count > 0:
                        plt.bar(
                            np.asarray(range(len(ctx))) + x,
                            muts / total_count * 100,
                            width=0.4,
                            color=colors_flat_list,
                            align="center",
                            zorder=1000,
                        )
                        ymax = np.max(muts / total_count * 100)
                    sig_probs = True
                else:
                    plt.bar(
                        np.asarray(range(len(ctx))) + x,
                        muts,
                        width=0.4,
                        color=colors_flat_list,
                        align="center",
                        zorder=1000,
                    )
                    ymax = np.max(muts)
                # for i in range(len(xlabels)):
                #     print(xlabels[i],muts[i])

                x = 0.043
                y3 = 0.87
                y = int(ymax * 1.25)
                y2 = y + 2
                i = 0

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count))
                        + " double subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_bottom:
                    if len(custom_text_bottom[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.98
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.75,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.7,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                if not percentage:
                    ylabels = spplt.getylabels(ylabels)

                labs = np.arange(0.44, 78.44, 1)
                panel1.set_xlim([0, 78])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                panel1.set_xticklabels(
                    xlabels,
                    rotation="vertical",
                    fontsize=30,
                    color="grey",
                    fontname="Courier New",
                    verticalalignment="top",
                    fontweight="bold",
                )

                panel1.set_yticklabels(ylabels, fontsize=25)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(
                    which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
                )
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                if percentage:
                    plt.ylabel(
                        "Percentage of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=True,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="lightgray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                [i.set_color("grey") for i in plt.gca().get_xticklabels()]
                sample_count += 1

            return output_results(
                savefig_format, output_path, project, figs, "DBS_78", dpi=dpi
            )

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "DBS_78_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    elif (
        plot_type == "312"
        or plot_type == "78SB"
        or plot_type == "SB78"
        or plot_type == "186"
    ):
        with open(matrix_path) as f:
            next(f)
            first_line = f.readline()
            first_line = first_line.strip().split()
            mutation_type = first_line[0]
            if len(mutation_type) != 7 and mutation_type[1] != ":":
                sys.exit(
                    "The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function."
                )

        file_path = os.path.join(output_path, f"DBS_186_plots_{project}.pdf")
        pp = PdfPages(file_path)

        dinucs = [
            "TT>GG",
            "TT>CG",
            "TT>AG",
            "TT>GC",
            "TT>CC",
            "TT>AC",
            "TT>GA",
            "TT>CA",
            "TT>AA",
            "CT>AA",
            "CT>AC",
            "CT>AG",
            "CT>GA",
            "CT>GC",
            "CT>GG",
            "CT>TG",
            "CT>TC",
            "CT>TA",
            "CC>AA",
            "CC>AG",
            "CC>AT",
            "CC>GA",
            "CC>GG",
            "CC>GT",
            "CC>TA",
            "CC>TG",
            "CC>TT",
            "TC>GT",
            "TC>CT",
            "TC>AT",
            "TC>GG",
            "TC>CG",
            "TC>AG",
            "TC>GA",
            "TC>CA",
            "TC>AA",
        ]

        revcompl = lambda x: "".join(
            [{"A": "T", "C": "G", "G": "C", "T": "A", ">": ">"}[B] for B in x][::-1]
        )
        mutations = OrderedDict()

        try:
            with open(matrix_path) as f:
                first_line = f.readline()
                samples = first_line.strip().split("\t")
                samples = samples[1:]
                for sample in samples:
                    mutations[sample] = OrderedDict()
                    mutations[sample]["CC"] = OrderedDict()
                    mutations[sample]["CT"] = OrderedDict()
                    mutations[sample]["TC"] = OrderedDict()
                    mutations[sample]["TT"] = OrderedDict()

                for lines in f:
                    line = lines.strip().split()
                    mut = line[0][2:]
                    nuc = line[0][5:]
                    mut_type = line[0][2:4]
                    bias = line[0][0]
                    if bias == "N" or bias == "B" or bias == "Q":
                        continue
                    else:
                        if mut not in dinucs:
                            if revcompl(mut) not in dinucs:
                                continue
                            nuc = revcompl(nuc)
                            mut_type = revcompl(mut_type)
                        sample_index = 1

                        for sample in samples:
                            if percentage:
                                mutCount = float(line[sample_index])
                                if mutCount < 1 and mutCount > 0:
                                    sig_probs = True
                            else:
                                # mutCount = int(line[sample_index])
                                try:
                                    mutCount = int(line[sample_index])
                                except:
                                    print(
                                        "It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ",
                                        end="",
                                    )

                            if nuc not in mutations[sample][mut_type]:
                                mutations[sample][mut_type][nuc] = [0, 0]
                            if bias == "T":
                                mutations[sample][mut_type][nuc][0] = mutCount
                            else:
                                mutations[sample][mut_type][nuc][1] = mutCount
                            sample_index += 1

            for sample in mutations.keys():
                total_count = sum(
                    sum(sum(tsb) for tsb in nuc.values())
                    for nuc in mutations[sample].values()
                )
                plt.rcParams["axes.linewidth"] = 2
                plot1 = plt.figure(figsize=(21, 9.92))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.07, 0.09, 0.92, 0.77])
                xlabels = []

                x = 0.3
                ymax = 0
                i = 0
                colors = [
                    [3 / 256, 189 / 256, 239 / 256],
                    [3 / 256, 102 / 256, 204 / 256],
                    [162 / 256, 207 / 256, 99 / 256],
                    [1 / 256, 102 / 256, 1 / 256],
                    [255 / 256, 153 / 256, 153 / 256],
                    [228 / 256, 41 / 256, 38 / 256],
                    [255 / 256, 178 / 256, 102 / 256],
                    [255 / 256, 128 / 256, 1 / 256],
                    [204 / 256, 153 / 256, 255 / 256],
                    [76 / 256, 1 / 256, 153 / 256],
                ]
                for key in mutations[sample]:
                    muts = mutations[sample][key].keys()
                    muts = sorted(muts)
                    for seq in muts:
                        xlabels.append(seq)
                        if percentage:
                            try:
                                trans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][0] / total_count * 100,
                                    width=0.2,
                                    color=[1 / 256, 70 / 256, 102 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-transcribed Strand",
                                )
                                x += 0.2
                                untrans = plt.bar(
                                    x,
                                    mutations[sample][key][seq][1] / total_count * 100,
                                    width=0.2,
                                    color=[228 / 256, 41 / 256, 38 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-untranscribed Strand",
                                )
                                x += 0.8
                                if (
                                    mutations[sample][key][seq][0] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][0]
                                        / total_count
                                        * 100
                                    )
                                if (
                                    mutations[sample][key][seq][1] / total_count * 100
                                    > ymax
                                ):
                                    ymax = (
                                        mutations[sample][key][seq][1]
                                        / total_count
                                        * 100
                                    )
                            except:
                                trans = plt.bar(
                                    x,
                                    0,
                                    width=0.2,
                                    color=[1 / 256, 70 / 256, 102 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-transcribed Strand",
                                )
                                untrans = plt.bar(
                                    x,
                                    0,
                                    width=0.2,
                                    color=[228 / 256, 41 / 256, 38 / 256],
                                    align="center",
                                    zorder=1000,
                                    label="Genic-untranscribed Strand",
                                )

                        else:
                            trans = plt.bar(
                                x,
                                mutations[sample][key][seq][0],
                                width=0.2,
                                color=[1 / 256, 70 / 256, 102 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-transcribed Strand",
                            )
                            x += 0.2
                            untrans = plt.bar(
                                x,
                                mutations[sample][key][seq][1],
                                width=0.2,
                                color=[228 / 256, 41 / 256, 38 / 256],
                                align="center",
                                zorder=1000,
                                label="Genic-untranscribed Strand",
                            )
                            x += 0.8
                            if mutations[sample][key][seq][0] > ymax:
                                ymax = mutations[sample][key][seq][0]
                            if mutations[sample][key][seq][1] > ymax:
                                ymax = mutations[sample][key][seq][1]
                    i += 1

                y3 = 0.87
                y = int(ymax * 1.25)

                panel1.add_patch(
                    plt.Rectangle(
                        (0.075, y3),
                        0.218,
                        0.05,
                        facecolor=colors[0],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (0.302, y3),
                        0.218,
                        0.05,
                        facecolor=colors[2],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (0.532, y3),
                        0.218,
                        0.05,
                        facecolor=colors[4],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )
                panel1.add_patch(
                    plt.Rectangle(
                        (0.765, y3),
                        0.218,
                        0.05,
                        facecolor=colors[7],
                        clip_on=False,
                        transform=plt.gcf().transFigure,
                    )
                )

                yText = y3 + 0.06
                plt.text(
                    0.13,
                    yText,
                    "CC>NN",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.37,
                    yText,
                    "CT>NN",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.59,
                    yText,
                    "TC>NN",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
                plt.text(
                    0.83,
                    yText,
                    "TT>NN",
                    fontsize=40,
                    fontweight="bold",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                x_shaded = 0
                panel1.add_patch(
                    plt.Rectangle(
                        (x_shaded, 0),
                        8.9,
                        y,
                        facecolor=colors[0],
                        zorder=0,
                        alpha=0.25,
                        edgecolor="grey",
                    )
                )
                x_shaded += 8.9
                panel1.add_patch(
                    plt.Rectangle(
                        (x_shaded, 0),
                        9,
                        y,
                        facecolor=colors[2],
                        zorder=0,
                        alpha=0.25,
                        edgecolor="grey",
                    )
                )
                x_shaded += 9
                panel1.add_patch(
                    plt.Rectangle(
                        (x_shaded, 0),
                        9,
                        y,
                        facecolor=colors[4],
                        zorder=0,
                        alpha=0.25,
                        edgecolor="grey",
                    )
                )
                x_shaded += 9
                panel1.add_patch(
                    plt.Rectangle(
                        (x_shaded, 0),
                        9.1,
                        y,
                        facecolor=colors[7],
                        zorder=0,
                        alpha=0.25,
                        edgecolor="grey",
                    )
                )

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    if ytick_offest == 0:
                        ytick_offest = 1
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                if sig_probs:
                    plt.text(
                        0.08,
                        0.8,
                        sample,
                        fontsize=35,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.08,
                        0.8,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count))
                        + " transcribed double subs",
                        fontsize=35,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                if not percentage:
                    ylabels = getylabels(ylabels)

                labs = np.arange(0.55, 36.44, 1)
                panel1.set_xlim([0, 36])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                panel1.set_xticklabels(
                    xlabels,
                    rotation="vertical",
                    fontsize=30,
                    color="grey",
                    fontname="Courier New",
                    verticalalignment="top",
                    fontweight="bold",
                )

                panel1.set_yticklabels(ylabels, fontsize=25)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(
                    which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
                )
                panel1.set_xlabel("")
                panel1.set_ylabel("")
                panel1.legend(handles=[trans, untrans], prop={"size": 30})

                if percentage:
                    plt.ylabel(
                        "Percentage of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=True,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="lightgray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                [i.set_color("grey") for i in plt.gca().get_xticklabels()]

                panel1.set_xlim([0, 36])
                pp.savefig(plot1)
                plt.close()
            pp.close()
        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "DBS_186_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    else:
        print(
            "Error: The function plotDBS does not support plot_type",
            plot_type,
            "so no plot has been generated."
        )

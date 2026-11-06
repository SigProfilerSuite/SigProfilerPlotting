# Plotting SBS4608

SigProfilerPlotting can render a strand-tagged SBS4608 matrix as either an
interactive HTML report or a static, multi-page PDF. Both formats combine a
96-context profile, pentanucleotide heatmaps, and strand/genic bias summaries.

!!! note
    `plotSBS(..., plot_type="4608", ...)` is the single static entry point
    for this format -- there is no separate `plotStatic()` function. Use
    `plotInteractive()` for the interactive HTML report.

## Input format

The input must be a tab-separated matrix with:

- a first column named exactly `MutationType`;
- one or more sample or signature columns containing finite, non-negative
  numeric values; and
- all 4,608 unique mutation rows: 6 substitution classes × 16 left
  dinucleotides × 16 right dinucleotides × 3 strand categories.

Mutation labels use the form `strand:NN[N>N]NN`, for example
`T:AA[C>A]AA`. The supported strand prefixes are:

| Prefix | Meaning |
| --- | --- |
| `T` | Transcribed strand |
| `U` | Untranscribed strand |
| `N` | Non-transcribed/intergenic region |

The beginning of a valid file looks like this:

```text
MutationType	Sample_1	Sample_2
T:AA[C>A]AA	12	0.0012
U:AA[C>A]AA	8	0.0008
N:AA[C>A]AA	4	0.0004
...
```

This is a signature/profile file -- rows are contributions, not raw sample
mutation counts (SigProfilerExtractor's own `*_Signatures.txt` output format).
`plotSBS(..., "4608", ...)` is file-path only (unlike every other `plotSBS`
context, which also accepts a DataFrame), since it reads this format
directly, the same way `plotInteractive()` does.

## Interactive HTML report

Use `plotInteractive()` to create a self-contained HTML report:

```python
import sigProfilerPlotting as sigPlt

html_path = sigPlt.plotInteractive(
    matrix_path="results/Project.SBS4608.all",
    project="Project",
    output_path="plots",
    signatures=["Sample_1", "Sample_2"],
    include_detailed_heatmaps=True,
    main_plot_type="default",
)

print(html_path)
```

Plotly is embedded in the generated file, so the report can be opened locally
and shared without an internet connection. The output in this example is
`plots/Project_interactive.html`.

The report provides these controls:

- **Main Plot Type** switches the upper panel among the default SBS96 profile,
  Transcribed/Untranscribed strand bars, and Genic/Intergenic bars.
- **Transcription bias** shows or hides the lower strand and genic summary
  bars.
- **5-base context (total heatmap)** shows or hides the combined
  pentanucleotide heatmap.
- **Detailed heatmaps (front/back)** shows or hides the two detailed
  pentanucleotide heatmaps. This control is available only when
  `include_detailed_heatmaps=True`.

When an option is cleared, its panel is removed and the remaining panels are
repacked automatically without leaving an empty section.

## Static PDF report

Use `plotSBS(..., "4608", ...)` when a PDF is more suitable for a manuscript,
presentation, or print workflow:

```python
import sigProfilerPlotting as sigPlt

pdf_path = sigPlt.plotSBS(
    "results/Project.SBS4608.all",  # matrix_path
    "plots",                        # output_path
    "Project",                      # project
    "4608",                         # plot_type
    signatures=["Sample_1", "Sample_2"],
    custom_text_upper=["Stability: 0.97", "Stability: 0.88"],
    custom_text_middle=["Mutations: 12,345", "Mutations: 9,021"],
)

print(pdf_path)
```

The output in this example is `plots/SBS_4608_plots_Project.pdf`, with one
page per selected sample or signature -- following the same
`<TYPE>_<context>_plots_<project>.pdf` naming every other `plotSBS`/`plotID`
context uses. The static report always shows the default SBS96-profile view
with front/back detailed heatmaps included (a PDF has no live controls to
switch views with, unlike the interactive report, so there is no
`main_plot_type`/`include_detailed_heatmaps` choice to make here).

`percentage` defaults to `True` for this context specifically -- the input
is already a signature/profile file, not raw mutation counts, matching how
the interactive report and the (now removed) `plotStatic()` always treated
it. Pass `percentage=False` to instead show each panel's raw contribution
values from the matrix (rescaled from the builder's native percent-of-total
representation, not a second independent calculation).

`custom_text_upper`/`custom_text_middle`/`custom_text_bottom` work the same
way they do for every other `plotSBS` context -- one string per selected
signature (in order), up to 40 characters each, drawn in the page's
top-right corner. This is the mechanism SigProfilerExtractor uses to attach
a signature's stability score and mutation count to its plots.

Static export uses Plotly and Kaleido and requires Chrome or Chromium. If
Kaleido cannot locate a browser, install one with:

```console
plotly_get_chrome
```

## Parameters

| Parameter | Applies to | Description |
| --- | --- | --- |
| `matrix_path` | both | Path to the tab-separated, strand-tagged SBS4608 matrix (file path only for the static report). |
| `output_path` | both | Output directory. |
| `project` | both | Project label used in the report title and output filename. |
| `signatures` | both | List of columns to plot. Defaults to every sample/signature column. |
| `percentage` | static only | `True` (default) shows percent-of-total; `False` shows raw contribution values. |
| `custom_text_upper`, `custom_text_middle`, `custom_text_bottom` | static only | Per-signature annotation text (e.g. stability, mutation count), one string each, <=40 characters, in selection order. |
| `include_detailed_heatmaps` | interactive only | Include front/back heatmaps in addition to the total heatmap. Defaults to `True`. Fixed on for the static report. |
| `main_plot_type` | interactive only | Initial view: `"default"`, `"strand_bias"`, or `"genic_bias"`. The static report always uses `"default"`. |

The static report returns the absolute path of the generated PDF; the
interactive report returns the absolute path of the generated HTML file.

## Common input errors

- **No sample columns are found:** confirm that the first header is exactly
  `MutationType` and that the file is tab-separated.
- **Incomplete matrix:** provide all 4,608 required `T`, `U`, and `N` rows,
  including rows whose value is zero.
- **Invalid mutation label:** use a label such as `T:AA[C>A]AA`; plain SBS96 or
  SBS1536 labels are not accepted by these functions.
- **Requested signature not found:** entries in `signatures` must match column
  headers exactly.
- **Static export fails:** install Chrome/Chromium or run `plotly_get_chrome`.
- **Unsupported `savefig_format`:** the static SBS4608 report only supports
  `savefig_format="pdf"` (the default); any other value raises a `ValueError`.

# Plotting SBS4608

SigProfilerPlotting can render a strand-tagged SBS4608 matrix as either an
interactive HTML report or a static, multi-page PDF. Both formats combine a
96-context profile, pentanucleotide heatmaps, and strand/genic bias summaries.

!!! note
    These functions expect the strand-tagged SBS4608 format described below.
    They are separate from `plotSBS(..., plot_type="4608")`, which produces the
    traditional static SBS bar-profile plot.

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

Rows may contain mutation counts or signature contributions. Values are
normalized to percentages for plotting.

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

Use `plotStatic()` when a PDF is more suitable for a manuscript, presentation,
or print workflow:

```python
import sigProfilerPlotting as sigPlt

pdf_path = sigPlt.plotStatic(
    matrix_path="results/Project.SBS4608.all",
    project="Project",
    output_path="plots",
    signatures=["Sample_1", "Sample_2"],
    include_detailed_heatmaps=True,
    main_plot_type="strand_bias",
)

print(pdf_path)
```

The output in this example is `plots/Project_static.pdf`, with one page per
selected sample or signature. A PDF cannot provide live controls, so
`main_plot_type` selects the single upper-panel view rendered on every page.

Static export uses Plotly and Kaleido and requires Chrome or Chromium. If
Kaleido cannot locate a browser, install one with:

```console
plotly_get_chrome
```

## Parameters

Both functions accept the same parameters:

| Parameter | Description |
| --- | --- |
| `matrix_path` | Path to the tab-separated, strand-tagged SBS4608 matrix. |
| `project` | Project label used in the report title and output filename. |
| `output_path` | Output directory. The default is `<project>_plots`. |
| `signatures` | List of columns to plot. The default plots every sample/signature column. |
| `include_detailed_heatmaps` | Include front/back heatmaps in addition to the total heatmap. Defaults to `True`. |
| `main_plot_type` | Initial interactive view or fixed static view: `"default"`, `"strand_bias"`, or `"genic_bias"`. |

Both functions return the absolute path of the generated report.

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

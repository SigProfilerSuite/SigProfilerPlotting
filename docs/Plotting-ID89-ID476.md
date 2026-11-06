# Plotting ID89 and ID476

SigProfilerPlotting can render the ID89 and ID476 indel classification
schemes described by [Koh et al. (2025)](https://doi.org/10.1038/s41588-025-02152-y)
and explored in the [Indel Signature Browser](https://indelsig.net/). ID89
groups small indels into 89 mutation categories; ID476 captures more detailed
sequence/repeat features across 476 categories. **Both are mutation
categories, not 89 or 476 extracted signatures** -- do not confuse them with
`plotSBS`/`plotID`-style signature plots.

Both contexts are supported by the static `plotID()` function (PDF, PNG, and
PIL image output) and by the [interactive dashboard](Plotting-Dashboard.md),
alongside the existing ID83 scheme. Requesting ID83 is unaffected either way.

## Generating the matrices

ID89 and ID476 matrices are produced by
[SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator)'s
own `indel_contexts=("89", "476")` option (see its own documentation for that
feature), not by SigProfilerPlotting. SigProfilerPlotting only plots the
resulting `project.ID89.all` / `project.ID476.all` matrices.

## Input format

Each matrix is a tab-separated file (or an equivalent `MutationType`-indexed
`pandas.DataFrame`) with:

- a first column named `MutationType` (or `index`/`Mutation Types`/
  `classification`, matching every other SigProfilerPlotting context);
- one or more sample/signature columns of finite, non-negative numeric
  values; and
- exactly 89 (ID89) or 476 (ID476) unique mutation-type rows, matching the
  exact label strings SigProfilerMatrixGenerator's Koh et al. classifier
  produces -- for example `[Del(C):R1]A` (ID89) or `A[Del(C):R1]A` (ID476).
  Rows may be supplied in any order; `plotID()` and the dashboard both
  canonicalize them, the same way they do for ID83.

Missing, duplicate, unexpected, non-numeric, negative, or non-finite rows are
all rejected with a message naming the exact row or column at fault, exactly
as for every other context.

## Static plots

```python
import sigProfilerPlotting as sigPlt

sigPlt.plotID(
    "project.ID89.all",   # or a DataFrame, or "project.ID476.all" with plot_type="476"
    "output/",
    "project",
    "89",                  # or "ID89" / "476" / "ID476"
    savefig_format="png",  # "pdf" and "pil_image" are also supported
)
```

Both contexts reproduce the grouping, header colors, and category order of
the updated `mSigPlot::plot_ID89()`/`plot_ID476()` reference implementations:

- **ID89** -- 12 groups: `Del 1 C`, `Del 1 T (1-4)`, `Del 1 T (5-7)`,
  `Del 1 T (8+)`, `Ins 1 C`, `Ins 1 T (0-4)`, `Ins 1 T (5-7)`,
  `Ins 1 T (8+)`, `Del >=2`, `Ins >=2`, `Del Mh`, and `Complex` (delins
  events, a category ID83 has no equivalent for).
- **ID476** -- 8 groups: `Del 1 C`, `Del 1 T`, `Ins 1 C`, `Ins 1 T`,
  `Del >=2`, `Ins >=2`, `Mh`, and `Complex`.

Each group's colored header band is sized proportionally to how many of the
89/476 rows belong to it (unlike ID83's fixed 16-category layout), and shows
that group's total count for the sample being plotted. `percentage=True`,
multiple samples/signatures, `custom_text_upper`/`custom_text_middle`/
`custom_text_bottom`, caller-selected `dpi`, and an all-zero sample all work
the same way they do for every other `plotID()` context.

## Interactive dashboard

Uploading a `project.ID89.all` or `project.ID476.all` matrix is detected
automatically, exactly like any other context -- no extra flag is needed.
Each gets its own Plotly renderer using the same category grouping, order,
and colors as the static plots above, with full mutation-type labels and
values on hover. Counts, percentage, and proportion units all work,
including an all-zero matrix.

## Attribution

The ID89 and ID476 classification schemes, category definitions, and
reference channel order are from Koh et al., *A redefined InDel taxonomy
provides insights into mutational signatures*, Nature Genetics 57, 1132-1141
(2025) ([10.1038/s41588-025-02152-y](https://doi.org/10.1038/s41588-025-02152-y)),
as implemented in SigProfilerMatrixGenerator and cross-validated there against
[`mSigSpectra`](https://github.com/steverozen/mSigSpectra) (GPL-3) and the
[Indel Signature Browser](https://indelsig.net/). SigProfilerPlotting's own
`reference_formats/ID89.txt` and `reference_formats/ID476.txt` channel-order
files are copied directly from SigProfilerMatrixGenerator's independently
checked fixtures so that matrices round-trip between the two packages without
any relabeling.

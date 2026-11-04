# Interactive Dashboard

SigProfilerPlotting includes a local, browser-based dashboard for exploring
any supported matrix interactively -- upload a file, pick samples or
signatures, and get real Plotly charts with hover, zoom, pan, and image
export, without writing any Python.

!!! note
    The dashboard runs entirely on your machine. It binds only to
    `127.0.0.1` and never sends your data anywhere else. The only way to
    load a matrix is through your browser's own file picker -- there is no
    server-side file browser and no network upload.

## Launching the dashboard

There are two ways to run the dashboard, depending on how you got the code.

### From a cloned repository

Open a terminal in the repository root (the folder containing
`launch_dashboard.sh`) and run the launcher for your platform:

| Platform | How to launch |
| --- | --- |
| macOS | `./launch_dashboard.sh`, or double-click `launch_dashboard.command` in Finder |
| Linux | `./launch_dashboard.sh` |
| Windows | Double-click `launch_dashboard.bat`, or run `launch_dashboard.bat` from a terminal/PowerShell |

The first time you run it, the launcher creates an isolated `.venv-dashboard`
virtual environment next to it and installs the package into it. This can
take a minute or two. Every later launch reuses that environment and starts
in a couple of seconds -- you'll see `Dependencies already installed;
skipping install.` printed instead of an install running.

### From a `pip install`

Install the package as usual and run the console script it provides. This
runs directly in whatever Python environment you installed it into -- no
separate virtual environment is created or managed for you:

```console
pip install SigProfilerPlotting
sigprofilerplotting-web
```

### What happens next

Either way, once the server is ready your default browser opens
automatically to the dashboard, and the terminal prints the URL it's running
at (`http://127.0.0.1:<port>/`) along with a reminder of how to stop it.
If a dashboard is already running, launching again does not start a second
copy -- it prints the existing one's URL and exits.

### Command-line options

All of these work with any of the launcher scripts, `sigprofilerplotting-web`,
or `SigProfilerPlotting dashboard` (the CLI subcommand, if you prefer
invoking it that way):

| Flag | Effect |
| --- | --- |
| `--no-launch` | Start the server without opening a browser tab. |
| `--reinstall` | Force a fresh dependency install even if nothing changed (launcher scripts only -- `sigprofilerplotting-web` never manages its own environment, so this has no effect there). |
| `--stop` | Stop a running dashboard server and exit. |
| `--help` | Print usage and exit. |

## Stopping the dashboard

Any of these will cleanly shut the server down:

* Click the **Quit** button in the top-right corner of the dashboard page.
* Press `Ctrl+C` in the terminal the server is running in.
* Run the launcher (or `sigprofilerplotting-web`) again with `--stop`, from
  another terminal.

## Using the dashboard

1. Click **Choose a matrix file** (or drag a file onto it) to upload a tab-
   or comma-delimited matrix with a mutation-context column and one column
   per sample or signature. The first column's header can be named
   `MutationType`, `Mutation Types`, `classification`, `index`, `context`,
   or simply left blank -- all of these are accepted, matching the different
   ways SigProfiler tools and `pandas` themselves name (or don't name) that
   column.
2. The dashboard detects the mutation context automatically from the
   complete set of row labels -- never from the filename or row count -- and
   lists every sample/signature column in the sidebar.
3. **Choosing what to view:**
    * Use the search box, **Select all**, or **Clear all** to pick which
      samples/signatures to draw.
    * If the matrix has more than 100 samples/signatures, only the first one
      is selected by default so the page stays responsive immediately --
      everyone else is still right there in the list to search for and
      select.
    * Selected items render as you scroll; the dashboard never tries to draw
      every selected sample up front, so selecting hundreds at once stays
      responsive.
4. **Values panel:** switch between **Counts** (raw numbers), **Percentage**
   (share of the sample's total, 0-100%), and **Proportion** (the same
   share, 0-1) at any time -- every open chart redraws immediately.
5. For `SBS4608` and `SBS1536` matrices, an extra panel controls which
   panels of the multi-part figure are shown (see below).
6. For `SBS96` matrices, selecting **exactly one** sample or signature
   enables a **Compare with COSMIC** panel (see below).
7. Every chart supports the standard Plotly toolbar: hover for exact values,
   zoom/pan, reset the view, and export as a PNG or SVG image.

## Comparing against COSMIC signatures

With an `SBS96`, `SBS288`, `SBS1536`, or `SBS4608` matrix loaded and exactly
one sample or signature selected, check **Compare with COSMIC signatures**
in the sidebar and pick a genome build (`GRCh37`, `GRCh38`, `mm9`, `mm10`,
`mm39`, `rn6`, or `rn7` -- COSMIC v3.6). The plot splits into two panels:

* **Upper panel:** your selected sample/signature, unchanged -- including
  `SBS4608`/`SBS1536`'s own multi-part view.
* **Lower panel:** the 5 COSMIC signatures from that build whose spectrum is
  most similar to it, each drawn full size (the same size as the upper
  panel) and labeled with its cosine similarity score. One is shown at a
  time, highest-similarity first; scroll the panel, or use the **‹**/**›**
  buttons, to move through the rest.

`SBS288` (strand-tagged trinucleotide), `SBS1536` (pentanucleotide), and
`SBS4608` (strand-tagged pentanucleotide) are compared by first collapsing
them down to the 96 trinucleotide contexts COSMIC's own catalogue uses --
summing out the strand tag and/or the wider flanking bases. An `SBS96`
matrix needs no collapsing.

The comparison is scale-invariant -- it gives the same ranking regardless of
whether the Values panel is set to Counts, Percentage, or Proportion (the
comparison panels themselves follow Percentage/Proportion, since a
reference signature has no raw count of its own; Counts falls back to
Proportion there). It deactivates automatically (and the checkbox unchecks
itself) as soon as more than one sample/signature is selected, since "most
similar to what" stops being well-defined at that point.

## Supported contexts

The dashboard has a purpose-built, COSMIC-style view for:

* **`SBS96`** -- the classic colored substitution-type header bands and
  trinucleotide-colored axis ticks.
* **`SBS4608`** -- the same multi-panel view as
  [`plotInteractive`](Plotting-SBS4608.md): a 96-bin bar chart, a 5-base
  context ("total") heatmap, front/back detailed heatmaps, and a
  transcription-bias panel (Transcribed/Untranscribed, Genic/Intergenic).
  A **Top panel** control switches the main chart between **Default**,
  **Strand bias**, and **Genic bias**; three **Show** checkboxes toggle the
  **Transcription bias**, **5-base context (total heatmap)**, and **Detailed
  heatmaps (front/back)** sections independently -- hiding a section
  reclaims its space instead of leaving it blank.
* **`SBS1536`** -- the same bar chart, 5-base context heatmap, and detailed
  front/back heatmaps as `SBS4608`'s default view, since a 1536-row matrix
  is exactly what `SBS4608` data becomes once its strand tag is summed away.
  There is no strand information in an `SBS1536` matrix, so the **Top
  panel** variants and the **Transcription bias** section do not apply and
  are not shown.
* **`SBS288`**, **`DBS78`**, **`ID83`**, **`CNV48`**, and **`SV32`** --
  each with its own grouped, correctly colored and labeled layout.

Any other context SigProfilerPlotting recognizes -- `SBS1536`'s sibling
sizes `SBS6`, `SBS18`, `SBS24`, `SBS384`, `SBS6144`, the other `DBS`/`ID`
sizes (`DBS150`, `DBS186`, `DBS1248`, `DBS2400`, `DBS2976`, `ID28`, `ID96`,
`ID332`, `ID415`, `ID8628`), and so on -- is never rejected; it renders with
a generic grouped bar view instead of a specialized one.

## Common input errors

The dashboard reports the same kind of actionable, specific errors as the
plotting functions: duplicate or unexpected mutation-context rows, malformed
or non-numeric cells, negative values, and duplicate or empty column names
are all rejected with a message naming the exact row or column at fault.

## Troubleshooting

- **The browser didn't open:** open the URL printed in the terminal manually
  (`http://127.0.0.1:<port>/`).
- **"Dashboard is already running" / it won't start a second copy:** this is
  expected -- the launcher refuses to start a second server rather than
  silently orphaning the first. Run it with `--stop` first if you actually
  want to restart it.
- **Uploading a large matrix seems slow / shows nothing for a while:** very
  wide matrices (many thousands of samples) take real time to parse and
  upload -- this is normal; the upload has succeeded once the sample list
  populates. There is no progress bar for the upload itself yet.
- **Changes to the package aren't showing up:** run the launcher with
  `--reinstall` to force a fresh dependency install (repository launcher
  scripts only); if you're running the server from a previous launch, stop
  it (`--stop`) and start it again -- code changes require a restart.

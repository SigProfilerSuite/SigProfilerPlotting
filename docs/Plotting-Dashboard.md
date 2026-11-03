# Interactive Dashboard

SigProfilerPlotting includes a local, browser-based dashboard for exploring
any supported matrix interactively -- upload a file, pick samples or
signatures, and get real Plotly charts with hover, zoom, pan, and image
export, without writing any Python.

!!! note
    The dashboard runs entirely on your machine. It binds only to
    `127.0.0.1` and never sends your data anywhere else.

## Installing and launching

There are two ways to run the dashboard, depending on how you got the code.

**From a cloned repository**, use one of the one-click launcher scripts in
the repository root. Each creates (or reuses) its own isolated `.venv-dashboard`
virtual environment and installs the package into it -- only the first time,
or whenever the package's dependencies change:

| Platform | Script |
| --- | --- |
| macOS | Double-click `launch_dashboard.command` |
| Linux / macOS (terminal) | `./launch_dashboard.sh` |
| Windows | Double-click `launch_dashboard.bat`, or run it from a terminal |

**From a `pip install`**, use the installed console command instead -- it
runs directly in your current environment, with no separate virtual
environment:

```console
pip install sigProfilerPlotting
sigprofilerplotting-web
```

Both open your default browser to the dashboard once the server is ready.

### Command-line options

| Flag | Effect |
| --- | --- |
| `--no-launch` | Start the server without opening a browser tab. |
| `--reinstall` | Force a fresh dependency install (launcher scripts only; ignored by `sigprofilerplotting-web`, which never manages its own environment). |
| `--stop` | Stop a running dashboard server and exit. |
| `--help` | Print usage and exit. |

## Stopping the dashboard

Press `Ctrl+C` in the terminal running the server, or run the same launcher
(or `sigprofilerplotting-web`) with `--stop` from another terminal.

## Using the dashboard

1. Click **Browse** (or drag a file in) to upload a tab- or comma-separated
   matrix. The file never leaves your machine -- there is no server-side file
   browser, only your browser's own file picker.
2. The dashboard detects the mutation context automatically from the
   complete set of row labels, not from the filename or row count, and lists
   every sample/signature column. Use the search box, **Select all**, or
   **Clear all** to choose which ones to view.
3. Toggle **Counts** / **Percentage** to switch how values are displayed.
4. Selected samples render as you scroll -- the dashboard does not try to
   render every selected sample up front, so selecting hundreds of samples
   stays responsive.
5. Every chart supports the standard Plotly toolbar: hover for exact values,
   zoom/pan, reset the view, and export as an image.

## Supported contexts

The dashboard has purpose-built views for `SBS96`, `SBS288`, `DBS78`, `ID83`,
`CNV48`, `SV32`, and the strand-tagged `SBS4608` format (the same
multi-panel view as [`plotInteractive`](Plotting-SBS4608.md), including its
default, transcription-strand, genic/intergenic, and heatmap panels).

Any other context SigProfilerPlotting recognizes -- `SBS1536`, `SBS6`,
`ID96`, `DBS1248`, and so on -- is never rejected; it renders with a generic
grouped bar view instead of a specialized one.

## Common input errors

The dashboard reports the same kind of actionable, specific errors as the
plotting functions: duplicate or unexpected mutation-context rows, malformed
or non-numeric cells, negative values, and duplicate or empty column names
are all rejected with a message naming the exact row or column at fault.

## Troubleshooting

- **The browser didn't open:** open the URL printed in the terminal manually
  (`http://127.0.0.1:<port>/`).
- **"Address already in use" or the dashboard won't start:** another
  dashboard instance may already be running -- run the launcher with
  `--stop` first.
- **Changes to the package aren't showing up:** run the launcher with
  `--reinstall` to force a fresh dependency install (repository launcher
  scripts only).

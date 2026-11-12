# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- A "Save plot as PDF" button in each plot's toolbar, next to the existing SVG export, via a new `/api/export/pdf` endpoint that wraps a browser-rasterized PNG in a one-page PDF (no Kaleido/Chrome needed, unlike the static SBS4608 PDF export).
- A redesigned dashboard home screen: an illustrated card per mutation class (substitutions, doublets, indels, copy-number, structural variants), each with a "Try `<context>`" button that loads a real bundled example matrix with a single click, via a new `/api/examples` and `/api/examples/load` endpoint pair. The bundled example matrices are now included in the package itself (`setup.py`/`MANIFEST.in`), so this also works from a plain `pip install`.
- A bundled `dbs78_example.all` example matrix, giving the home screen's "Doublet substitutions" card a working "Try DBS78" button.
- Made the logo/title in the dashboard's top bar a working "Home" button that clears the loaded matrix and returns to the home screen.

### Fixed
- Reject incomplete, duplicated, negative, or non-finite user-supplied COSMIC catalogues and ambiguous catalogue versions.
- Keep COSMIC catalogue metadata synchronized when a file is replaced in place.
- Prevent uploaded matrix replacement from deleting a backing file during an active dashboard render.
- Preserve dashboard socket errors during startup and support IPv6 loopback hosts correctly.
- Parse repeated CLI custom-text options as per-sample values instead of individual characters.
- Accept R row-name matrix exports whose header omits the mutation-context column, inferred only when the data proves the shift, so a genuinely wrong first column is still rejected.
- Canonicalize the legacy ID89 `Del(C):R(6,9)` repeat-interval label to the reference `Del(C):R(6,)`.
- Omit, with a visible dashboard warning, ID476 signature columns R-padded with `NA` for an unresolved or entirely missing context instead of zero-filling them; reject any other `NA` pattern as a corrupted matrix.
- Stopped the dashboard's "Samples and signatures" sidebar panel from overlapping the citation footer when the COSMIC and SBS4608 view panels were both visible at once and their combined height left the sample list no room to lay out.
- Fixed a dead dashboard CSS rule that could never hide the sidebar at narrow widths (a same-specificity override declared before, instead of after, the base rule it was meant to win against); replaced it with a rule that keeps the sidebar's matrix/sample controls reachable by stacking it above the plots instead of hiding it, since nothing in the main content duplicates them.
- Fixed DBS1248's dashboard renderer grouping bars by flank prefix instead of by reference dinucleotide, and DBS78's renderer reusing colors across different reference dinucleotides.
- Fixed the dashboard's "Compare with COSMIC signatures" checkbox staying stuck disabled, with a stale hint, if a matrix was uploaded and narrowed to one sample before the one-time COSMIC catalogue scan at page load had finished.
- Removed a stray trailing tab from every line of the bundled `breast_cancer_samples_example.SBS96.all`, `.DBS78.all`, and `.SBS192.all` example matrices, which made the dashboard's strict parser reject them with "every column needs a non-empty name".
- Fixed the dashboard launcher's Windows process-liveness check calling `os.kill(pid, 0)`, which on Windows sends Ctrl+C to every process on the console (the caller included) instead of probing the pid, since signal 0 is `CTRL_C_EVENT` there; it could raise `KeyboardInterrupt` out of `--stop` whenever it found the startup lock still held by the server it had just stopped.

### Changed
- Updated repository, maintainer, and contact metadata.
- The dashboard's plots now default to Percentage on the y-axis instead of Counts, for every mutation context.
- Restyled the interactive dashboard to match the shared SigProfiler Suite workbench design: four selectable themes (UC San Diego, UAlbany, SigProfilerSuite, Night), a glossy sidebar, card-based plot layout, status pills, and "?" explainer panels.

## [1.5.0] - 2026-11-06

### Added
- Added static and interactive SBS4608 plotting through `plotSBS`.
- Added a local interactive dashboard with one-click launchers for macOS, Linux, and Windows.
- Added counts, percentage, and proportion views, plot export, and optional comparison with user-supplied COSMIC SBS catalogues.
- Added static plotting support for the ID89 and ID476 indel classifications.
- Added tagged TestPyPI and GitHub Release-to-PyPI publication workflows.

### Changed
- Unified SBS4608 static plotting under `plotSBS` and removed the transitional `plotStatic` entry point.
- Streamed dashboard matrix uploads to disk and made sample rendering lazy for large matrices.
- Expanded MkDocs and README guidance for SBS4608, the dashboard, indel plots, and releases.

### Fixed
- Corrected SBS4608 PDF export, scaling, labels, heatmap colorbars, and compact interactive layouts.
- Hardened dashboard upload, process identity, shutdown, and relaunch behavior.
- Corrected ID89/ID476 labels, peak annotations, complex-indel handling, and plot widths.

## [1.4.4] - 2026-10-02

### Added
- Migrated CI, documentation deployment, tagging, and publication automation to GitHub Actions.

### Fixed
- Made image-baseline comparisons portable across supported runners.
- Added input-validation and plotting regression coverage.
- Updated project organization links and contact information.

## [1.4.3] - 2026-01-22

### Fixed
- Fixed KeyError when accessing groupby groups in `tmbplot.py` for pandas 3.12 compatibility
- Fixed integer indexing on pandas Series in `tmbplot.py` by converting to lists for pandas 3.12 compatibility
- Improved error handling in `plotTMB()` function to handle group access failures gracefully

## [1.4.2] - 2025-10-13

### Added
- Implemented a CI/CD pipeline with Travis CI to automate the building and publishing of Docker images to Docker Hub.
- Added a `pyproject.toml` file to standardize the package build system.

## [1.4.1] - 2025-03-10

### Fixed
- The plotting function for structural variants has had syntax to reflect pandas 2.0.0 changes.

## [1.4.0] - 2025-02-11

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**

## [1.3.24] - 2024-08-01

### Added
- Added environmental variable `SIGPROFILERPLOTTING_VOLUME` to enhance configurability.

### Fixed
- Ensured correct file path construction with `os.path.join()`.
- Updated CLI to handle boolean arguments using `str2bool`, resolving an issue where boolean strings like `False` were incorrectly evaluated as `True`.
- Stopped unecessary output directory creation when the output format is PIL_Image.

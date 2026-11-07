# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- Reject incomplete, duplicated, negative, or non-finite user-supplied COSMIC catalogues and ambiguous catalogue versions.
- Keep COSMIC catalogue metadata synchronized when a file is replaced in place.
- Prevent uploaded matrix replacement from deleting a backing file during an active dashboard render.
- Preserve dashboard socket errors during startup and support IPv6 loopback hosts correctly.
- Parse repeated CLI custom-text options as per-sample values instead of individual characters.

### Changed
- Updated repository, maintainer, and contact metadata.

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

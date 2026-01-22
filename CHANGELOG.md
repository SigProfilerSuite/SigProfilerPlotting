# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

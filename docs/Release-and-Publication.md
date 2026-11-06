# Release and publication

SigProfilerPlotting uses Git tags as its version source. `setuptools-scm`
derives the package version from tags matching `vX.Y.Z`; maintainers should
not manually edit a version constant for a release.

## Before tagging

1. Merge the intended release changes into the target branch, normally
   `master`.
2. Confirm CI passes for Python 3.9 and Python 3.12.
3. Confirm the documentation builds without warnings.
4. For changes to SBS4608 static output, run the real renderer check:

   ```console
   SPP_RUN_KALEIDO_INTEGRATION=1 pytest -q tests/test_sbs4608_plotting.py -k real_kaleido
   ```

5. Choose a version that has not already been used on GitHub, TestPyPI, or
   PyPI. Published package versions cannot be replaced.

## Create the version tag and TestPyPI release

1. Open **Actions > Create version tag > Run workflow**.
2. Enter the version without a leading `v`, for example `1.5.0`.
3. Select the branch to tag, normally `master`.
4. Run the workflow.

The workflow creates and pushes the annotated `vX.Y.Z` tag, then dispatches
the **Release** workflow. That workflow:

- builds the source distribution and wheel;
- validates both with Twine;
- installs the wheel and verifies the derived version;
- stores the distributions as a workflow artifact; and
- publishes the tagged build to TestPyPI.

TestPyPI publishing uses the protected GitHub environment named `testpypi`
and PyPI trusted publishing (OIDC), not a stored API token.

## Publish to PyPI

After validating installation and basic plotting from TestPyPI, create and
publish a GitHub Release for the same `vX.Y.Z` tag. The published-release
event runs the verified build again and publishes it to PyPI through the
protected GitHub environment named `pypi` and PyPI trusted publishing.

The PyPI trusted publisher must match exactly:

| Setting | Value |
| --- | --- |
| Repository owner | `SigProfilerSuite` |
| Repository name | `SigProfilerPlotting` |
| Workflow | `release.yml` |
| Environment | `pypi` |

Configure the equivalent TestPyPI publisher with environment `testpypi`.

## Post-release verification

1. Confirm the GitHub Actions release run is green.
2. Install the release into a clean environment:

   ```console
   python -m pip install --upgrade SigProfilerPlotting==X.Y.Z
   python -c "import sigProfilerPlotting as spp; print(spp.__version__)"
   ```

3. Generate at least one standard plot and, when applicable, one SBS4608
   static PDF or dashboard plot.
4. Confirm the GitHub Release, PyPI project, and documentation links point to
   the intended version.

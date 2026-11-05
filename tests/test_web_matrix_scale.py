"""Memory budget for SigProfilerPlotting's own largest matrix shape.

plotting_A's review measured the original (pre-fix) upload path at ~20x
the input size in peak memory for a 4.67 MB synthetic SBS4608 payload
(4,608 rows x 500 samples), and projected that scaling linearly to the
user's real matrix (~51,021 samples) could need 9-10+ GB and was likely to
freeze an ordinary laptop. The fix (sigProfilerPlotting/web/matrix.py:
MatrixData reading lazily from disk, one sample at a time, instead of
holding every sample's floats in memory at once) is meant to make peak
memory roughly constant in the number of samples, not linear in it.

This test is deliberately larger than the dashboard's other matrix tests
(2,000 SBS4608 columns, ~27 MB of text) and takes real wall-clock time
(several seconds) to write and parse -- that is the point: it is the one
place in the suite that actually exercises the shape of matrix the review
was concerned about, rather than a handful of columns. The documented
budget: peak traced Python memory for loading and for reading any single
sample must stay under 20 MB, regardless of how many sample columns the
matrix has -- i.e. it must not scale with sample count at all.
"""

import random
import tracemalloc

import pytest

from sigProfilerPlotting.web.matrix import load_matrix

_MEMORY_BUDGET_BYTES = 20 * 1024 * 1024  # 20 MB, documented above
_COLUMN_COUNT = 2000


@pytest.fixture(scope="module")
def large_sbs4608_path(tmp_path_factory):
    import sigProfilerPlotting.web.matrix as matrix_module

    labels = list(matrix_module.REFERENCE_CONTEXTS["SBS4608"])
    rng = random.Random(0)
    path = tmp_path_factory.mktemp("scale") / "large_sbs4608.tsv"
    columns = [f"S{i}" for i in range(_COLUMN_COUNT)]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("MutationType\t" + "\t".join(columns) + "\n")
        for label in labels:
            row = [str(rng.randint(0, 100)) for _ in columns]
            # Make exactly one cell per row identifiable so a later test can
            # confirm the right column was actually read, not just that
            # *some* value came back.
            row[-1] = "999"
            handle.write(label + "\t" + "\t".join(row) + "\n")
    return path, columns


def test_loading_a_large_sbs4608_matrix_stays_within_the_memory_budget(large_sbs4608_path):
    path, _columns = large_sbs4608_path
    tracemalloc.start()
    try:
        matrix = load_matrix(path)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert matrix.context == "SBS4608"
    assert len(matrix.columns) == _COLUMN_COUNT
    assert peak < _MEMORY_BUDGET_BYTES, (
        f"loading a {_COLUMN_COUNT}-column SBS4608 matrix used {peak / 1e6:.1f} MB "
        f"(budget: {_MEMORY_BUDGET_BYTES / 1e6:.0f} MB) -- memory must not scale with "
        "sample count."
    )


def test_rendering_one_sample_from_a_large_matrix_stays_within_the_memory_budget(large_sbs4608_path):
    path, columns = large_sbs4608_path
    matrix = load_matrix(path)

    tracemalloc.start()
    try:
        values = matrix.sample_values(columns[-1])  # the "999"-marked column
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert len(values) == len(matrix.labels)
    assert all(v == 999.0 for v in values)
    assert peak < _MEMORY_BUDGET_BYTES, (
        f"rendering one sample from a {_COLUMN_COUNT}-column SBS4608 matrix used "
        f"{peak / 1e6:.1f} MB (budget: {_MEMORY_BUDGET_BYTES / 1e6:.0f} MB) -- this must "
        "cost roughly the same whether the matrix has 2,000 or 51,000 columns."
    )


def test_rendering_a_different_sample_returns_its_own_distinct_values(large_sbs4608_path):
    path, columns = large_sbs4608_path
    matrix = load_matrix(path)
    first_column_values = matrix.sample_values(columns[0])
    last_column_values = matrix.sample_values(columns[-1])
    assert all(v == 999.0 for v in last_column_values)
    assert not all(v == 999.0 for v in first_column_values)

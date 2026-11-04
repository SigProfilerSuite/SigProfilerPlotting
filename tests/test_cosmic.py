import pytest

from sigProfilerPlotting.web import cosmic
from sigProfilerPlotting.web.matrix import MatrixError


def test_available_builds_lists_the_bundled_genome_builds():
    builds = cosmic.available_builds()
    assert builds == ["GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6", "rn7"]


def test_cosine_similarity_of_identical_vectors_is_one():
    assert cosmic._cosine_similarity([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == pytest.approx(1.0)


def test_cosine_similarity_is_scale_invariant():
    a = [1.0, 2.0, 3.0, 0.0]
    b = [2.0, 4.0, 6.0, 0.0]  # a scaled by 2
    assert cosmic._cosine_similarity(a, b) == pytest.approx(1.0)


def test_cosine_similarity_of_orthogonal_vectors_is_zero():
    assert cosmic._cosine_similarity([1.0, 0.0], [0.0, 1.0]) == pytest.approx(0.0)


def test_cosine_similarity_handles_an_all_zero_vector():
    assert cosmic._cosine_similarity([0.0, 0.0], [1.0, 2.0]) == 0.0


def test_top_matches_ranks_a_signature_as_its_own_best_match():
    # A query identical to one bundled signature must rank that signature
    # first, with a cosine similarity of (numerically) 1.0.
    signatures = cosmic._load_build("GRCh37")
    query = signatures["SBS4"]

    matches = cosmic.top_matches(query, "GRCh37", n=5)

    assert len(matches) == 5
    assert matches[0]["name"] == "SBS4"
    assert matches[0]["score"] == pytest.approx(1.0, abs=1e-9)
    # Descending order.
    scores = [m["score"] for m in matches]
    assert scores == sorted(scores, reverse=True)


def test_top_matches_respects_n():
    signatures = cosmic._load_build("GRCh37")
    query = signatures["SBS1"]
    assert len(cosmic.top_matches(query, "GRCh37", n=3)) == 3


def test_top_matches_rejects_an_unknown_build():
    with pytest.raises(MatrixError, match="Unknown COSMIC genome build"):
        cosmic.top_matches({"A[C>A]A": 1.0}, "not-a-real-build")


def test_load_build_caches_across_calls():
    first = cosmic._load_build("mm10")
    second = cosmic._load_build("mm10")
    assert first is second

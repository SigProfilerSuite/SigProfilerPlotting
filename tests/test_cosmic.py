import importlib

import pytest

from sigProfilerPlotting.web import cosmic
from sigProfilerPlotting.web.matrix import MatrixError

# SigProfilerPlotting does not bundle or redistribute COSMIC's signature
# catalogues (see cosmic.py's module docstring) -- every test here builds
# its own small, synthetic, user-supplied-style catalogue directory rather
# than relying on real bundled data.

_HEADER = ["Type", "SBS1", "SBS4", "SBS5"]
_ROWS = [
    _HEADER,
    ["A[C>A]A", "1", "0", "2"],
    ["A[C>G]A", "0", "1", "2"],
    ["A[C>T]A", "0", "0", "2"],
    ["A[T>A]A", "0", "0", "2"],
]


def _write_catalogue(path, rows=_ROWS):
    content = "\n".join("\t".join(row) for row in rows) + "\n"
    path.write_text(content, encoding="utf-8")


@pytest.fixture
def catalogue_dir(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(directory / "COSMIC_v3.4_SBS_GRCh37.txt")
    _write_catalogue(directory / "COSMIC_v3.4_SBS_mm10.txt")
    cosmic.configure_catalogue_dir(directory)
    yield directory
    cosmic.configure_catalogue_dir(None)


def test_available_builds_is_empty_when_nothing_is_configured():
    cosmic.configure_catalogue_dir(None)
    assert cosmic.available_builds() == []


def test_available_builds_lists_builds_found_in_the_configured_directory(catalogue_dir):
    assert cosmic.available_builds() == ["GRCh37", "mm10"]


def test_available_builds_is_empty_for_a_directory_that_does_not_exist(tmp_path):
    cosmic.configure_catalogue_dir(tmp_path / "does-not-exist")
    assert cosmic.available_builds() == []
    cosmic.configure_catalogue_dir(None)


def test_unrecognized_filenames_in_the_directory_are_ignored(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(directory / "COSMIC_v3.4_SBS_GRCh37.txt")
    (directory / "readme.txt").write_text("not a catalogue\n", encoding="utf-8")
    (directory / "COSMIC_v3.4_SBS_notabuild.txt").write_text("Type\tSBS1\n", encoding="utf-8")
    cosmic.configure_catalogue_dir(directory)
    try:
        assert cosmic.available_builds() == ["GRCh37"]
    finally:
        cosmic.configure_catalogue_dir(None)


def test_configure_catalogue_dir_none_clears_it(catalogue_dir):
    assert cosmic.available_builds() != []
    cosmic.configure_catalogue_dir(None)
    assert cosmic.available_builds() == []
    cosmic.configure_catalogue_dir(catalogue_dir)  # restore for fixture teardown parity


def test_catalogue_metadata_reports_version_path_checksum_and_signature_count(catalogue_dir):
    metadata = {info.build: info for info in cosmic.catalogue_metadata()}
    info = metadata["GRCh37"]
    assert info.version == "3.4"
    assert info.source_path == str(catalogue_dir / "COSMIC_v3.4_SBS_GRCh37.txt")
    assert info.signature_count == 3
    import hashlib

    expected = hashlib.sha256((catalogue_dir / "COSMIC_v3.4_SBS_GRCh37.txt").read_bytes()).hexdigest()
    assert info.checksum == expected


def test_load_build_rejects_a_header_with_no_signature_columns(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(directory / "COSMIC_v3.4_SBS_GRCh37.txt", rows=[["Type"]])
    cosmic.configure_catalogue_dir(directory)
    try:
        with pytest.raises(MatrixError, match="header row"):
            cosmic._load_build("GRCh37")
    finally:
        cosmic.configure_catalogue_dir(None)


def test_load_build_rejects_an_unnamed_signature_column(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(
        directory / "COSMIC_v3.4_SBS_GRCh37.txt",
        rows=[["Type", "SBS1", ""], ["A[C>A]A", "1", "2"]],
    )
    cosmic.configure_catalogue_dir(directory)
    try:
        with pytest.raises(MatrixError, match="unnamed signature column"):
            cosmic._load_build("GRCh37")
    finally:
        cosmic.configure_catalogue_dir(None)


def test_load_build_rejects_a_row_with_the_wrong_number_of_fields(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(
        directory / "COSMIC_v3.4_SBS_GRCh37.txt",
        rows=[["Type", "SBS1"], ["A[C>A]A", "1", "2"]],
    )
    cosmic.configure_catalogue_dir(directory)
    try:
        with pytest.raises(MatrixError, match="has 3 fields"):
            cosmic._load_build("GRCh37")
    finally:
        cosmic.configure_catalogue_dir(None)


def test_load_build_rejects_a_non_numeric_value(tmp_path):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(
        directory / "COSMIC_v3.4_SBS_GRCh37.txt",
        rows=[["Type", "SBS1"], ["A[C>A]A", "not-a-number"]],
    )
    cosmic.configure_catalogue_dir(directory)
    try:
        with pytest.raises(MatrixError, match="non-numeric"):
            cosmic._load_build("GRCh37")
    finally:
        cosmic.configure_catalogue_dir(None)


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


def test_top_matches_ranks_a_signature_as_its_own_best_match(catalogue_dir):
    # A query identical to one catalogue signature must rank that signature
    # first, with a cosine similarity of (numerically) 1.0.
    signatures = cosmic._load_build("GRCh37")
    query = signatures["SBS4"]

    matches = cosmic.top_matches(query, "GRCh37", n=3)

    assert len(matches) == 3
    assert matches[0]["name"] == "SBS4"
    assert matches[0]["score"] == pytest.approx(1.0, abs=1e-9)
    # Descending order.
    scores = [m["score"] for m in matches]
    assert scores == sorted(scores, reverse=True)


def test_top_matches_respects_n(catalogue_dir):
    signatures = cosmic._load_build("GRCh37")
    query = signatures["SBS1"]
    assert len(cosmic.top_matches(query, "GRCh37", n=2)) == 2


def test_top_matches_rejects_an_unknown_build(catalogue_dir):
    with pytest.raises(MatrixError, match="Unknown or unavailable COSMIC genome build"):
        cosmic.top_matches({"A[C>A]A": 1.0}, "not-a-real-build")


def test_top_matches_rejects_a_build_removed_since_available_builds_was_checked(catalogue_dir):
    with pytest.raises(MatrixError, match="Unknown or unavailable COSMIC genome build"):
        cosmic.top_matches({"A[C>A]A": 1.0}, "rn7")  # never configured in this fixture


def test_load_build_caches_across_calls(catalogue_dir):
    first = cosmic._load_build("mm10")
    second = cosmic._load_build("mm10")
    assert first is second


def test_load_build_cache_invalidates_when_the_file_content_changes(catalogue_dir):
    first = cosmic._load_build("GRCh37")
    assert first["SBS1"]["A[C>A]A"] == 1.0
    _write_catalogue(
        catalogue_dir / "COSMIC_v3.4_SBS_GRCh37.txt",
        rows=[["Type", "SBS1"], ["A[C>A]A", "99"]],
    )
    second = cosmic._load_build("GRCh37")
    assert second["SBS1"]["A[C>A]A"] == 99.0


def test_collapse_to_sbs96_is_a_no_op_for_already_96_context_labels():
    raw = {"A[C>A]A": 3.0, "C[C>G]T": 5.0}
    assert cosmic.collapse_to_sbs96(raw) == raw


def test_collapse_to_sbs96_strips_a_strand_tag():
    raw = {"T:A[C>A]A": 4.0, "U:A[C>A]A": 2.0, "N:A[C>A]A": 1.0}
    assert cosmic.collapse_to_sbs96(raw) == {"A[C>A]A": 7.0}


def test_collapse_to_sbs96_sums_wider_pentanucleotide_flanks():
    raw = {"AA[C>A]AA": 10.0, "CA[C>A]AA": 5.0, "AA[C>A]CC": 3.0}
    assert cosmic.collapse_to_sbs96(raw) == {"A[C>A]A": 15.0, "A[C>A]C": 3.0}


def test_collapse_to_sbs96_strips_strand_tag_and_pentanucleotide_flanks_together():
    raw = {"T:AA[C>A]AA": 10.0, "U:CA[C>A]AA": 5.0, "N:AA[C>A]CC": 3.0}
    assert cosmic.collapse_to_sbs96(raw) == {"A[C>A]A": 15.0, "A[C>A]C": 3.0}


def test_environment_variable_configures_the_default_catalogue_dir(tmp_path, monkeypatch):
    directory = tmp_path / "cosmic"
    directory.mkdir()
    _write_catalogue(directory / "COSMIC_v3.4_SBS_GRCh37.txt")
    monkeypatch.setenv("SIGPROFILERPLOTTING_COSMIC_DIR", str(directory))
    importlib.reload(cosmic)
    try:
        assert cosmic.available_builds() == ["GRCh37"]
    finally:
        # Unset before reloading again so the module genuinely goes back to
        # its un-configured state, rather than reloading while the env var
        # (undone only later, by monkeypatch's own fixture teardown) is
        # still in effect.
        monkeypatch.delenv("SIGPROFILERPLOTTING_COSMIC_DIR", raising=False)
        importlib.reload(cosmic)

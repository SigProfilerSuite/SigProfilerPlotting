"""The bundled "Try an example" menu's own registry (sigProfilerPlotting.web.
examples), independent of the HTTP layer tested in test_web_server.py.
"""

from sigProfilerPlotting.web import examples


def test_every_bundled_example_is_present_and_loads_cleanly():
    # A guard against silent drift: if a filename here were ever renamed, or
    # a newly added entry's file didn't actually pass the dashboard's own
    # strict parser, this is where that would be caught, rather than a
    # reader simply seeing fewer (or zero) "Try an example" buttons.
    found = examples.available()
    names = {entry["name"] for entry in found}
    assert names == {"sbs96", "sbs288", "dbs78", "id83", "cnv48", "sv32"}
    contexts = {entry["name"]: entry["context"] for entry in found}
    assert contexts == {
        "sbs96": "SBS96",
        "sbs288": "SBS288",
        "dbs78": "DBS78",
        "id83": "ID83",
        "cnv48": "CNV48",
        "sv32": "SV32",
    }
    for entry in found:
        assert entry["label"]


def test_example_path_resolves_known_names_to_real_files():
    for name in ("sbs96", "sbs288", "dbs78", "id83", "cnv48", "sv32"):
        path = examples.example_path(name)
        assert path is not None
        assert path.is_file()


def test_example_path_is_none_for_an_unknown_name():
    assert examples.example_path("not-a-real-example") is None

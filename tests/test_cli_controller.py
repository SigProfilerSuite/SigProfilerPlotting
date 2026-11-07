import pytest

from sigProfilerPlotting.controllers import cli_controller


BASE_ARGS = ["matrix.tsv", "plots", "Project", "4608"]


def test_sbs_cli_percentage_omission_is_none():
    assert cli_controller.parse_arguments_sbs(BASE_ARGS).percentage is None


@pytest.mark.parametrize(
    ("value", "expected"),
    [("false", False), ("true", True)],
)
def test_sbs_cli_preserves_explicit_percentage(value, expected):
    parsed = cli_controller.parse_arguments_sbs([*BASE_ARGS, "--percentage", value])
    assert parsed.percentage is expected


def test_id_and_dbs_cli_percentage_defaults_remain_false():
    args = ["matrix.tsv", "plots", "Project", "83"]
    assert cli_controller.parse_arguments_id(args).percentage is False
    assert cli_controller.parse_arguments_dbs(args).percentage is False


def test_sbs4608_cli_omission_reaches_percentage_true(monkeypatch):
    import sigProfilerPlotting.sbs4608_static as static_module

    captured = {}

    def fake_build(*args, **kwargs):
        captured["positional"] = args
        captured.update(kwargs)
        return "plots/SBS_4608_plots_Project.pdf"

    monkeypatch.setattr(static_module, "_build_sbs4608_pdf", fake_build)
    parsed = cli_controller.parse_arguments_sbs(BASE_ARGS)
    cli_controller.dispatch_plot_sbs(parsed)

    assert captured["percentage"] is True


def test_cli_custom_text_is_a_per_sample_list_instead_of_characters():
    parsed = cli_controller.parse_arguments_sbs(
        [
            *BASE_ARGS,
            "--custom_text_upper",
            "Stability: 0.97",
            "--custom_text_upper",
            "Stability: 0.88",
            "--custom_text_middle",
            "Mutations: 12,345",
        ]
    )

    assert parsed.custom_text_upper == ["Stability: 0.97", "Stability: 0.88"]
    assert parsed.custom_text_middle == ["Mutations: 12,345"]
    assert parsed.custom_text_bottom is None

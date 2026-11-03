#!/usr/bin/env python3

import sys
from sigProfilerPlotting import sigProfilerPlotting
from sigProfilerPlotting.controllers import cli_controller


def main_function():
    commands = {
        "plotSBS": "Plot Single Base Substitutions.",
        "plotID": "Plot Small Insertions and Deletions.",
        "plotDBS": "Plot Doublet Base Substitutions.",
        "plotSV": "Plot Structural Variations.",
        "plotCNV": "Plot Copy Number Variations.",
        "dashboard": "Launch the local interactive dashboard.",
    }

    if len(sys.argv) < 2 or sys.argv[1] not in commands.keys():
        print_usage(commands)
        sys.exit(1)

    command = sys.argv[1]
    args = sys.argv[1:]

    controller = cli_controller.CliController()
    valid_commands = {"plotSBS", "plotID", "plotDBS", "plotSV", "plotCNV", "dashboard"}

    if command in valid_commands:
        controller.dispatch(args)
    else:
        print_usage(commands)
        sys.exit(1)


def print_usage(commands):
    """Prints the usage message."""
    print("Usage: SigProfilerPlotting <command> [<args>]\n")
    print("Commands:")
    for cmd, desc in commands.items():
        print(f"  {cmd}: {desc}")
    print(
        "\nUse 'SigProfilerPlotting <command> --help' for more information on a specific command."
    )


if __name__ == "__main__":
    main_function()

#!/usr/bin/env python3
"""Public PHap wrapper for the v2 phase-reads workflow."""

import os
import subprocess
import sys


def main():
    script_directory = os.path.dirname(os.path.realpath(__file__))
    utility = os.path.join(
        script_directory, "..", "utils", "phase_reads_assemble_anchor.py"
    )
    subprocess.run([sys.executable, utility, *sys.argv[1:]], check=True)


if __name__ == "__main__":
    main()

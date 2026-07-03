#!/usr/bin/env python3
"""Compatibility entry point for the installable ``phap dosage`` command."""

from __future__ import annotations

import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.phap_dosage import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())

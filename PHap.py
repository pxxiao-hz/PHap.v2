#!/usr/bin/env python3
"""Backward-compatible source-checkout entry point for PHap."""

from phap_core.cli import main


if __name__ == "__main__":
    raise SystemExit(main())

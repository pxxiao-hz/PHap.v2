"""Command-line wrapper for window-level dosage inference."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path
from typing import Iterator, Optional, Sequence, TextIO

from phap_core.dosage import (
    DosageError,
    DosageModel,
    UnitigDosageCall,
    WindowDosageCall,
    classify_windows,
    fit_dosage_model,
    parse_window_depths,
    summarize_unitigs,
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse dosage command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Infer unitig dosage from window-level read depth."
    )
    parser.add_argument(
        "--input-file",
        "--input_file",
        dest="input_file",
        type=Path,
        required=True,
        help="Tab-separated window-depth table.",
    )
    parser.add_argument(
        "--ploidy",
        type=int,
        required=True,
        help="Maximum biological copy state/haplotype count.",
    )
    parser.add_argument(
        "--haploid-depth",
        default="auto",
        metavar="auto|FLOAT",
        help="Estimate from all windows (auto) or use a positive depth.",
    )
    parser.add_argument(
        "--relative-sigma",
        type=float,
        help="Override the fitted relative depth standard deviation (0, 0.5].",
    )
    parser.add_argument(
        "--min-confidence",
        type=float,
        default=0.8,
        help="Minimum posterior probability for a window dosage call (default: 0.8).",
    )
    parser.add_argument(
        "--min-unitig-support",
        type=float,
        default=0.8,
        help=(
            "Minimum fraction of all windows supporting the dominant unitig "
            "class (default: 0.8)."
        ),
    )
    parser.add_argument(
        "--pandepth",
        action="store_true",
        help="Use PanDepth defaults: contig/start/end columns 1/2/3 and depth column 8.",
    )
    parser.add_argument(
        "--contig-column",
        type=int,
        default=1,
        help="One-based contig ID column (default: 1).",
    )
    parser.add_argument(
        "--start-column",
        type=int,
        default=2,
        help="One-based window start column (default: 2).",
    )
    parser.add_argument(
        "--end-column",
        type=int,
        default=3,
        help="One-based window end column (default: 3).",
    )
    parser.add_argument(
        "--depth-column",
        type=int,
        help="One-based depth column (default: 8 with --pandepth, otherwise 6).",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Skip the first non-comment line as a header.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory).",
    )
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """Run dosage inference and return its three output paths."""

    input_path: Path = args.input_file
    if not input_path.is_file():
        raise DosageError(f"missing input file: {input_path}")
    if args.ploidy < 1:
        raise DosageError("ploidy must be at least one")
    for option_name in ("contig_column", "start_column", "end_column"):
        if getattr(args, option_name) < 1:
            raise DosageError(f"{option_name.replace('_', '-')} must be at least one")
    depth_column = args.depth_column
    if depth_column is None:
        depth_column = 8 if args.pandepth else 6
    if depth_column < 1:
        raise DosageError("depth-column must be at least one")

    haploid_depth = _parse_haploid_depth(args.haploid_depth)
    with input_path.open(encoding="utf-8") as input_handle:
        windows = parse_window_depths(
            input_handle,
            contig_column=args.contig_column - 1,
            start_column=args.start_column - 1,
            end_column=args.end_column - 1,
            depth_column=depth_column - 1,
            skip_header=args.header,
        )
    model = fit_dosage_model(
        [window.depth for window in windows],
        ploidy=args.ploidy,
        haploid_depth=haploid_depth,
        relative_sigma=args.relative_sigma,
    )
    if haploid_depth is None and model.estimation_status == "weakly_identified":
        alternatives = ", ".join(
            f"{alternative.haploid_depth:.6g}" for alternative in model.alternatives[:3]
        )
        raise DosageError(
            "automatic haploid depth is weakly identified"
            f" (best {model.haploid_depth:.6g}; alternatives: {alternatives});"
            " rerun with --haploid-depth FLOAT"
        )
    window_calls = classify_windows(
        windows,
        model,
        min_confidence=args.min_confidence,
    )
    unitig_calls = summarize_unitigs(
        window_calls,
        min_support=args.min_unitig_support,
    )

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    unitig_path = output_dir / "contig_depth.txt"
    window_path = output_dir / "dosage_windows.tsv"
    model_path = output_dir / "dosage_model.json"
    _write_unitig_calls(unitig_path, unitig_calls)
    _write_window_calls(window_path, window_calls)
    _write_model(model_path, model, args)
    return unitig_path, window_path, model_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the dosage CLI."""

    args = parse_args(argv)
    try:
        unitig_path, window_path, model_path = run(args)
    except (DosageError, OSError) as exc:
        print(f"phap dosage: error: {exc}", file=sys.stderr)
        return 1
    print(f"Unitig calls: {unitig_path}")
    print(f"Window audit: {window_path}")
    print(f"Model audit: {model_path}")
    return 0


def _parse_haploid_depth(value: str) -> Optional[float]:
    if value.lower() == "auto":
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise DosageError("haploid-depth must be 'auto' or a positive number") from exc
    if parsed <= 0:
        raise DosageError("haploid-depth must be positive")
    return parsed


@contextmanager
def _atomic_text_writer(path: Path) -> Iterator[TextIO]:
    temporary_name: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def _write_unitig_calls(path: Path, calls: Sequence[UnitigDosageCall]) -> None:
    with _atomic_text_writer(path) as handle:
        handle.write(
            "contig_ID\taverage_depth\tcontig_type\tdosage\tstatus\t"
            "dominant_class\tdominant_fraction\tmixed_dosage\twindow_count\tclass_counts\n"
        )
        for call in calls:
            dosage = "." if call.dosage is None else str(call.dosage)
            counts = json.dumps(
                dict(call.class_counts),
                ensure_ascii=True,
                sort_keys=True,
                separators=(",", ":"),
            )
            handle.write(
                f"{call.contig_id}\t{call.average_depth:.6f}\t{call.contig_type}\t"
                f"{dosage}\t{call.status}\t{call.dominant_class}\t"
                f"{call.dominant_fraction:.6f}\t"
                f"{str(call.mixed_dosage).lower()}\t{call.window_count}\t{counts}\n"
            )


def _write_window_calls(path: Path, calls: Sequence[WindowDosageCall]) -> None:
    with _atomic_text_writer(path) as handle:
        handle.write(
            "contig_ID\tstart\tend\tdepth\tnormalized_depth\t"
            "classification\tdosage\tconfidence\n"
        )
        for call in calls:
            dosage = "." if call.dosage is None else str(call.dosage)
            handle.write(
                f"{call.contig_id}\t{call.start}\t{call.end}\t{call.depth:.6f}\t"
                f"{call.normalized_depth:.6f}\t{call.classification}\t"
                f"{dosage}\t{call.confidence:.6f}\n"
            )


def _write_model(path: Path, model: DosageModel, args: argparse.Namespace) -> None:
    payload = asdict(model)
    payload["input"] = {
        "path": str(args.input_file),
        "pandepth": bool(args.pandepth),
        "contig_column": args.contig_column,
        "start_column": args.start_column,
        "end_column": args.end_column,
        "depth_column": args.depth_column or (8 if args.pandepth else 6),
        "header": bool(args.header),
    }
    payload["classification"] = {
        "low_depth_policy": "dosage_1",
        "min_confidence": args.min_confidence,
        "min_unitig_support": args.min_unitig_support,
        "unitig_summary_policy": "dominant_class",
    }
    with _atomic_text_writer(path) as handle:
        json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())

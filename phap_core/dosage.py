"""Deterministic dosage inference from window-level read depth.

The functions in this module are pure: they do not read files, write files, or
invoke external programs.  Window coordinates are preserved exactly as they
appear in the input; no coordinate-system conversion is performed here.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple


class DosageError(ValueError):
    """Raised when dosage input or model parameters are invalid."""


@dataclass(frozen=True)
class WindowDepth:
    """Depth observed for one unitig window."""

    contig_id: str
    start: int
    end: int
    depth: float


@dataclass(frozen=True)
class DosageAlternative:
    """A separated alternative solution considered during automatic fitting."""

    haploid_depth: float
    score: float


@dataclass(frozen=True)
class DosageModel:
    """Parameters of a constrained integer-copy depth model."""

    ploidy: int
    haploid_depth: float
    relative_sigma: float
    estimation_method: str
    estimation_status: str
    objective_score: float
    separation: float
    fitted_window_count: int
    alternatives: Tuple[DosageAlternative, ...]


@dataclass(frozen=True)
class WindowDosageCall:
    """Auditable classification of one depth window."""

    contig_id: str
    start: int
    end: int
    depth: float
    normalized_depth: float
    classification: str
    dosage: Optional[int]
    confidence: float


@dataclass(frozen=True)
class UnitigDosageCall:
    """Unitig summary which retains evidence of heterogeneous windows."""

    contig_id: str
    average_depth: float
    contig_type: str
    dosage: Optional[int]
    status: str
    dominant_class: str
    dominant_fraction: float
    mixed_dosage: bool
    window_count: int
    class_counts: Tuple[Tuple[str, int], ...]


def parse_window_depths(
    lines: Iterable[str],
    *,
    contig_column: int,
    start_column: int,
    end_column: int,
    depth_column: int,
    skip_header: bool = False,
) -> list[WindowDepth]:
    """Parse tabular window depths using zero-based column indexes.

    Blank lines and lines beginning with ``#`` are ignored.  Identifiers and
    coordinates are preserved, while invalid or negative depth is rejected.
    """

    columns = (contig_column, start_column, end_column, depth_column)
    if any(column < 0 for column in columns):
        raise DosageError("column indexes must be non-negative")
    required_columns = max(columns) + 1
    windows: list[WindowDepth] = []
    skipped_header = not skip_header

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        if not skipped_header:
            skipped_header = True
            continue

        fields = line.split("\t")
        if len(fields) < required_columns:
            raise DosageError(
                f"line {line_number}: expected at least {required_columns} tab-separated columns"
            )
        contig_id = fields[contig_column]
        if not contig_id:
            raise DosageError(f"line {line_number}: empty contig identifier")
        try:
            start = int(fields[start_column])
            end = int(fields[end_column])
            depth = float(fields[depth_column])
        except ValueError as exc:
            raise DosageError(f"line {line_number}: invalid coordinate or depth") from exc
        if start < 0 or end < start:
            raise DosageError(
                f"line {line_number}: invalid window coordinates {start}..{end}"
            )
        if not math.isfinite(depth) or depth < 0:
            raise DosageError(f"line {line_number}: depth must be finite and non-negative")
        windows.append(WindowDepth(contig_id, start, end, depth))

    if not windows:
        raise DosageError("no window-depth records were found")
    return windows


def dosage_class_name(dosage: int) -> str:
    """Return the canonical, ploidy-independent name for a copy count."""

    if dosage < 1:
        raise DosageError("dosage must be at least one")
    return f"dosage_{dosage}"


def legacy_contig_type(dosage: int) -> str:
    """Return a legacy label where one exists, otherwise a numeric label."""

    names = {
        1: "haplotig",
        2: "diplotig",
        3: "triplotig",
        4: "tetraplotig",
    }
    return names.get(dosage, dosage_class_name(dosage))


def fit_dosage_model(
    depths: Sequence[float],
    *,
    ploidy: int,
    haploid_depth: Optional[float] = None,
    relative_sigma: Optional[float] = None,
) -> DosageModel:
    """Fit integer-spaced copy-depth states to window depths.

    Automatic fitting minimizes a robust loss to the nearest expected copy
    state (``1 * haploid_depth`` through ``ploidy * haploid_depth``).  Candidate
    values are derived only from deterministic quantiles of the observed
    positive depths.  Zero-depth windows do not inform the fit, but are called
    ``dosage_1`` by the current low-depth classification policy.
    """

    _validate_model_parameters(ploidy, haploid_depth, relative_sigma)
    numeric_depths = [float(depth) for depth in depths]
    if any(not math.isfinite(depth) or depth < 0 for depth in numeric_depths):
        raise DosageError("window depths must be finite and non-negative")
    positive_depths = sorted(depth for depth in numeric_depths if depth > 0)
    if not positive_depths:
        raise DosageError("at least one finite positive window depth is required")
    if len(positive_depths) < 2 and haploid_depth is None:
        raise DosageError(
            "automatic haploid-depth estimation requires at least two positive windows"
        )

    sample = _deterministic_quantile_sample(positive_depths, limit=257)
    alternatives: Tuple[DosageAlternative, ...]
    if haploid_depth is None:
        candidate_scores = _score_haploid_depth_candidates(sample, ploidy)
        best_depth, best_score = candidate_scores[0]
        separated = _separated_alternatives(candidate_scores, best_depth)
        alternatives = tuple(
            DosageAlternative(candidate_depth, score)
            for candidate_depth, score in separated[:5]
        )
        if alternatives:
            second_score = alternatives[0].score
            separation = max(0.0, (second_score - best_score) / max(second_score, 1e-12))
        else:
            separation = 1.0
        estimation_status = "well_separated" if separation >= 0.1 else "weakly_identified"
        fitted_haploid_depth = best_depth
        method = "robust_integer_spacing"
        objective_score = best_score
    else:
        fitted_haploid_depth = float(haploid_depth)
        objective_score = _candidate_score(fitted_haploid_depth, sample, ploidy)
        alternatives = ()
        separation = 1.0
        estimation_status = "user_supplied"
        method = "user_supplied"

    fitted_sigma = (
        float(relative_sigma)
        if relative_sigma is not None
        else _estimate_relative_sigma(sample, fitted_haploid_depth, ploidy)
    )
    return DosageModel(
        ploidy=ploidy,
        haploid_depth=fitted_haploid_depth,
        relative_sigma=fitted_sigma,
        estimation_method=method,
        estimation_status=estimation_status,
        objective_score=objective_score,
        separation=separation,
        fitted_window_count=len(positive_depths),
        alternatives=alternatives,
    )


def classify_window(
    window: WindowDepth,
    model: DosageModel,
    *,
    min_confidence: float = 0.8,
) -> WindowDosageCall:
    """Classify one window as a copy state or an explicit exceptional state."""

    if not 0.5 <= min_confidence < 1.0:
        raise DosageError("min_confidence must be in [0.5, 1.0)")
    ratio = window.depth / model.haploid_depth
    if ratio < 0.5:
        return WindowDosageCall(
            window.contig_id,
            window.start,
            window.end,
            window.depth,
            ratio,
            dosage_class_name(1),
            1,
            1.0,
        )
    if ratio > model.ploidy + 0.5:
        return WindowDosageCall(
            window.contig_id,
            window.start,
            window.end,
            window.depth,
            ratio,
            "high_copy",
            None,
            1.0,
        )

    probabilities = _copy_probabilities(ratio, model.ploidy, model.relative_sigma)
    best_index = max(range(model.ploidy), key=lambda index: probabilities[index])
    dosage = best_index + 1
    confidence = probabilities[best_index]
    if confidence < min_confidence:
        classification = "ambiguous"
        called_dosage: Optional[int] = None
    else:
        classification = dosage_class_name(dosage)
        called_dosage = dosage
    return WindowDosageCall(
        window.contig_id,
        window.start,
        window.end,
        window.depth,
        ratio,
        classification,
        called_dosage,
        confidence,
    )


def classify_windows(
    windows: Sequence[WindowDepth],
    model: DosageModel,
    *,
    min_confidence: float = 0.8,
) -> list[WindowDosageCall]:
    """Classify windows in a stable identifier/coordinate order."""

    ordered = sorted(
        windows,
        key=lambda window: (window.contig_id, window.start, window.end, window.depth),
    )
    return [
        classify_window(window, model, min_confidence=min_confidence) for window in ordered
    ]


def summarize_unitigs(
    calls: Sequence[WindowDosageCall],
    model: DosageModel,
    *,
    min_confidence: float = 0.8,
) -> list[UnitigDosageCall]:
    """Classify each unitig from its mean window depth using ``model``.

    Window-class composition remains visible through ``dominant_class``,
    ``dominant_fraction``, ``class_counts``, and ``mixed_dosage`` for audit
    purposes.  Those fields do not affect the final unitig dosage call.
    """

    by_contig: dict[str, list[WindowDosageCall]] = {}
    for call in calls:
        by_contig.setdefault(call.contig_id, []).append(call)

    summaries: list[UnitigDosageCall] = []
    for contig_id in sorted(by_contig):
        contig_calls = by_contig[contig_id]
        counts = Counter(call.classification for call in contig_calls)
        class_counts = tuple(sorted(counts.items()))
        dominant_class, dominant_count = sorted(
            counts.items(), key=lambda item: (-item[1], item[0])
        )[0]
        dominant_fraction = dominant_count / len(contig_calls)
        confident_dosages = sorted(
            {call.dosage for call in contig_calls if call.dosage is not None}
        )
        mixed_dosage = len(confident_dosages) > 1
        average_depth = math.fsum(call.depth for call in contig_calls) / len(contig_calls)
        average_call = classify_window(
            WindowDepth(contig_id, 0, 0, average_depth),
            model,
            min_confidence=min_confidence,
        )

        if average_call.dosage is not None:
            dosage = average_call.dosage
            status = "assigned"
            contig_type = legacy_contig_type(dosage)
        elif average_call.classification == "high_copy":
            dosage = None
            status = "high_copy"
            contig_type = "high_copy"
        else:
            dosage = None
            status = "ambiguous"
            contig_type = "ambiguous"

        summaries.append(
            UnitigDosageCall(
                contig_id=contig_id,
                average_depth=average_depth,
                contig_type=contig_type,
                dosage=dosage,
                status=status,
                dominant_class=dominant_class,
                dominant_fraction=dominant_fraction,
                mixed_dosage=mixed_dosage,
                window_count=len(contig_calls),
                class_counts=class_counts,
            )
        )
    return summaries


def _validate_model_parameters(
    ploidy: int,
    haploid_depth: Optional[float],
    relative_sigma: Optional[float],
) -> None:
    if ploidy < 1:
        raise DosageError("ploidy must be at least one")
    if haploid_depth is not None and (
        not math.isfinite(haploid_depth) or haploid_depth <= 0
    ):
        raise DosageError("haploid_depth must be finite and positive")
    if relative_sigma is not None and (
        not math.isfinite(relative_sigma) or not 0.0 < relative_sigma <= 0.5
    ):
        raise DosageError("relative_sigma must be in (0, 0.5]")


def _deterministic_quantile_sample(values: Sequence[float], *, limit: int) -> list[float]:
    if len(values) <= limit:
        return list(values)
    return [_quantile(values, index / (limit - 1)) for index in range(limit)]


def _quantile(sorted_values: Sequence[float], fraction: float) -> float:
    position = fraction * (len(sorted_values) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def _score_haploid_depth_candidates(
    sample: Sequence[float], ploidy: int
) -> list[Tuple[float, float]]:
    candidates = {
        round(depth / copy_count, 12)
        for depth in sample
        for copy_count in range(1, ploidy + 1)
        if depth > 0
    }
    scored = [
        (candidate, _candidate_score(candidate, sample, ploidy))
        for candidate in candidates
        if candidate > 0
    ]
    # On an exact tie, prefer the larger fundamental spacing.  Such ties are
    # reported as weakly identified so the user can provide an explicit value.
    scored.sort(key=lambda item: (item[1], -item[0]))
    return scored


def _candidate_score(candidate: float, sample: Sequence[float], ploidy: int) -> float:
    losses = []
    for depth in sample:
        ratio = depth / candidate
        nearest = min(ploidy, max(1, int(math.floor(ratio + 0.5))))
        residual = min(abs(ratio - nearest), 0.5)
        losses.append(residual * residual)
    losses.sort()
    median_loss = _median(losses)
    mean_loss = math.fsum(losses) / len(losses)
    return median_loss + 0.1 * mean_loss


def _separated_alternatives(
    scored: Sequence[Tuple[float, float]], best_depth: float
) -> list[Tuple[float, float]]:
    alternatives: list[Tuple[float, float]] = []
    accepted_depths = [best_depth]
    for candidate_depth, score in scored[1:]:
        if all(
            abs(candidate_depth - accepted) / max(candidate_depth, accepted) >= 0.1
            for accepted in accepted_depths
        ):
            alternatives.append((candidate_depth, score))
            accepted_depths.append(candidate_depth)
    return alternatives


def _estimate_relative_sigma(
    sample: Sequence[float], haploid_depth: float, ploidy: int
) -> float:
    residuals = []
    for depth in sample:
        ratio = depth / haploid_depth
        if 0.5 <= ratio <= ploidy + 0.5:
            nearest = min(ploidy, max(1, int(math.floor(ratio + 0.5))))
            residuals.append(abs(ratio - nearest))
    if not residuals:
        return 0.2
    robust_sigma = 1.4826 * _median(sorted(residuals))
    return min(0.35, max(0.05, robust_sigma))


def _copy_probabilities(ratio: float, ploidy: int, sigma: float) -> list[float]:
    log_weights = [
        -0.5 * ((ratio - copy_count) / sigma) ** 2
        for copy_count in range(1, ploidy + 1)
    ]
    maximum = max(log_weights)
    weights = [math.exp(log_weight - maximum) for log_weight in log_weights]
    total = math.fsum(weights)
    return [weight / total for weight in weights]


def _median(sorted_values: Sequence[float]) -> float:
    midpoint = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return sorted_values[midpoint]
    return (sorted_values[midpoint - 1] + sorted_values[midpoint]) / 2.0

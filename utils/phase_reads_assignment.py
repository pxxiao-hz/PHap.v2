#!/usr/bin/env python3
"""Auditable read assignment and one-pass FASTQ dispatch for PHap v2."""

from __future__ import annotations

import csv
import gzip
import hashlib
import heapq
import itertools
import json
import re
import sqlite3
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Sequence, TextIO, Tuple

import pysam


DOSAGE_BY_TYPE = {
    "haplotig": 1,
    "diplotig": 2,
    "triplotig": 3,
    "tetraplotig": 4,
}
GROUP_PATTERN = re.compile(r"^(.+)_group([1-9][0-9]*)$")


@dataclass(frozen=True)
class GroupModel:
    groups: Dict[str, Tuple[str, ...]]
    unitig_groups: Dict[str, Tuple[str, ...]]
    group_chromosome: Dict[str, str]
    group_number: Dict[str, int]
    ploidy: int


@dataclass(frozen=True)
class AlignmentHit:
    unitig: str
    query_start: int
    query_end: int
    identity: float
    mapq: int
    groups: Tuple[str, ...]

    @property
    def weight(self) -> float:
        mapq_weight = min(max(self.mapq, 0), 60) / 60.0
        return self.identity * mapq_weight


@dataclass
class ReadObservation:
    read_id: str
    read_length: int = 0
    hits: list[AlignmentHit] = field(default_factory=list)


@dataclass(frozen=True)
class EvidenceEvaluation:
    best_group: str
    best_score: float
    second_score: float
    margin: float
    best_aligned_bp: int
    candidate_groups: Tuple[str, ...]
    confident_group: Optional[str]
    unique_group_support: bool
    retained_alignments: int
    supporting_unitigs: Tuple[str, ...]
    reason: str


@dataclass(frozen=True)
class AssignmentDecision:
    read_id: str
    read_length: int
    status: str
    group: Optional[str]
    basis: str
    candidate_groups: Tuple[str, ...]
    best_score: float
    second_score: float
    margin: float
    best_aligned_bp: int
    retained_alignments: int
    supporting_unitigs: Tuple[str, ...]


@dataclass(frozen=True)
class AssignmentParameters:
    min_mapq: int = 20
    min_alignment_length: int = 1000
    min_identity: float = 0.80
    min_total_aligned_bp: int = 1000
    min_group_margin: float = 0.10
    balanced_score_tolerance: float = 0.02
    collapsed_policy: str = "balanced"
    include_supplementary: bool = True


def group_sort_key(group: str):
    match = GROUP_PATTERN.match(group)
    if match is None:
        return group, 0
    return match.group(1), int(match.group(2))


def normalize_read_name(name: str) -> str:
    value = name.strip().split()[0]
    if value.startswith("@") or value.startswith(">"):
        value = value[1:]
    if value.endswith("/1") or value.endswith("/2"):
        value = value[:-2]
    return value


def parse_group_file(path: Path) -> GroupModel:
    groups = {}
    unitig_groups = defaultdict(set)
    chromosomes = defaultdict(set)
    group_chromosome = {}
    group_number = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            group = fields[0]
            match = GROUP_PATTERN.match(group)
            if match is None:
                raise ValueError(
                    f"Invalid final group name at {path}:{line_number}: {group}"
                )
            if group in groups:
                raise ValueError(f"Duplicate group at {path}:{line_number}: {group}")
            members = tuple(fields[1:])
            if len(members) != len(set(members)):
                raise ValueError(f"Duplicate unitig in {group} at {path}:{line_number}")
            chromosome = match.group(1)
            number = int(match.group(2))
            groups[group] = members
            group_chromosome[group] = chromosome
            group_number[group] = number
            chromosomes[chromosome].add(number)
            for unitig in members:
                unitig_groups[unitig].add(group)
    if not groups:
        raise ValueError(f"No groups found in {path}")
    ploidies = {max(numbers) for numbers in chromosomes.values()}
    if len(ploidies) != 1:
        raise ValueError(f"Chromosomes do not have a common ploidy in {path}")
    ploidy = ploidies.pop()
    expected = set(range(1, ploidy + 1))
    for chromosome, numbers in sorted(chromosomes.items()):
        if numbers != expected:
            missing = sorted(expected - numbers)
            raise ValueError(f"{chromosome} is missing haplotype groups: {missing}")
    for unitig, assigned_groups in unitig_groups.items():
        assigned_chromosomes = {
            group_chromosome[group] for group in assigned_groups
        }
        if len(assigned_chromosomes) != 1:
            raise ValueError(
                f"Unitig {unitig!r} occurs in groups from multiple chromosomes: "
                f"{', '.join(sorted(assigned_chromosomes))}"
            )
        if len(assigned_groups) > ploidy:
            raise ValueError(
                f"Unitig {unitig!r} occurs in {len(assigned_groups)} groups, "
                f"exceeding the inferred ploidy ({ploidy})"
            )
    return GroupModel(
        groups=dict(sorted(groups.items(), key=lambda item: group_sort_key(item[0]))),
        unitig_groups={
            unitig: tuple(sorted(values, key=group_sort_key))
            for unitig, values in unitig_groups.items()
        },
        group_chromosome=group_chromosome,
        group_number=group_number,
        ploidy=ploidy,
    )


def parse_contig_types(path: Path) -> Dict[str, str]:
    result = {}
    with path.open() as handle:
        header = handle.readline().split()
        try:
            id_index = header.index("contig_ID")
            type_index = header.index("contig_type")
        except ValueError as exc:
            raise ValueError(
                "contig type file must contain contig_ID and contig_type columns"
            ) from exc
        for line_number, line in enumerate(handle, 2):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= max(id_index, type_index):
                raise ValueError(f"Malformed contig type row at {path}:{line_number}")
            unitig = fields[id_index]
            if unitig in result:
                raise ValueError(f"Duplicate contig type for {unitig} at line {line_number}")
            result[unitig] = fields[type_index]
    return result


def validate_group_dosage(
    model: GroupModel, contig_types: Dict[str, str], unknown_policy: str = "group"
):
    if unknown_policy not in {"group", "error"}:
        raise ValueError("unknown dosage policy must be group or error")
    counts = Counter()
    examples = defaultdict(list)
    for unitig, assigned_groups in model.unitig_groups.items():
        contig_type = contig_types.get(unitig)
        dosage = DOSAGE_BY_TYPE.get(contig_type)
        if dosage is None:
            category = "missing" if contig_type is None else "unsupported"
            counts[category] += 1
            if len(examples[category]) < 5:
                examples[category].append(unitig)
            if unknown_policy == "error":
                raise ValueError(
                    f"Cannot validate dosage for {unitig}: contig_type={contig_type}"
                )
            continue
        if dosage != len(assigned_groups):
            raise ValueError(
                f"Group membership dosage mismatch for {unitig}: "
                f"contig_type={contig_type} expects {dosage}, observed {len(assigned_groups)}"
            )
        counts["validated"] += 1
    return {
        "unitigs": len(model.unitig_groups),
        "validated": counts["validated"],
        "missing_contig_type": counts["missing"],
        "unsupported_contig_type": counts["unsupported"],
        "examples": dict(sorted(examples.items())),
        "unknown_policy": unknown_policy,
    }


def _alignment_identity(alignment) -> float:
    aligned = max(alignment.query_alignment_length or 0, 1)
    if alignment.has_tag("NM"):
        return max(0.0, 1.0 - alignment.get_tag("NM") / aligned)
    return 1.0


def _alignment_filter_reason(alignment, parameters: AssignmentParameters):
    if alignment.is_unmapped:
        return "unmapped"
    if alignment.is_secondary:
        return "secondary"
    if alignment.is_duplicate:
        return "duplicate"
    if alignment.is_qcfail:
        return "qcfail"
    if alignment.is_supplementary and not parameters.include_supplementary:
        return "supplementary"
    if alignment.mapping_quality < parameters.min_mapq:
        return "low_mapq"
    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    if query_start is None or query_end is None:
        return "missing_query_interval"
    if query_end - query_start < parameters.min_alignment_length:
        return "short_alignment"
    if _alignment_identity(alignment) < parameters.min_identity:
        return "low_identity"
    return None


def _alignment_is_usable(alignment, parameters: AssignmentParameters) -> bool:
    return _alignment_filter_reason(alignment, parameters) is None


def _alignment_hit(alignment, bam, model: GroupModel) -> Optional[AlignmentHit]:
    if alignment.reference_id < 0:
        return None
    unitig = bam.get_reference_name(alignment.reference_id)
    groups = model.unitig_groups.get(unitig)
    if not groups:
        return None
    return AlignmentHit(
        unitig=unitig,
        query_start=alignment.query_alignment_start,
        query_end=alignment.query_alignment_end,
        identity=_alignment_identity(alignment),
        mapq=alignment.mapping_quality,
        groups=groups,
    )


def collect_long_read_evidence(
    bam_path: Path, model: GroupModel, parameters: AssignmentParameters
):
    observations = {}
    stats = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reference_lengths = dict(zip(bam.references, bam.lengths))
        for alignment in bam.fetch(until_eof=True):
            stats["alignment_records"] += 1
            if not _alignment_is_usable(alignment, parameters):
                stats["filtered_alignment_records"] += 1
                continue
            hit = _alignment_hit(alignment, bam, model)
            if hit is None:
                stats["non_group_alignment_records"] += 1
                continue
            read_id = normalize_read_name(alignment.query_name)
            read_length = alignment.infer_read_length() or alignment.query_length or 0
            observation = observations.setdefault(read_id, ReadObservation(read_id))
            observation.read_length = max(observation.read_length, read_length)
            observation.hits.append(hit)
            stats["retained_alignment_records"] += 1
    stats["reads_with_group_evidence"] = len(observations)
    return observations, reference_lengths, dict(stats)


def collect_hic_pair_evidence(
    bam_path: Path, model: GroupModel, parameters: AssignmentParameters
):
    pairs = defaultdict(lambda: {1: None, 2: None})
    stats = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reference_lengths = dict(zip(bam.references, bam.lengths))
        for alignment in bam.fetch(until_eof=True):
            stats["alignment_records"] += 1
            if not (alignment.is_read1 or alignment.is_read2):
                stats["non_pair_alignment_records"] += 1
                continue
            if not _alignment_is_usable(alignment, parameters):
                stats["filtered_alignment_records"] += 1
                continue
            hit = _alignment_hit(alignment, bam, model)
            if hit is None:
                stats["non_group_alignment_records"] += 1
                continue
            read_id = normalize_read_name(alignment.query_name)
            mate = 1 if alignment.is_read1 else 2
            observation = pairs[read_id][mate]
            if observation is None:
                observation = ReadObservation(read_id)
                pairs[read_id][mate] = observation
            read_length = alignment.infer_read_length() or alignment.query_length or 0
            observation.read_length = max(observation.read_length, read_length)
            observation.hits.append(hit)
            stats["retained_alignment_records"] += 1
    stats["pairs_with_group_evidence"] = len(pairs)
    return dict(pairs), reference_lengths, dict(stats)


def group_target_lengths(model: GroupModel, reference_lengths: Dict[str, int]):
    missing = sorted(set(model.unitig_groups) - set(reference_lengths))
    if missing:
        raise ValueError(
            f"{len(missing)} grouped unitigs are absent from BAM references: "
            f"{', '.join(missing[:5])}"
        )
    return {
        group: sum(reference_lengths[unitig] for unitig in members)
        for group, members in model.groups.items()
    }


def _interval_score(intervals: Sequence[Tuple[int, int, float]]) -> Tuple[float, int]:
    if not intervals:
        return 0.0, 0
    events = defaultdict(list)
    for start, end, weight in intervals:
        events[start].append((1, weight))
        events[end].append((-1, weight))
    active = Counter()
    heap = []
    score = 0.0
    aligned_bp = 0
    previous = None
    for position in sorted(events):
        while heap and active[-heap[0]] == 0:
            heapq.heappop(heap)
        if previous is not None and heap:
            length = position - previous
            aligned_bp += length
            score += length * (-heap[0])
        for direction, weight in events[position]:
            if direction > 0:
                active[weight] += 1
                heapq.heappush(heap, -weight)
            else:
                active[weight] -= 1
        previous = position
    return score, aligned_bp


def evaluate_observation(
    observation: Optional[ReadObservation],
    model: GroupModel,
    parameters: AssignmentParameters,
) -> Optional[EvidenceEvaluation]:
    if observation is None or not observation.hits:
        return None
    intervals = defaultdict(list)
    for hit in observation.hits:
        for group in hit.groups:
            intervals[group].append((hit.query_start, hit.query_end, hit.weight))
    scored = []
    for group, values in intervals.items():
        score, aligned_bp = _interval_score(values)
        scored.append((score, aligned_bp, group))
    scored.sort(key=lambda item: (-item[0], -item[1], group_sort_key(item[2])))
    best_score, best_aligned_bp, best_group = scored[0]
    second_score = scored[1][0] if len(scored) > 1 else 0.0
    margin = (best_score - second_score) / best_score if best_score > 0 else 0.0
    candidate_groups = tuple(
        sorted(
            (
                group
                for score, _, group in scored
                if best_score > 0
                and (best_score - score) / best_score
                <= parameters.balanced_score_tolerance
            ),
            key=group_sort_key,
        )
    )
    unique_support = any(len(hit.groups) == 1 for hit in observation.hits)
    confident_group = None
    if best_aligned_bp < parameters.min_total_aligned_bp:
        reason = "insufficient_aligned_bp"
    elif margin >= parameters.min_group_margin:
        confident_group = best_group
        reason = "alignment_confident"
    elif (
        not unique_support
        and len(candidate_groups) > 1
        and len({model.group_chromosome[group] for group in candidate_groups}) == 1
    ):
        reason = "collapsed_balanced_candidate"
    elif unique_support:
        reason = "ambiguous_unique_group_evidence"
    else:
        reason = "ambiguous_or_cross_chromosome"
    return EvidenceEvaluation(
        best_group=best_group,
        best_score=best_score,
        second_score=second_score,
        margin=margin,
        best_aligned_bp=best_aligned_bp,
        candidate_groups=candidate_groups,
        confident_group=confident_group,
        unique_group_support=unique_support,
        retained_alignments=len(observation.hits),
        supporting_unitigs=tuple(sorted({hit.unitig for hit in observation.hits})),
        reason=reason,
    )


def _decision(
    observation: ReadObservation,
    evaluation: EvidenceEvaluation,
    status: str,
    group: Optional[str],
    basis: str,
    candidate_groups: Optional[Tuple[str, ...]] = None,
):
    return AssignmentDecision(
        read_id=observation.read_id,
        read_length=observation.read_length,
        status=status,
        group=group,
        basis=basis,
        candidate_groups=(
            evaluation.candidate_groups
            if candidate_groups is None
            else candidate_groups
        ),
        best_score=evaluation.best_score,
        second_score=evaluation.second_score,
        margin=evaluation.margin,
        best_aligned_bp=evaluation.best_aligned_bp,
        retained_alignments=evaluation.retained_alignments,
        supporting_unitigs=evaluation.supporting_unitigs,
    )


def _stable_hash(seed: int, read_id: str) -> str:
    return hashlib.sha256(f"{seed}\0{read_id}".encode()).hexdigest()


def _balanced_assign(
    pending,
    decisions,
    load,
    target_lengths,
    seed,
    weight_getter,
    decision_builder,
):
    for item in sorted(
        pending,
        key=lambda value: (
            -weight_getter(value),
            _stable_hash(seed, value[0]),
            value[0],
        ),
    ):
        read_id, candidates = item[0], item[1]
        chosen = min(
            candidates,
            key=lambda group: (
                load[group] / max(target_lengths[group], 1),
                group_sort_key(group),
            ),
        )
        decisions[read_id] = decision_builder(item, chosen)
        load[chosen] += weight_getter(item)


def assign_long_reads(
    observations: Dict[str, ReadObservation],
    model: GroupModel,
    target_lengths: Dict[str, int],
    parameters: AssignmentParameters,
    seed: int,
):
    decisions = {}
    load = Counter()
    pending = []
    for read_id in sorted(observations):
        observation = observations[read_id]
        evaluation = evaluate_observation(observation, model, parameters)
        if evaluation is None:
            continue
        if evaluation.confident_group is not None:
            group = evaluation.confident_group
            decisions[read_id] = _decision(
                observation, evaluation, "assigned", group, "alignment_confident"
            )
            load[group] += observation.read_length
        elif (
            evaluation.reason == "collapsed_balanced_candidate"
            and parameters.collapsed_policy == "balanced"
        ):
            pending.append((read_id, evaluation.candidate_groups, observation, evaluation))
        else:
            decisions[read_id] = _decision(
                observation, evaluation, "unassigned", None, evaluation.reason
            )

    def build(item, chosen):
        _, candidates, observation, evaluation = item
        return _decision(
            observation,
            evaluation,
            "assigned",
            chosen,
            "collapsed_depth_balanced",
            candidates,
        )

    _balanced_assign(
        pending,
        decisions,
        load,
        target_lengths,
        seed,
        lambda item: item[2].read_length,
        build,
    )
    return decisions, assignment_summary(decisions, target_lengths)


def assign_hic_pairs(
    pairs,
    model: GroupModel,
    target_lengths: Dict[str, int],
    parameters: AssignmentParameters,
    seed: int,
):
    decisions = {}
    load = Counter()
    pending = []
    for read_id in sorted(pairs):
        mate1 = pairs[read_id].get(1)
        mate2 = pairs[read_id].get(2)
        eval1 = evaluate_observation(mate1, model, parameters)
        eval2 = evaluate_observation(mate2, model, parameters)
        read_length = sum(
            observation.read_length for observation in (mate1, mate2) if observation
        )
        combined = ReadObservation(read_id, read_length)
        for observation in (mate1, mate2):
            if observation:
                combined.hits.extend(observation.hits)
        combined_eval = evaluate_observation(combined, model, parameters)
        if combined_eval is None:
            continue
        group = None
        basis = None
        candidates = ()
        if eval1 is None or eval2 is None:
            basis = "one_mate_without_group_evidence"
        elif eval1.confident_group and eval2.confident_group:
            if eval1.confident_group == eval2.confident_group:
                group = eval1.confident_group
                basis = "hic_mates_concordant"
            else:
                basis = "hic_mates_conflict"
        elif eval1.confident_group or eval2.confident_group:
            anchor = eval1.confident_group or eval2.confident_group
            other = eval2 if eval1.confident_group else eval1
            if anchor in other.candidate_groups:
                group = anchor
                basis = "hic_anchor_compatible_mate"
            else:
                basis = "hic_anchor_incompatible_mate"
        else:
            candidates = tuple(
                sorted(
                    set(eval1.candidate_groups) & set(eval2.candidate_groups),
                    key=group_sort_key,
                )
            )
            if (
                parameters.collapsed_policy == "balanced"
                and not eval1.unique_group_support
                and not eval2.unique_group_support
                and len(candidates) > 1
                and len({model.group_chromosome[value] for value in candidates}) == 1
            ):
                pending.append(
                    (read_id, candidates, combined, combined_eval, eval1, eval2)
                )
                continue
            basis = "hic_pair_ambiguous"
        if group is not None:
            decisions[read_id] = _decision(
                combined, combined_eval, "assigned", group, basis
            )
            load[group] += 1
        else:
            decisions[read_id] = _decision(
                combined,
                combined_eval,
                "unassigned",
                None,
                basis,
                candidates or combined_eval.candidate_groups,
            )

    def build(item, chosen):
        _, candidates, combined, combined_eval, _, _ = item
        return _decision(
            combined,
            combined_eval,
            "assigned",
            chosen,
            "hic_collapsed_pair_balanced",
            candidates,
        )

    _balanced_assign(
        pending,
        decisions,
        load,
        target_lengths,
        seed,
        lambda item: 1,
        build,
    )
    return decisions, assignment_summary(decisions, target_lengths)


def assignment_summary(decisions, target_lengths):
    basis_counts = Counter()
    group_counts = Counter()
    group_bp = Counter()
    assigned = 0
    unassigned = 0
    for decision in decisions.values():
        basis_counts[decision.basis] += 1
        if decision.status == "assigned":
            assigned += 1
            group_counts[decision.group] += 1
            group_bp[decision.group] += decision.read_length
        else:
            unassigned += 1
    groups = {}
    for group in sorted(target_lengths, key=group_sort_key):
        target_bp = target_lengths[group]
        groups[group] = {
            "reads": group_counts[group],
            "read_bp": group_bp[group],
            "target_bp": target_bp,
            "estimated_depth": group_bp[group] / max(target_bp, 1),
        }
    return {
        "reads_with_group_evidence": len(decisions),
        "assigned": assigned,
        "unassigned": unassigned,
        "basis_counts": dict(sorted(basis_counts.items())),
        "groups": groups,
    }


ASSIGNMENT_FIELDS = (
    "read_id",
    "read_length",
    "status",
    "group",
    "basis",
    "candidate_groups",
    "best_score",
    "second_score",
    "margin",
    "best_aligned_bp",
    "retained_alignments",
    "supporting_unitigs",
)


def _sqlite_connect(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path), timeout=120)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-65536")
    return connection


def _meta_get(connection, key, default=None):
    row = connection.execute(
        "SELECT value FROM metadata WHERE key = ?", (key,)
    ).fetchone()
    return default if row is None else json.loads(row[0])


def _meta_set(connection, key, value):
    connection.execute(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        (key, json.dumps(value, sort_keys=True)),
    )


def _reset_evidence_database(connection, checkpoint_key):
    connection.executescript(
        """
        DROP TABLE IF EXISTS evidence;
        DROP TABLE IF EXISTS metadata;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE evidence (
            read_id TEXT NOT NULL,
            mate INTEGER NOT NULL,
            read_length INTEGER NOT NULL,
            unitig TEXT NOT NULL,
            query_start INTEGER NOT NULL,
            query_end INTEGER NOT NULL,
            identity REAL NOT NULL,
            mapq INTEGER NOT NULL
        );
        """
    )
    _meta_set(connection, "checkpoint_key", checkpoint_key)
    _meta_set(connection, "last_bam_offset", 0)
    _meta_set(connection, "stats", {})
    _meta_set(connection, "collection_complete", False)
    connection.commit()


def collect_evidence_disk(
    bam_path: Path,
    model: GroupModel,
    parameters: AssignmentParameters,
    database_path: Path,
    checkpoint_key: str,
    hic: bool = False,
    progress=None,
    progress_every: int = 1_000_000,
):
    """Stream a BAM into a resumable, read-keyed on-disk evidence store."""
    connection = _sqlite_connect(database_path)
    connection.execute(
        "CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
    )
    if _meta_get(connection, "checkpoint_key") != checkpoint_key:
        _reset_evidence_database(connection, checkpoint_key)
    if _meta_get(connection, "collection_complete", False):
        result = {
            "reference_lengths": _meta_get(connection, "reference_lengths", {}),
            "stats": _meta_get(connection, "stats", {}),
            "resumed_from_record": _meta_get(connection, "stats", {}).get(
                "alignment_records", 0
            ),
        }
        connection.close()
        return result

    stats = Counter(_meta_get(connection, "stats", {}))
    resumed_from_record = stats["alignment_records"]
    last_offset = int(_meta_get(connection, "last_bam_offset", 0))
    checkpoint_interval = min(max(progress_every, 1), 50_000)
    bam_size = bam_path.stat().st_size
    started = time.monotonic()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reference_lengths = dict(zip(bam.references, bam.lengths))
        previous_lengths = _meta_get(connection, "reference_lengths")
        if previous_lengths is not None and previous_lengths != reference_lengths:
            connection.close()
            raise ValueError(f"BAM header changed while resuming {bam_path}")
        _meta_set(connection, "reference_lengths", reference_lengths)
        if last_offset:
            bam.seek(last_offset)
        iterator = bam.fetch(until_eof=True)
        pending_rows = []
        since_checkpoint = 0
        for alignment in iterator:
            stats["alignment_records"] += 1
            since_checkpoint += 1
            if hic and not (alignment.is_read1 or alignment.is_read2):
                stats["filtered_non_pair"] += 1
            else:
                reason = _alignment_filter_reason(alignment, parameters)
                if reason is not None:
                    stats[f"filtered_{reason}"] += 1
                else:
                    hit = _alignment_hit(alignment, bam, model)
                    if hit is None:
                        stats["non_group_alignment_records"] += 1
                    else:
                        read_id = normalize_read_name(alignment.query_name)
                        read_length = (
                            alignment.infer_read_length()
                            or alignment.query_length
                            or 0
                        )
                        pending_rows.append(
                            (
                                read_id,
                                1 if hic and alignment.is_read1 else (2 if hic else 0),
                                read_length,
                                hit.unitig,
                                hit.query_start,
                                hit.query_end,
                                hit.identity,
                                hit.mapq,
                            )
                        )
                        stats["retained_alignment_records"] += 1
            if since_checkpoint >= checkpoint_interval:
                connection.executemany(
                    "INSERT INTO evidence VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    pending_rows,
                )
                pending_rows.clear()
                last_offset = bam.tell()
                _meta_set(connection, "last_bam_offset", last_offset)
                _meta_set(connection, "stats", dict(stats))
                connection.commit()
                if progress:
                    byte_percent = 100.0 * (last_offset >> 16) / max(bam_size, 1)
                    progress(
                        stats["alignment_records"], None, started,
                        f"BAM bytes {min(byte_percent, 100.0):.1f}%",
                    )
                since_checkpoint = 0
        if pending_rows:
            connection.executemany(
                "INSERT INTO evidence VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                pending_rows,
            )
        _meta_set(connection, "last_bam_offset", bam.tell())
        _meta_set(connection, "stats", dict(stats))
        connection.commit()
    connection.execute(
        "CREATE INDEX IF NOT EXISTS evidence_read_mate "
        "ON evidence(read_id, mate)"
    )
    _meta_set(connection, "collection_complete", True)
    connection.commit()
    if progress:
        progress(
            stats["alignment_records"], stats["alignment_records"], started,
            "BAM bytes 100.0%",
        )
    result = {
        "reference_lengths": reference_lengths,
        "stats": dict(stats),
        "resumed_from_record": resumed_from_record,
    }
    connection.close()
    return result


def _hit_from_disk_row(row, model):
    unitig = row[3]
    return AlignmentHit(
        unitig=unitig,
        query_start=row[4],
        query_end=row[5],
        identity=row[6],
        mapq=row[7],
        groups=model.unitig_groups[unitig],
    )


def iter_long_observations_disk(database_path: Path, model: GroupModel):
    connection = _sqlite_connect(database_path)
    cursor = connection.execute(
        "SELECT read_id, mate, read_length, unitig, query_start, query_end, "
        "identity, mapq FROM evidence ORDER BY read_id, mate"
    )
    current = None
    observation = None
    try:
        for row in cursor:
            if row[0] != current:
                if observation is not None:
                    yield observation
                current = row[0]
                observation = ReadObservation(current)
            observation.read_length = max(observation.read_length, row[2])
            observation.hits.append(_hit_from_disk_row(row, model))
        if observation is not None:
            yield observation
    finally:
        connection.close()


def iter_hic_pairs_disk(database_path: Path, model: GroupModel):
    connection = _sqlite_connect(database_path)
    cursor = connection.execute(
        "SELECT read_id, mate, read_length, unitig, query_start, query_end, "
        "identity, mapq FROM evidence ORDER BY read_id, mate"
    )
    current = None
    pair = None
    try:
        for row in cursor:
            if row[0] != current:
                if pair is not None:
                    yield current, pair
                current = row[0]
                pair = {1: None, 2: None}
            mate = row[1]
            if pair[mate] is None:
                pair[mate] = ReadObservation(current)
            pair[mate].read_length = max(pair[mate].read_length, row[2])
            pair[mate].hits.append(_hit_from_disk_row(row, model))
        if pair is not None:
            yield current, pair
    finally:
        connection.close()


_DECISION_SQL_FIELDS = (
    "read_id", "read_length", "status", "group_name", "basis",
    "candidate_groups", "best_score", "second_score", "margin",
    "best_aligned_bp", "retained_alignments", "supporting_unitigs",
    "stable_hash",
)


def _reset_decision_database(path: Path):
    if path.exists():
        path.unlink()
    connection = _sqlite_connect(path)
    schema = """
        read_id TEXT PRIMARY KEY,
        read_length INTEGER NOT NULL,
        status TEXT NOT NULL,
        group_name TEXT,
        basis TEXT NOT NULL,
        candidate_groups TEXT NOT NULL,
        best_score REAL NOT NULL,
        second_score REAL NOT NULL,
        margin REAL NOT NULL,
        best_aligned_bp INTEGER NOT NULL,
        retained_alignments INTEGER NOT NULL,
        supporting_unitigs TEXT NOT NULL,
        stable_hash TEXT NOT NULL
    """
    connection.execute(f"CREATE TABLE decisions ({schema})")
    connection.execute(f"CREATE TABLE pending ({schema})")
    connection.execute("CREATE INDEX decisions_status_group ON decisions(status, group_name)")
    return connection


def _decision_values(decision: AssignmentDecision, seed: int):
    return (
        decision.read_id,
        decision.read_length,
        decision.status,
        decision.group,
        decision.basis,
        ",".join(decision.candidate_groups),
        decision.best_score,
        decision.second_score,
        decision.margin,
        decision.best_aligned_bp,
        decision.retained_alignments,
        ",".join(decision.supporting_unitigs),
        _stable_hash(seed, decision.read_id),
    )


def _insert_disk_decision(connection, table, decision, seed):
    placeholders = ",".join("?" for _ in _DECISION_SQL_FIELDS)
    connection.execute(
        f"INSERT OR REPLACE INTO {table} ({','.join(_DECISION_SQL_FIELDS)}) "
        f"VALUES ({placeholders})",
        _decision_values(decision, seed),
    )


def _pending_row_to_decision(row, chosen, basis):
    return AssignmentDecision(
        read_id=row[0], read_length=row[1], status="assigned", group=chosen,
        basis=basis,
        candidate_groups=tuple(filter(None, row[5].split(","))),
        best_score=row[6], second_score=row[7], margin=row[8],
        best_aligned_bp=row[9], retained_alignments=row[10],
        supporting_unitigs=tuple(filter(None, row[11].split(","))),
    )


def _finish_disk_balancing(
    connection, load, target_lengths, seed, weight_mode, basis, progress=None
):
    total = connection.execute("SELECT COUNT(*) FROM pending").fetchone()[0]
    cursor = connection.execute(
        f"SELECT {','.join(_DECISION_SQL_FIELDS)} FROM pending "
        "ORDER BY read_length DESC, stable_hash, read_id"
    )
    started = time.monotonic()
    processed = 0
    for row in cursor:
        candidates = tuple(filter(None, row[5].split(",")))
        chosen = min(
            candidates,
            key=lambda group: (
                load[group] / max(target_lengths[group], 1), group_sort_key(group)
            ),
        )
        decision = _pending_row_to_decision(row, chosen, basis)
        _insert_disk_decision(connection, "decisions", decision, seed)
        load[chosen] += row[1] if weight_mode == "bp" else 1
        processed += 1
        if processed % 10000 == 0:
            connection.commit()
            if progress:
                progress(processed, total, started)
    connection.commit()
    connection.execute("DELETE FROM pending")
    connection.commit()
    if progress and total:
        progress(total, total, started)


def _disk_assignment_summary(connection, target_lengths):
    basis_counts = dict(
        connection.execute(
            "SELECT basis, COUNT(*) FROM decisions GROUP BY basis ORDER BY basis"
        ).fetchall()
    )
    status_counts = dict(
        connection.execute(
            "SELECT status, COUNT(*) FROM decisions GROUP BY status"
        ).fetchall()
    )
    group_values = {
        row[0]: (row[1], row[2])
        for row in connection.execute(
            "SELECT group_name, COUNT(*), COALESCE(SUM(read_length), 0) "
            "FROM decisions WHERE status='assigned' GROUP BY group_name"
        )
    }
    groups = {}
    for group in sorted(target_lengths, key=group_sort_key):
        reads, read_bp = group_values.get(group, (0, 0))
        groups[group] = {
            "reads": reads,
            "read_bp": read_bp,
            "target_bp": target_lengths[group],
            "estimated_depth": read_bp / max(target_lengths[group], 1),
        }
    return {
        "reads_with_group_evidence": sum(status_counts.values()),
        "assigned": status_counts.get("assigned", 0),
        "unassigned": status_counts.get("unassigned", 0),
        "basis_counts": basis_counts,
        "groups": groups,
    }


def write_assignments_database(path: Path, database_path: Path):
    connection = _sqlite_connect(database_path)
    with _open_text(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in connection.execute(
            f"SELECT {','.join(_DECISION_SQL_FIELDS[:-1])} "
            "FROM decisions ORDER BY read_id"
        ):
            writer.writerow(
                {
                    "read_id": row[0], "read_length": row[1], "status": row[2],
                    "group": row[3] or "", "basis": row[4],
                    "candidate_groups": row[5], "best_score": f"{row[6]:.6f}",
                    "second_score": f"{row[7]:.6f}", "margin": f"{row[8]:.6f}",
                    "best_aligned_bp": row[9], "retained_alignments": row[10],
                    "supporting_unitigs": row[11],
                }
            )
    connection.close()


def assign_long_reads_disk(
    evidence_database: Path, decision_database: Path, model: GroupModel,
    target_lengths, parameters, seed, progress=None,
):
    connection = _reset_decision_database(decision_database)
    load = Counter()
    started = time.monotonic()
    processed = 0
    for observation in iter_long_observations_disk(evidence_database, model):
        evaluation = evaluate_observation(observation, model, parameters)
        if evaluation.confident_group is not None:
            decision = _decision(
                observation, evaluation, "assigned", evaluation.confident_group,
                "alignment_confident",
            )
            load[evaluation.confident_group] += observation.read_length
            _insert_disk_decision(connection, "decisions", decision, seed)
        elif (
            evaluation.reason == "collapsed_balanced_candidate"
            and parameters.collapsed_policy == "balanced"
        ):
            decision = _decision(
                observation, evaluation, "pending", None,
                "collapsed_balanced_candidate",
            )
            _insert_disk_decision(connection, "pending", decision, seed)
        else:
            decision = _decision(
                observation, evaluation, "unassigned", None, evaluation.reason
            )
            _insert_disk_decision(connection, "decisions", decision, seed)
        processed += 1
        if processed % 10000 == 0:
            connection.commit()
            if progress:
                progress(processed, None, started)
    connection.commit()
    _finish_disk_balancing(
        connection, load, target_lengths, seed, "bp", "collapsed_depth_balanced",
        progress,
    )
    summary = _disk_assignment_summary(connection, target_lengths)
    connection.close()
    return summary


def assign_hic_pairs_disk(
    evidence_database: Path, decision_database: Path, model: GroupModel,
    target_lengths, parameters, seed, progress=None,
):
    connection = _reset_decision_database(decision_database)
    load = Counter()
    started = time.monotonic()
    processed = 0
    for read_id, pair in iter_hic_pairs_disk(evidence_database, model):
        mate1, mate2 = pair.get(1), pair.get(2)
        eval1 = evaluate_observation(mate1, model, parameters)
        eval2 = evaluate_observation(mate2, model, parameters)
        read_length = sum(
            value.read_length for value in (mate1, mate2) if value is not None
        )
        combined = ReadObservation(read_id, read_length)
        for value in (mate1, mate2):
            if value is not None:
                combined.hits.extend(value.hits)
        combined_eval = evaluate_observation(combined, model, parameters)
        group = None
        basis = None
        candidates = ()
        pending = False
        if eval1 is None or eval2 is None:
            basis = "one_mate_without_group_evidence"
        elif eval1.confident_group and eval2.confident_group:
            if eval1.confident_group == eval2.confident_group:
                group, basis = eval1.confident_group, "hic_mates_concordant"
            else:
                basis = "hic_mates_conflict"
        elif eval1.confident_group or eval2.confident_group:
            anchor = eval1.confident_group or eval2.confident_group
            other = eval2 if eval1.confident_group else eval1
            if anchor in other.candidate_groups:
                group, basis = anchor, "hic_anchor_compatible_mate"
            else:
                basis = "hic_anchor_incompatible_mate"
        else:
            candidates = tuple(
                sorted(
                    set(eval1.candidate_groups) & set(eval2.candidate_groups),
                    key=group_sort_key,
                )
            )
            pending = (
                parameters.collapsed_policy == "balanced"
                and not eval1.unique_group_support
                and not eval2.unique_group_support
                and len(candidates) > 1
                and len({model.group_chromosome[value] for value in candidates}) == 1
            )
            basis = "hic_collapsed_balanced_candidate" if pending else "hic_pair_ambiguous"
        if group is not None:
            decision = _decision(combined, combined_eval, "assigned", group, basis)
            load[group] += 1
            _insert_disk_decision(connection, "decisions", decision, seed)
        elif pending:
            decision = _decision(
                combined, combined_eval, "pending", None, basis, candidates
            )
            _insert_disk_decision(connection, "pending", decision, seed)
        else:
            decision = _decision(
                combined, combined_eval, "unassigned", None, basis,
                candidates or combined_eval.candidate_groups,
            )
            _insert_disk_decision(connection, "decisions", decision, seed)
        processed += 1
        if processed % 10000 == 0:
            connection.commit()
            if progress:
                progress(processed, None, started)
    connection.commit()
    _finish_disk_balancing(
        connection, load, target_lengths, seed, "pairs",
        "hic_collapsed_pair_balanced", progress,
    )
    summary = _disk_assignment_summary(connection, target_lengths)
    connection.close()
    return summary


def _open_text(path: Path, mode: str, compresslevel: int = 1):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", compresslevel=compresslevel)
    return path.open(mode)


def write_assignments(path: Path, decisions):
    path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t")
        writer.writeheader()
        for read_id in sorted(decisions):
            decision = decisions[read_id]
            writer.writerow(
                {
                    "read_id": decision.read_id,
                    "read_length": decision.read_length,
                    "status": decision.status,
                    "group": decision.group or "",
                    "basis": decision.basis,
                    "candidate_groups": ",".join(decision.candidate_groups),
                    "best_score": f"{decision.best_score:.6f}",
                    "second_score": f"{decision.second_score:.6f}",
                    "margin": f"{decision.margin:.6f}",
                    "best_aligned_bp": decision.best_aligned_bp,
                    "retained_alignments": decision.retained_alignments,
                    "supporting_unitigs": ",".join(decision.supporting_unitigs),
                }
            )


def read_assignments(path: Path):
    decisions = {}
    with _open_text(path, "r") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decision = AssignmentDecision(
                read_id=row["read_id"],
                read_length=int(row["read_length"]),
                status=row["status"],
                group=row["group"] or None,
                basis=row["basis"],
                candidate_groups=tuple(filter(None, row["candidate_groups"].split(","))),
                best_score=float(row["best_score"]),
                second_score=float(row["second_score"]),
                margin=float(row["margin"]),
                best_aligned_bp=int(row["best_aligned_bp"]),
                retained_alignments=int(row["retained_alignments"]),
                supporting_unitigs=tuple(filter(None, row["supporting_unitigs"].split(","))),
            )
            if decision.read_id in decisions:
                raise ValueError(f"Duplicate read assignment for {decision.read_id}")
            decisions[decision.read_id] = decision
    return decisions


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


@dataclass(frozen=True)
class FastqRecord:
    name: str
    text: str
    sequence: str
    quality: str


def iter_fastq(path: Path) -> Iterator[FastqRecord]:
    with _open_text(path, "r") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}: {header.strip()}")
            yield FastqRecord(
                normalize_read_name(header),
                header + sequence + plus + quality,
                sequence.strip(),
                quality.strip(),
            )


def _fastq_chunks(records, max_records=2000, max_bases=20_000_000):
    chunk = []
    bases = 0
    for record in records:
        chunk.append(record)
        if isinstance(record, tuple):
            bases += sum(
                len(value.sequence) for value in record if value is not None
            )
        else:
            bases += len(record.sequence)
        if len(chunk) >= max_records or bases >= max_bases:
            yield chunk
            chunk = []
            bases = 0
    if chunk:
        yield chunk


def _lookup_assigned_groups(connection, read_ids, allowed_groups):
    result = {}
    allowed = set(allowed_groups)
    unique_ids = list(dict.fromkeys(read_ids))
    for start in range(0, len(unique_ids), 800):
        values = unique_ids[start:start + 800]
        placeholders = ",".join("?" for _ in values)
        for read_id, group in connection.execute(
            f"SELECT read_id, group_name FROM decisions "
            f"WHERE status='assigned' AND read_id IN ({placeholders})",
            values,
        ):
            if group in allowed:
                result[read_id] = group
    return result


def _prepare_extraction_audit(connection):
    connection.execute("DROP TABLE IF EXISTS extraction_found")
    connection.execute(
        "CREATE TABLE extraction_found ("
        "read_id TEXT PRIMARY KEY, disposition TEXT NOT NULL)"
    )
    connection.commit()


def _requested_assignment_count(connection, groups):
    groups = tuple(groups)
    placeholders = ",".join("?" for _ in groups)
    return connection.execute(
        f"SELECT COUNT(*) FROM decisions WHERE status='assigned' "
        f"AND group_name IN ({placeholders})",
        groups,
    ).fetchone()[0]


def _missing_assignment_values(connection, groups):
    groups = tuple(groups)
    placeholders = ",".join("?" for _ in groups)
    arguments = list(groups)
    missing = connection.execute(
        f"SELECT COUNT(*) FROM decisions d LEFT JOIN extraction_found f "
        f"ON d.read_id=f.read_id WHERE d.status='assigned' "
        f"AND d.group_name IN ({placeholders}) AND f.read_id IS NULL",
        arguments,
    ).fetchone()[0]
    examples = [
        row[0]
        for row in connection.execute(
            f"SELECT d.read_id FROM decisions d LEFT JOIN extraction_found f "
            f"ON d.read_id=f.read_id WHERE d.status='assigned' "
            f"AND d.group_name IN ({placeholders}) AND f.read_id IS NULL "
            f"ORDER BY d.read_id LIMIT 10",
            arguments,
        )
    ]
    return missing, examples


def dispatch_single_fastq_database(
    input_path: Path,
    assignment_database: Path,
    groups,
    output_directory: Path,
    suffix: str,
    min_length: int = 0,
    min_mean_quality: float = 0.0,
    progress=None,
):
    writers = _group_fastq_writers(output_directory, groups, suffix)
    connection = _sqlite_connect(assignment_database)
    _prepare_extraction_audit(connection)
    requested = _requested_assignment_count(connection, groups)
    records = written = filtered = 0
    group_written = Counter()
    started = time.monotonic()
    try:
        for chunk in _fastq_chunks(iter_fastq(input_path)):
            assignments = _lookup_assigned_groups(
                connection, [record.name for record in chunk], groups
            )
            for record in chunk:
                records += 1
                group = assignments.get(record.name)
                if group is None:
                    continue
                mean_quality = (
                    sum(ord(value) - 33 for value in record.quality)
                    / len(record.quality)
                    if record.quality else 0.0
                )
                disposition = "written"
                if len(record.sequence) < min_length or mean_quality < min_mean_quality:
                    filtered += 1
                    disposition = "filtered"
                else:
                    writers[group].write(record.text)
                    written += 1
                    group_written[group] += 1
                try:
                    connection.execute(
                        "INSERT INTO extraction_found VALUES (?, ?)",
                        (record.name, disposition),
                    )
                except sqlite3.IntegrityError as exc:
                    raise ValueError(
                        f"Duplicate assigned FASTQ ID in {input_path}: {record.name}"
                    ) from exc
            connection.commit()
            if progress:
                progress(records, None, started)
    finally:
        _close_all(writers.values())
    missing, examples = _missing_assignment_values(connection, groups)
    connection.execute("DROP TABLE extraction_found")
    connection.commit()
    connection.close()
    return {
        "input": str(input_path), "input_records": records,
        "requested_assigned_reads": requested, "written_reads": written,
        "filtered_assigned_reads": filtered, "missing_reads": missing,
        "missing_examples": examples,
        "group_written_reads": {
            group: group_written[group] for group in sorted(groups, key=group_sort_key)
        },
        "full_input_scans": 1,
        "filters": {"min_length": min_length, "min_mean_quality": min_mean_quality},
    }


def dispatch_paired_fastq_database(
    input1: Path,
    input2: Path,
    assignment_database: Path,
    groups,
    output_directory: Path,
    progress=None,
):
    writers1 = _group_fastq_writers(output_directory, groups, "Hi-C.1")
    writers2 = _group_fastq_writers(output_directory, groups, "Hi-C.2")
    connection = _sqlite_connect(assignment_database)
    _prepare_extraction_audit(connection)
    requested = _requested_assignment_count(connection, groups)
    pairs = written = 0
    group_written = Counter()
    started = time.monotonic()
    paired_records = itertools.zip_longest(iter_fastq(input1), iter_fastq(input2))
    try:
        for chunk in _fastq_chunks(paired_records, max_records=20_000):
            if any(pair[0] is None or pair[1] is None for pair in chunk):
                raise ValueError("Hi-C FASTQ files contain different record counts")
            assignments = _lookup_assigned_groups(
                connection, [pair[0].name for pair in chunk], groups
            )
            for record1, record2 in chunk:
                pairs += 1
                if record1.name != record2.name:
                    raise ValueError(
                        f"Hi-C mate name mismatch: {record1.name} != {record2.name}"
                    )
                group = assignments.get(record1.name)
                if group is None:
                    continue
                try:
                    connection.execute(
                        "INSERT INTO extraction_found VALUES (?, 'written')",
                        (record1.name,),
                    )
                except sqlite3.IntegrityError as exc:
                    raise ValueError(
                        f"Duplicate assigned Hi-C pair ID: {record1.name}"
                    ) from exc
                writers1[group].write(record1.text)
                writers2[group].write(record2.text)
                written += 1
                group_written[group] += 1
            connection.commit()
            if progress:
                progress(pairs, None, started)
    finally:
        _close_all(writers1.values())
        _close_all(writers2.values())
    missing, examples = _missing_assignment_values(connection, groups)
    connection.execute("DROP TABLE extraction_found")
    connection.commit()
    connection.close()
    return {
        "input1": str(input1), "input2": str(input2), "input_pairs": pairs,
        "requested_assigned_pairs": requested, "written_pairs": written,
        "missing_pairs": missing, "missing_examples": examples,
        "group_written_pairs": {
            group: group_written[group] for group in sorted(groups, key=group_sort_key)
        },
        "full_input_scans_per_mate": 1,
    }


def _group_fastq_writers(output_directory: Path, groups, suffix: str):
    output_directory.mkdir(parents=True, exist_ok=True)
    return {
        group: gzip.open(
            output_directory / f"{group}.{suffix}.fq.gz",
            "wt",
            compresslevel=1,
        )
        for group in sorted(groups, key=group_sort_key)
    }


def _close_all(handles: Iterable[TextIO]):
    first_error = None
    for handle in handles:
        try:
            handle.close()
        except Exception as exc:  # pragma: no cover - defensive close path
            first_error = first_error or exc
    if first_error:
        raise first_error


def dispatch_single_fastq(
    input_path: Path,
    decisions,
    groups,
    output_directory: Path,
    suffix: str,
    min_length: int = 0,
    min_mean_quality: float = 0.0,
):
    writers = _group_fastq_writers(output_directory, groups, suffix)
    requested = {
        read_id for read_id, decision in decisions.items() if decision.status == "assigned"
    }
    found = set()
    records = 0
    written = 0
    filtered = 0
    group_written = Counter()
    try:
        for record in iter_fastq(input_path):
            records += 1
            decision = decisions.get(record.name)
            if decision is None or decision.status != "assigned":
                continue
            if record.name in found:
                raise ValueError(f"Duplicate assigned FASTQ ID in {input_path}: {record.name}")
            found.add(record.name)
            mean_quality = (
                sum(ord(value) - 33 for value in record.quality) / len(record.quality)
                if record.quality else 0.0
            )
            if (
                len(record.sequence) < min_length
                or mean_quality < min_mean_quality
            ):
                filtered += 1
                continue
            writers[decision.group].write(record.text)
            written += 1
            group_written[decision.group] += 1
    finally:
        _close_all(writers.values())
    missing = sorted(requested - found)
    return {
        "input": str(input_path),
        "input_records": records,
        "requested_assigned_reads": len(requested),
        "written_reads": written,
        "filtered_assigned_reads": filtered,
        "missing_reads": len(missing),
        "missing_examples": missing[:10],
        "group_written_reads": {
            group: group_written[group] for group in sorted(groups, key=group_sort_key)
        },
        "full_input_scans": 1,
        "filters": {
            "min_length": min_length,
            "min_mean_quality": min_mean_quality,
        },
    }


def dispatch_paired_fastq(
    input1: Path,
    input2: Path,
    decisions,
    groups,
    output_directory: Path,
):
    writers1 = _group_fastq_writers(output_directory, groups, "Hi-C.1")
    writers2 = _group_fastq_writers(output_directory, groups, "Hi-C.2")
    requested = {
        read_id for read_id, decision in decisions.items() if decision.status == "assigned"
    }
    found = set()
    pairs = 0
    written = 0
    group_written = Counter()
    try:
        for record1, record2 in itertools.zip_longest(
            iter_fastq(input1), iter_fastq(input2)
        ):
            if record1 is None or record2 is None:
                raise ValueError("Hi-C FASTQ files contain different record counts")
            pairs += 1
            if record1.name != record2.name:
                raise ValueError(
                    f"Hi-C mate name mismatch: {record1.name} != {record2.name}"
                )
            decision = decisions.get(record1.name)
            if decision is None or decision.status != "assigned":
                continue
            if record1.name in found:
                raise ValueError(f"Duplicate assigned Hi-C pair ID: {record1.name}")
            writers1[decision.group].write(record1.text)
            writers2[decision.group].write(record2.text)
            found.add(record1.name)
            written += 1
            group_written[decision.group] += 1
    finally:
        _close_all(writers1.values())
        _close_all(writers2.values())
    missing = sorted(requested - found)
    return {
        "input1": str(input1),
        "input2": str(input2),
        "input_pairs": pairs,
        "requested_assigned_pairs": len(requested),
        "written_pairs": written,
        "missing_pairs": len(missing),
        "missing_examples": missing[:10],
        "group_written_pairs": {
            group: group_written[group] for group in sorted(groups, key=group_sort_key)
        },
        "full_input_scans_per_mate": 1,
    }

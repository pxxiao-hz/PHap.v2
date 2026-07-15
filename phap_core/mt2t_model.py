"""Typed data model for conflict-aware mT2T overlap assembly."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from .paf import PafRecord


@dataclass(frozen=True)
class Mt2tOverlapParameters:
    """Thresholds for conservative overlap classification."""

    min_contig_length: int = 100_000
    min_alignment_block_length: int = 200
    min_identity: float = 0.8
    min_shorter_coverage: float = 0.1
    min_longer_coverage: float = 0.05
    internal_margin_ratio: float = 0.05
    max_query_gap: int = 3_000_000
    max_target_gap: int = 3_000_000

    def __post_init__(self) -> None:
        if self.min_contig_length < 1:
            raise ValueError("min_contig_length must be positive")
        if self.min_alignment_block_length < 1:
            raise ValueError("min_alignment_block_length must be positive")
        for name, value in (
            ("min_identity", self.min_identity),
            ("min_shorter_coverage", self.min_shorter_coverage),
            ("min_longer_coverage", self.min_longer_coverage),
            ("internal_margin_ratio", self.internal_margin_ratio),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        if self.max_query_gap < 0 or self.max_target_gap < 0:
            raise ValueError("chain gap thresholds must be non-negative")


@dataclass(frozen=True, order=True)
class OrientedContig:
    contig_id: str
    orientation: str

    def __post_init__(self) -> None:
        if self.orientation not in {"+", "-"}:
            raise ValueError("orientation must be '+' or '-'")

    def flipped(self) -> "OrientedContig":
        return OrientedContig(
            self.contig_id,
            "-" if self.orientation == "+" else "+",
        )


@dataclass(frozen=True)
class PairChainEvidence:
    query_id: str
    target_id: str
    strand: str
    query_length: int
    target_length: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    query_coverage: float
    target_coverage: float
    query_union_bases: int
    target_union_bases: int
    identity: float
    score: float
    collinear: bool
    record_count: int
    records: Tuple[PafRecord, ...]


@dataclass(frozen=True)
class OrientedOverlapEdge:
    left: OrientedContig
    right: OrientedContig
    left_length: int
    right_length: int
    left_interval: Tuple[int, int]
    right_interval: Tuple[int, int]
    query_id: str
    target_id: str
    identity: float
    query_coverage: float
    target_coverage: float
    score: float
    source_records: Tuple[Tuple[str, int], ...]

    @property
    def left_endpoint(self) -> Tuple[str, str]:
        return (
            self.left.contig_id,
            "R" if self.left.orientation == "+" else "L",
        )

    @property
    def right_endpoint(self) -> Tuple[str, str]:
        return (
            self.right.contig_id,
            "L" if self.right.orientation == "+" else "R",
        )

    @property
    def key(self) -> Tuple[object, ...]:
        return (
            self.left.contig_id,
            self.left.orientation,
            self.right.contig_id,
            self.right.orientation,
            self.left_interval,
            self.right_interval,
        )

    def reversed(self) -> "OrientedOverlapEdge":
        return OrientedOverlapEdge(
            left=self.right.flipped(),
            right=self.left.flipped(),
            left_length=self.right_length,
            right_length=self.left_length,
            left_interval=(
                self.right_length - self.right_interval[1],
                self.right_length - self.right_interval[0],
            ),
            right_interval=(
                self.left_length - self.left_interval[1],
                self.left_length - self.left_interval[0],
            ),
            query_id=self.query_id,
            target_id=self.target_id,
            identity=self.identity,
            query_coverage=self.query_coverage,
            target_coverage=self.target_coverage,
            score=self.score,
            source_records=self.source_records,
        )


@dataclass(frozen=True)
class AlignmentAudit:
    source: str
    line_number: int
    query_id: str
    target_id: str
    initial_status: str
    initial_reason: str
    final_status: str
    final_reason: str


@dataclass(frozen=True)
class ChainAudit:
    query_id: str
    target_id: str
    strand: str
    record_count: int
    identity: float
    query_coverage: float
    target_coverage: float
    query_union_bases: int
    target_union_bases: int
    query_span_bases: int
    target_span_bases: int
    score: float
    collinear: bool
    classification: str
    reason: str


@dataclass(frozen=True)
class EdgeAudit:
    edge: OrientedOverlapEdge
    status: str
    reason: str


@dataclass(frozen=True)
class ContainmentAudit:
    child_id: str
    host_id: str
    strand: str
    identity: float
    child_coverage: float
    host_start: int
    host_end: int
    source_records: Tuple[Tuple[str, int], ...]
    score: float
    status: str
    reason: str


@dataclass(frozen=True)
class OrientedPath:
    nodes: Tuple[OrientedContig, ...]
    edges: Tuple[OrientedOverlapEdge, ...]


@dataclass(frozen=True)
class SequenceOutput:
    identifier: str
    header: str
    sequence: str
    source_ids: Tuple[str, ...]


@dataclass(frozen=True)
class RoutingAudit:
    source_id: str
    destination_id: str
    status: str
    reason: str
    evidence: str


@dataclass(frozen=True)
class JoinAudit:
    destination_id: str
    step: int
    left_id: str
    left_orientation: str
    right_id: str
    right_orientation: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    right_trim_bases: int
    right_retained_bases: int
    score: float


@dataclass(frozen=True)
class Mt2tAssemblyResult:
    alignment_audit: Tuple[AlignmentAudit, ...]
    chain_audit: Tuple[ChainAudit, ...]
    edge_audit: Tuple[EdgeAudit, ...]
    containment_audit: Tuple[ContainmentAudit, ...]
    paths: Tuple[OrientedPath, ...]
    outputs: Tuple[SequenceOutput, ...]
    routing: Tuple[RoutingAudit, ...]
    joins: Tuple[JoinAudit, ...]
    contained_ids: Tuple[str, ...]


@dataclass(frozen=True)
class _ClassifiedChain:
    chain: PairChainEvidence
    classification: str
    reason: str
    edge: Optional[OrientedOverlapEdge] = None
    child_id: Optional[str] = None
    host_id: Optional[str] = None


def reverse_complement(sequence: str) -> str:
    """Return an IUPAC-aware reverse complement, preserving letter case."""

    alphabet = "ACGTRYMKBDHVNacgtrymkbdhvn"
    complements = "TGCAYRKMVHDBNtgcayrkmv hdbn".replace(" ", "")
    table = str.maketrans(alphabet, complements)
    translated = sequence.translate(table)
    invalid = sorted({base for base in translated if base not in alphabet})
    if invalid:
        raise ValueError(
            "sequence contains unsupported symbols: " + ", ".join(invalid)
        )
    return translated[::-1]

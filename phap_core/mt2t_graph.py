"""Resolve oriented overlap graphs and materialize unambiguous mT2T paths."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Dict, Tuple

from .fasta import FastaRecord
from .mt2t_model import (
    _ClassifiedChain,
    ContainmentAudit,
    EdgeAudit,
    JoinAudit,
    Mt2tOverlapParameters,
    OrientedOverlapEdge,
    OrientedPath,
    RoutingAudit,
    SequenceOutput,
    reverse_complement,
)


def _resolve_graph(
    edges: Sequence[OrientedOverlapEdge],
) -> Tuple[
    Tuple[OrientedOverlapEdge, ...],
    Tuple[EdgeAudit, ...],
    Dict[str, str],
]:
    components = _edge_components(edges)
    accepted: list[OrientedOverlapEdge] = []
    audits: list[EdgeAudit] = []
    node_conflicts: Dict[str, str] = {}
    for component_nodes, component_edges in components:
        endpoint_counts: Dict[Tuple[str, str], int] = defaultdict(int)
        pair_geometries: Dict[Tuple[str, str], set[Tuple[object, ...]]] = defaultdict(set)
        for edge in component_edges:
            endpoint_counts[edge.left_endpoint] += 1
            endpoint_counts[edge.right_endpoint] += 1
            pair_ids = sorted((edge.left.contig_id, edge.right.contig_id))
            pair = (pair_ids[0], pair_ids[1])
            pair_geometries[pair].add(
                tuple(sorted((edge.left_endpoint, edge.right_endpoint)))
            )
        if any(len(geometries) > 1 for geometries in pair_geometries.values()):
            reason = "orientation_conflict"
        elif any(count > 1 for count in endpoint_counts.values()):
            reason = "branch_conflict"
        elif len(component_edges) >= len(component_nodes):
            reason = "cycle_conflict"
        else:
            reason = "accepted_unambiguous_path"
        if reason == "accepted_unambiguous_path":
            try:
                _walk_component(component_nodes, component_edges)
            except ValueError:
                reason = "orientation_conflict"
        if reason == "accepted_unambiguous_path":
            accepted.extend(component_edges)
            status = "accepted"
        else:
            status = "rejected"
            for node in component_nodes:
                node_conflicts[node] = reason
        audits.extend(EdgeAudit(edge, status, reason) for edge in component_edges)
    return (
        tuple(sorted(accepted, key=lambda edge: edge.key)),
        tuple(sorted(audits, key=lambda row: row.edge.key)),
        node_conflicts,
    )


def _extract_paths(edges: Sequence[OrientedOverlapEdge]) -> Tuple[OrientedPath, ...]:
    paths = [
        _walk_component(nodes, component_edges)
        for nodes, component_edges in _edge_components(edges)
    ]
    return tuple(sorted(paths, key=lambda path: tuple(path.nodes)))


def _walk_component(
    nodes: Sequence[str],
    edges: Sequence[OrientedOverlapEdge],
) -> OrientedPath:
    adjacency: Dict[str, list[OrientedOverlapEdge]] = defaultdict(list)
    for edge in edges:
        adjacency[edge.left.contig_id].append(edge)
        adjacency[edge.right.contig_id].append(edge)
    endpoints = sorted(node for node in nodes if len(adjacency[node]) == 1)
    if len(endpoints) != 2:
        raise ValueError("an unambiguous path must have exactly two endpoints")
    current_id = endpoints[0]
    first = adjacency[current_id][0]
    oriented = first if first.left.contig_id == current_id else first.reversed()
    path_nodes = [oriented.left, oriented.right]
    path_edges = [oriented]
    used = {first.key}
    while True:
        current = path_nodes[-1]
        remaining = [edge for edge in adjacency[current.contig_id] if edge.key not in used]
        if not remaining:
            break
        if len(remaining) != 1:
            raise ValueError("path traversal encountered a branch")
        source_edge = remaining[0]
        if (
            source_edge.left.contig_id == current.contig_id
            and source_edge.left.orientation == current.orientation
        ):
            next_edge = source_edge
        elif (
            source_edge.right.contig_id == current.contig_id
            and source_edge.right.flipped().orientation == current.orientation
        ):
            next_edge = source_edge.reversed()
        else:
            raise ValueError("path traversal encountered an orientation conflict")
        used.add(source_edge.key)
        path_edges.append(next_edge)
        path_nodes.append(next_edge.right)
    if len(used) != len(edges):
        raise ValueError("path traversal did not consume every edge")
    return OrientedPath(tuple(path_nodes), tuple(path_edges))


def _materialize_outputs(
    fasta_records: Sequence[FastaRecord],
    paths: Sequence[OrientedPath],
    contained: Mapping[str, str],
    node_conflicts: Mapping[str, str],
    parameters: Mt2tOverlapParameters,
) -> Tuple[Tuple[SequenceOutput, ...], Tuple[RoutingAudit, ...], Tuple[JoinAudit, ...]]:
    records_by_id = {record.identifier: record for record in fasta_records}
    source_ids = set(records_by_id)
    path_destinations: Dict[str, Tuple[str, str]] = {}
    outputs = []
    joins = []
    for index, path in enumerate(paths, start=1):
        destination = f"PHap_mT2T_path_{index:06d}"
        if destination in source_ids:
            raise ValueError(f"generated mT2T identifier collides with source {destination!r}")
        first = path.nodes[0]
        sequence = _oriented_sequence(records_by_id[first.contig_id].sequence, first.orientation)
        for step, edge in enumerate(path.edges, start=1):
            right_sequence = _oriented_sequence(
                records_by_id[edge.right.contig_id].sequence,
                edge.right.orientation,
            )
            trim = edge.right_interval[1]
            if not 0 <= trim <= len(right_sequence):
                raise AssertionError("validated overlap produced an invalid trim")
            sequence += right_sequence[trim:]
            joins.append(
                JoinAudit(
                    destination_id=destination,
                    step=step,
                    left_id=edge.left.contig_id,
                    left_orientation=edge.left.orientation,
                    right_id=edge.right.contig_id,
                    right_orientation=edge.right.orientation,
                    left_start=edge.left_interval[0],
                    left_end=edge.left_interval[1],
                    right_start=edge.right_interval[0],
                    right_end=edge.right_interval[1],
                    right_trim_bases=trim,
                    right_retained_bases=len(right_sequence) - trim,
                    score=edge.score,
                )
            )
        path_ids = tuple(node.contig_id for node in path.nodes)
        outputs.append(SequenceOutput(destination, destination, sequence, path_ids))
        for node in path.nodes:
            path_destinations[node.contig_id] = (destination, node.orientation)

    def containment_root(contig_id: str) -> str:
        seen = set()
        current = contig_id
        while current in contained:
            if current in seen:
                raise ValueError("containment cycle detected")
            seen.add(current)
            current = contained[current]
        return current

    for record in sorted(fasta_records, key=lambda row: row.identifier):
        if record.identifier in contained or record.identifier in path_destinations:
            continue
        outputs.append(
            SequenceOutput(
                record.identifier,
                record.header,
                record.sequence,
                (record.identifier,),
            )
        )
    destination_by_source = {
        source_id: destination
        for source_id, (destination, _) in path_destinations.items()
    }
    for child_id in contained:
        root = containment_root(child_id)
        destination_by_source[child_id] = destination_by_source.get(root, root)
    routing = []
    for record in sorted(fasta_records, key=lambda row: row.identifier):
        source_id = record.identifier
        if source_id in contained:
            root = containment_root(source_id)
            destination = destination_by_source[source_id]
            status = "contained"
            reason = "unique_internal_containment"
            evidence = f"host={root}"
        elif source_id in path_destinations:
            destination, orientation = path_destinations[source_id]
            status = "merged"
            reason = "accepted_unambiguous_overlap_path"
            evidence = f"orientation={orientation}"
        else:
            destination = source_id
            status = "preserved"
            if len(record.sequence) < parameters.min_contig_length:
                reason = "below_min_contig_length_preserved"
            else:
                reason = node_conflicts.get(source_id, "no_supported_terminal_overlap")
            evidence = "original_sequence"
        routing.append(RoutingAudit(source_id, destination, status, reason, evidence))
    return (
        tuple(sorted(outputs, key=lambda row: row.identifier)),
        tuple(routing),
        tuple(sorted(joins, key=lambda row: (row.destination_id, row.step))),
    )


def _final_alignment_reasons(
    classified: Sequence[_ClassifiedChain],
    edge_audits: Sequence[EdgeAudit],
    containment_audits: Sequence[ContainmentAudit],
) -> Dict[Tuple[str, int], str]:
    edge_reason = {
        source_record: audit.reason
        for audit in edge_audits
        for source_record in audit.edge.source_records
    }
    containment_reason = {
        (row.child_id, row.host_id): row.reason for row in containment_audits
    }
    final = {}
    for item in classified:
        if item.edge is not None:
            reason = edge_reason.get(item.edge.source_records[0], item.reason)
        elif item.child_id is not None and item.host_id is not None:
            reason = containment_reason.get((item.child_id, item.host_id), item.reason)
        else:
            reason = item.reason
        for record in item.chain.records:
            final[(record.source, record.line_number)] = reason
    return final


def _edge_components(
    edges: Sequence[OrientedOverlapEdge],
) -> Tuple[Tuple[Tuple[str, ...], Tuple[OrientedOverlapEdge, ...]], ...]:
    adjacency: Dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        left = edge.left.contig_id
        right = edge.right.contig_id
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(adjacency)
    components = []
    while unseen:
        start = min(unseen)
        stack = [start]
        nodes = set()
        while stack:
            node = stack.pop()
            if node in nodes:
                continue
            nodes.add(node)
            stack.extend(sorted(adjacency[node] - nodes, reverse=True))
        unseen -= nodes
        component_edges = tuple(
            sorted(
                (
                    edge
                    for edge in edges
                    if edge.left.contig_id in nodes and edge.right.contig_id in nodes
                ),
                key=lambda edge: edge.key,
            )
        )
        components.append((tuple(sorted(nodes)), component_edges))
    return tuple(components)



def _oriented_sequence(sequence: str, orientation: str) -> str:
    return sequence if orientation == "+" else reverse_complement(sequence)

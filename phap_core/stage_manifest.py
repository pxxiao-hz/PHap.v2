"""Validated, content-addressed manifests for resumable PHap stages."""

from __future__ import annotations

import hashlib
import json
import os
import fcntl
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Tuple, Union

from . import __version__
from .atomic_io import atomic_write_lines


PathLike = Union[str, os.PathLike[str]]
MANIFEST_SCHEMA = "phap-stage-manifest-v1"


@dataclass(frozen=True)
class StageInput:
    """One named input whose path and content are part of a cache key."""

    name: str
    path: PathLike


@dataclass(frozen=True)
class StageSpec:
    """Complete non-output state that determines one stage result."""

    stage: str
    inputs: Tuple[StageInput, ...]
    parameters: Tuple[Tuple[str, str], ...]
    tool_versions: Tuple[Tuple[str, str], ...]
    phap_version: str = __version__

    def __post_init__(self) -> None:
        if not self.stage:
            raise ValueError("stage must not be empty")
        if not self.phap_version:
            raise ValueError("phap_version must not be empty")
        _require_unique_names("input", tuple(row.name for row in self.inputs))
        _require_unique_names("parameter", tuple(name for name, _ in self.parameters))
        _require_unique_names("tool", tuple(name for name, _ in self.tool_versions))
        if any(not name or not value for name, value in self.parameters):
            raise ValueError("parameter names and values must not be empty")
        if any(not name or not version for name, version in self.tool_versions):
            raise ValueError("tool names and versions must not be empty")


@dataclass(frozen=True)
class CacheValidation:
    """Result of validating a manifest and all declared outputs."""

    hit: bool
    reason: str


@dataclass(frozen=True)
class StageSnapshot:
    """A stage signature frozen before external work starts."""

    spec: StageSpec
    signature_json: str


def capture_stage(spec: StageSpec) -> StageSnapshot:
    """Fingerprint all inputs once, before any output is generated."""

    return StageSnapshot(spec, _canonical_json(_signature(spec)))


def validate_stage_cache(
    snapshot: StageSnapshot,
    manifest_path: PathLike,
    outputs: Mapping[str, PathLike],
) -> CacheValidation:
    """Return a cache hit only when manifest, inputs, and outputs all match."""

    _validate_output_mapping(outputs)
    manifest = Path(manifest_path)
    if not manifest.is_file():
        return CacheValidation(False, "manifest_missing")
    try:
        observed = json.loads(
            manifest.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        return CacheValidation(False, "manifest_invalid")
    if not isinstance(observed, dict):
        return CacheValidation(False, "manifest_invalid")

    expected_signature = json.loads(snapshot.signature_json)
    if observed.get("schema") != MANIFEST_SCHEMA:
        return CacheValidation(False, "manifest_schema_mismatch")
    if observed.get("signature") != expected_signature:
        return CacheValidation(False, "stage_signature_mismatch")

    declared_outputs = observed.get("outputs")
    if not isinstance(declared_outputs, list):
        return CacheValidation(False, "manifest_outputs_invalid")
    expected_names = tuple(sorted(outputs))
    observed_names = tuple(
        row.get("name") for row in declared_outputs if isinstance(row, dict)
    )
    if observed_names != expected_names or len(declared_outputs) != len(outputs):
        return CacheValidation(False, "output_set_mismatch")
    by_name = {
        row["name"]: row for row in declared_outputs if isinstance(row, dict)
    }
    for name in expected_names:
        path = Path(outputs[name]).resolve()
        row = by_name[name]
        if row.get("path") != str(path):
            return CacheValidation(False, "output_path_mismatch")
        if not path.is_file():
            return CacheValidation(False, "output_missing")
        if row.get("size") != path.stat().st_size:
            return CacheValidation(False, "output_size_mismatch")
        if row.get("sha256") != sha256_file(path):
            return CacheValidation(False, "output_hash_mismatch")
    return CacheValidation(True, "cache_valid")


def write_stage_manifest(
    snapshot: StageSnapshot,
    manifest_path: PathLike,
    outputs: Mapping[str, PathLike],
) -> None:
    """Fingerprint successful outputs and atomically write the manifest last."""

    _validate_output_mapping(outputs)
    output_rows = []
    for name in sorted(outputs):
        path = Path(outputs[name]).resolve()
        if not path.is_file():
            raise ValueError(f"stage output is missing: {path}")
        output_rows.append(
            {
                "name": name,
                "path": str(path),
                "sha256": sha256_file(path),
                "size": path.stat().st_size,
            }
        )
    current_signature = _canonical_json(_signature(snapshot.spec))
    if current_signature != snapshot.signature_json:
        raise ValueError("stage inputs changed while the stage was running")
    document = {
        "outputs": output_rows,
        "schema": MANIFEST_SCHEMA,
        "signature": json.loads(snapshot.signature_json),
    }
    atomic_write_lines(
        manifest_path,
        (_canonical_json(document),),
    )


def invalidate_stage_manifest(manifest_path: PathLike) -> None:
    """Remove a completion marker before rebuilding a stage."""

    Path(manifest_path).unlink(missing_ok=True)


@contextmanager
def stage_lock(manifest_path: PathLike) -> Iterator[None]:
    """Serialize builders for one stage; callers recheck the cache inside."""

    manifest = Path(manifest_path)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    lock_path = manifest.with_name(f".{manifest.name}.lock")
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def sha256_file(path: PathLike) -> str:
    """Return the SHA-256 digest of one regular file."""

    source = Path(path)
    if not source.is_file():
        raise ValueError(f"cannot fingerprint non-file input: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(spec: StageSpec) -> dict[str, object]:
    inputs = []
    for row in sorted(spec.inputs, key=lambda item: item.name):
        path = Path(row.path).resolve()
        if not path.is_file():
            raise ValueError(f"stage input is missing: {path}")
        inputs.append(
            {
                "name": row.name,
                "path": str(path),
                "sha256": sha256_file(path),
                "size": path.stat().st_size,
            }
        )
    return {
        "inputs": inputs,
        "parameters": [
            {"name": name, "value": value}
            for name, value in sorted(spec.parameters)
        ],
        "phap_version": spec.phap_version,
        "stage": spec.stage,
        "tools": [
            {"name": name, "version": version}
            for name, version in sorted(spec.tool_versions)
        ],
    }


def _validate_output_mapping(outputs: Mapping[str, PathLike]) -> None:
    if not outputs:
        raise ValueError("at least one stage output is required")
    _require_unique_names("output", tuple(outputs))
    if any(not name for name in outputs):
        raise ValueError("output names must not be empty")


def _require_unique_names(kind: str, names: Tuple[str, ...]) -> None:
    if len(names) != len(set(names)):
        raise ValueError(f"duplicate {kind} name")
    if any(not name for name in names):
        raise ValueError(f"{kind} names must not be empty")


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result

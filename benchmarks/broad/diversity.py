"""Diversity metrics and strain stratification for the broad experiment.

The metric is a bottom-k MinHash sketch of all canonical 21-mers. Jaccard
similarity is converted to a Mash distance and an approximate ANI. Genomes are
clustered with average-linkage hierarchical clustering, then three disjoint
fixed-size strata are selected:

* low    -- the most central members of the tightest cluster;
* high   -- a farthest-point sample of the remainder;
* medium -- the remaining records.

The selection is deterministic and the measured pairwise distances are stored
alongside the strain lists, so a stratum is never defined by its label alone.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform

K = 21
STRIDE = 1
SKETCH_SIZE = 20_000
HASH_SEED = 12345

_TABLE = np.zeros(256, dtype=np.uint8)
for _base, _value in ((b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3)):
    _TABLE[_base[0]] = _value
    _TABLE[_base[0] + 32] = _value

_SPLIT_1 = np.uint64(0x9E3779B97F4A7C15)
_SPLIT_2 = np.uint64(0xBF58476D1CE4E5B9)
_SPLIT_3 = np.uint64(0x94D049BB133111EB)
_SHIFT_30 = np.uint64(30)
_SHIFT_27 = np.uint64(27)
_SHIFT_31 = np.uint64(31)


def _splitmix64(value: np.ndarray) -> np.ndarray:
    value = value + _SPLIT_1
    value = (value ^ (value >> _SHIFT_30)) * _SPLIT_2
    value = (value ^ (value >> _SHIFT_27)) * _SPLIT_3
    return value ^ (value >> _SHIFT_31)


def sketch_sequence(sequence: bytes, k: int = K, stride: int = STRIDE,
                    size: int = SKETCH_SIZE, chunk: int = 100_000) -> np.ndarray:
    """Bottom-k hashes over all canonical k-mers (chunked to bound memory)."""
    if len(sequence) < k:
        return np.empty(0, dtype=np.uint64)
    encoded = _TABLE[np.frombuffer(sequence, dtype=np.uint8)]
    shifts = np.arange(k, dtype=np.uint64) * np.uint64(2)
    total = len(encoded) - k + 1
    best = np.empty(0, dtype=np.uint64)
    for start in range(0, total, chunk):
        stop = min(start + chunk, total)
        windows = np.lib.stride_tricks.sliding_window_view(
            encoded[start : stop + k - 1], k
        )
        if stride > 1:
            windows = windows[::stride]
        forward = (windows.astype(np.uint64) << shifts).sum(axis=1)
        reverse = (3 - windows[:, ::-1]).astype(np.uint64)
        backward = (reverse << shifts).sum(axis=1)
        hashes = np.minimum(_splitmix64(forward), _splitmix64(backward))
        merged = np.unique(np.concatenate([best, hashes]))
        if len(merged) > size:
            merged = np.partition(merged, size - 1)[:size]
        best = merged
    return np.sort(best)


def load_sketch_cache(path: Path, manifest: dict) -> list[dict] | None:
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as data:
        if str(data["source_sha256"][0]) != manifest["fasta"]["sha256"]:
            return None
        if list(data["accessions"]) != [row["accession"] for row in manifest["records"]]:
            return None
        if len(set(int(value) for value in data["sketch_lengths"])) != 1:
            return None
        matrix = data["sketches"]
        return [
            {
                "accession": str(data["accessions"][index]),
                "sequence_length": int(data["sequence_length"][index]),
                "gc_fraction": float(data["gc_fraction"][index]),
                "non_acgt_bases": int(data["non_acgt_bases"][index]),
                "sketch": matrix[index].copy(),
            }
            for index in range(len(matrix))
        ]


def save_sketch_cache(path: Path, sketches: list[dict], source_sha256: str) -> bool:
    lengths = [len(entry["sketch"]) for entry in sketches]
    if len(set(lengths)) != 1:
        return False
    np.savez(
        path,
        sketches=np.stack([entry["sketch"] for entry in sketches]),
        accessions=np.array([entry["accession"] for entry in sketches]),
        sequence_length=np.array([entry["sequence_length"] for entry in sketches], dtype=np.int64),
        gc_fraction=np.array([entry["gc_fraction"] for entry in sketches], dtype=np.float64),
        non_acgt_bases=np.array([entry["non_acgt_bases"] for entry in sketches], dtype=np.int64),
        sketch_lengths=np.array(lengths, dtype=np.int64),
        source_sha256=np.array([source_sha256]),
    )
    return True


def sketch_manifest(pool_path: Path, manifest: dict, cache_path: Path | None = None,
                    refresh: bool = False) -> list[dict]:
    if cache_path is not None and cache_path.exists() and not refresh:
        cached = load_sketch_cache(cache_path, manifest)
        if cached is not None:
            print(f"loaded sketches from {cache_path.name}", flush=True)
            return cached
    sketches = []
    with pool_path.open("rb") as stream:
        for index, record in enumerate(manifest["records"], start=1):
            stream.seek(record["file_offset"])
            payload = stream.read(record["file_length"])
            sequence = b"".join(
                line.strip() for line in payload.splitlines()
                if not line.startswith(b">")
            ).upper()
            sketches.append(
                {
                    "accession": record["accession"],
                    "sequence_length": record["sequence_length"],
                    "gc_fraction": record["gc_fraction"],
                    "non_acgt_bases": record["non_acgt_bases"],
                    "sketch": sketch_sequence(sequence),
                }
            )
            if index % 12 == 0:
                print(f"  sketched {index}/{len(manifest['records'])}", flush=True)
    if cache_path is not None:
        if save_sketch_cache(cache_path, sketches, manifest["fasta"]["sha256"]):
            print(f"wrote {cache_path.name}", flush=True)
    return sketches


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    if not len(a) or not len(b):
        return 0.0
    intersection = np.intersect1d(a, b, assume_unique=True).size
    return intersection / (len(a) + len(b) - intersection)


def mash_distance(jaccard_value: float, k: int = K) -> float:
    if jaccard_value <= 0.0:
        return 1.0
    return min(1.0, max(0.0, -math.log(2.0 * jaccard_value / (1.0 + jaccard_value)) / k))


def distance_matrix(sketches: list[dict]) -> np.ndarray:
    count = len(sketches)
    matrix = np.zeros((count, count), dtype=np.float64)
    for i in range(count):
        for j in range(i + 1, count):
            value = mash_distance(jaccard(sketches[i]["sketch"], sketches[j]["sketch"]))
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix


def _merge_to_size(labels: np.ndarray, distances: np.ndarray, minimum: int) -> np.ndarray:
    labels = labels.copy()
    while True:
        unique = [value for value in np.unique(labels)]
        sizes = {value: int((labels == value).sum()) for value in unique}
        if max(sizes.values()) >= minimum or len(unique) == 1:
            break
        smallest = min(unique, key=lambda value: sizes[value])
        others = [value for value in unique if value != smallest]
        nearest = min(
            others,
            key=lambda value: distances[
                np.ix_(np.where(labels == smallest)[0], np.where(labels == value)[0])
            ].mean(),
        )
        labels[labels == smallest] = nearest
    remap = {value: index for index, value in enumerate(np.unique(labels))}
    return np.array([remap[value] for value in labels], dtype=np.int64)


def cluster_and_label(distances: np.ndarray, target_n: int) -> tuple[np.ndarray, int]:
    count = len(distances)
    k = min(4, max(2, count // 18))
    tree = linkage(squareform(distances, checks=False), method="average")
    labels = cut_tree(tree, n_clusters=k).ravel()
    labels = _merge_to_size(labels, distances, minimum=6)
    return labels, len(np.unique(labels))


def _central(indices: np.ndarray, n: int, distances: np.ndarray) -> list[int]:
    if len(indices) <= n:
        return [int(value) for value in indices]
    block = distances[np.ix_(indices, indices)]
    order = np.argsort(block.sum(axis=1))
    return [int(indices[value]) for value in order[:n]]


def _farthest(indices: np.ndarray, n: int, distances: np.ndarray) -> list[int]:
    if len(indices) <= n:
        return [int(value) for value in indices]
    block = distances[np.ix_(indices, indices)]
    medoid = int(indices[int(np.argmin(block.sum(axis=1)))])
    selected = [medoid]
    remaining = set(int(value) for value in indices) - {medoid}
    while len(selected) < n and remaining:
        best = max(remaining, key=lambda value: min(distances[value, chosen] for chosen in selected))
        selected.append(best)
        remaining.remove(best)
    return selected


def _unique(values: list[int], limit: int) -> list[int]:
    seen: set[int] = set()
    result: list[int] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            result.append(value)
        if len(result) == limit:
            break
    return result


def _mean_pairwise(distances: np.ndarray, indices: list[int]) -> float:
    if len(indices) < 2:
        return 0.0
    block = distances[np.ix_(indices, indices)]
    upper = block[np.triu_indices(len(indices), k=1)]
    return float(upper.mean())


def _stratum(distances: np.ndarray, sketches: list[dict], indices: list[int],
             labels: np.ndarray, name: str) -> dict:
    unique_indices = sorted(set(indices))
    block = distances[np.ix_(unique_indices, unique_indices)]
    upper = block[np.triu_indices(len(unique_indices), k=1)] if len(unique_indices) > 1 else np.array([0.0])
    covered = sorted(set(int(labels[value]) for value in unique_indices))
    bases = sum(sketches[value]["sequence_length"] for value in unique_indices)
    gc = sum(
        sketches[value]["gc_fraction"] * sketches[value]["sequence_length"]
        for value in unique_indices
    ) / bases if bases else 0.0
    return {
        "name": name,
        "indices": unique_indices,
        "records": [sketches[value]["accession"] for value in unique_indices],
        "record_count": len(unique_indices),
        "total_bases": bases,
        "gc_fraction": gc,
        "non_acgt_bases": sum(sketches[value]["non_acgt_bases"] for value in unique_indices),
        "clusters_represented": covered,
        "mean_pairwise_distance": float(upper.mean()),
        "median_pairwise_distance": float(np.median(upper)),
        "min_pairwise_distance": float(upper.min()),
        "max_pairwise_distance": float(upper.max()),
        "mean_mash_ani": float(1.0 - upper.mean()),
    }


def qc_filter(sketches: list[dict], distances: np.ndarray,
              minimum: int = 48) -> tuple[np.ndarray, dict]:
    """Drop length/composition outliers and non-con-specific records."""
    count = len(sketches)
    offset = distances + np.eye(count) * 9.0
    median_distance = np.median(offset, axis=1)
    lengths = np.array([entry["sequence_length"] for entry in sketches], dtype=np.float64)
    median_length = float(np.median(lengths))
    non_acgt = np.array([entry["non_acgt_bases"] for entry in sketches], dtype=np.float64)

    def mask(distance_floor: float, distance_factor: float,
             low_length: float, high_length: float, max_non_acgt: float) -> np.ndarray:
        threshold = max(distance_floor, distance_factor * float(np.median(median_distance)))
        return (
            (median_distance <= threshold)
            & (lengths >= low_length * median_length)
            & (lengths <= high_length * median_length)
            & (non_acgt <= max_non_acgt * lengths)
        )

    keep = np.where(mask(0.05, 3.0, 0.6, 1.4, 0.05))[0]
    relaxed = False
    if len(keep) < minimum:
        keep = np.where(mask(0.08, 4.0, 0.5, 1.5, 0.10))[0]
        relaxed = True
    if len(keep) < 6:
        keep = np.arange(count)
        relaxed = True
    dropped = sorted(set(range(count)) - set(int(value) for value in keep))
    return keep, {
        "kept": len(keep),
        "dropped": len(dropped),
        "dropped_accessions": [sketches[index]["accession"] for index in dropped],
        "relaxed": relaxed,
        "median_length": median_length,
        "median_distance_min": float(median_distance.min()),
        "median_distance_max": float(median_distance.max()),
        "threshold": max(0.05, 3.0 * float(np.median(median_distance))),
    }


def select_strata(distances: np.ndarray, sketches: list[dict],
                  target_n: int = 24) -> dict:
    """Partition the pool into disjoint low/medium/high diversity strata."""
    count = len(distances)
    n = min(target_n, count // 3)
    if n < 2:
        raise ValueError("pool is too small for three strata")
    labels, cluster_count = cluster_and_label(distances, n)
    clusters = [np.where(labels == value)[0] for value in range(cluster_count)]
    within = {
        index: float(distances[np.ix_(cluster, cluster)][np.triu_indices(len(cluster), k=1)].mean())
        if len(cluster) > 1 else 0.0
        for index, cluster in enumerate(clusters)
    }
    base = min(
        (index for index, cluster in enumerate(clusters) if len(cluster) >= n),
        key=lambda index: within[index],
        default=max(range(cluster_count), key=lambda index: len(clusters[index])),
    )

    low = _central(clusters[base], n, distances)
    remaining = [index for index in range(count) if index not in set(low)]
    high = _farthest(np.array(remaining, dtype=np.int64), n, distances)
    medium = [index for index in remaining if index not in set(high)]
    if len(medium) > n:
        # Larger pools: keep the n records closest to the midpoint between the
        # low and high distributions rather than an arbitrary contiguity.
        midpoint = 0.5 * (
            _mean_pairwise(distances, low) + _mean_pairwise(distances, high)
        )
        ordered = sorted(
            medium,
            key=lambda index: abs(
                float(distances[index, low].mean()) - midpoint
            ),
        )
        medium = ordered[:n]
    medium = medium[:n]

    strata = {
        name: _stratum(distances, sketches, indices, labels, name)
        for name, indices in (("low", low), ("medium", medium), ("high", high))
    }
    means = [strata[name]["mean_pairwise_distance"] for name in ("low", "medium", "high")]
    return {
        "method": (
            "canonical 21-mer bottom-20000 MinHash over all k-mers; "
            "Jaccard -> Mash distance; average-linkage clustering"
        ),
        "k": int(cluster_count),
        "strategy": (
            "low = central members of the tightest cluster; high = farthest-point "
            "sample of the remainder; medium = remaining records"
        ),
        "stratum_size": n,
        "monotonic": bool(means[0] <= means[1] <= means[2]),
        "means": means,
        "span": means[2] - means[0],
        "clusters": [
            {
                "id": index,
                "size": len(cluster),
                "within_mean_distance": within[index],
                "mean_distance_from_base": (
                    0.0 if index == base else float(
                        distances[np.ix_(clusters[base], cluster)].mean()
                    )
                ),
            }
            for index, cluster in enumerate(clusters)
        ],
        "strata": strata,
    }


def natural_selection(distances: np.ndarray, sketches: list[dict],
                      target_bases: int = 120_000_000) -> dict:
    order = _farthest(np.arange(len(distances)), len(distances), distances)
    selected: list[int] = []
    total = 0
    for index in order:
        size = sketches[index]["sequence_length"]
        if selected and total >= target_bases:
            break
        selected.append(index)
        total += size
        if total >= target_bases and total <= target_bases * 1.05:
            break
    block = distances[np.ix_(selected, selected)]
    upper = block[np.triu_indices(len(selected), k=1)] if len(selected) > 1 else np.array([0.0])
    return {
        "records": [sketches[index]["accession"] for index in selected],
        "indices": selected,
        "record_count": len(selected),
        "total_bases": total,
        "mean_pairwise_distance": float(upper.mean()),
        "median_pairwise_distance": float(np.median(upper)),
        "mean_mash_ani": float(1.0 - upper.mean()),
        "selection": (
            "farthest-point order from the medoid; stop at the first record "
            f"at or above {target_bases:,} bases"
        ),
    }


def write_subset(pool_path: Path, records: list[dict], indices: list[int],
                 destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    bases = 0
    count = 0
    with pool_path.open("rb") as source, destination.open("wb") as out:
        for index in indices:
            record = records[index]
            source.seek(record["file_offset"])
            payload = source.read(record["file_length"])
            out.write(payload)
            digest.update(payload)
            bases += record["sequence_length"]
            count += 1
    return {
        "path": str(destination),
        "size_bytes": destination.stat().st_size,
        "sha256": digest.hexdigest(),
        "records": count,
        "bases": bases,
    }


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))

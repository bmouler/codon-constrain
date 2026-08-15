"""Deterministic end-to-end benchmark for exact constrained optimization."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

_PROTEIN = "MKTIIALSYIFCLVFAQKTLGAVLGKDSTNVGDEGGFAPNILENKEALELLK" * 2
_GC_RANGE = (0.42, 0.58)
_MOTIFS = ("GAATTC", "GGATCC", "GGTCTC", "GAGACC")
_HOMOPOLYMER_MAX = 3
_FLANK_5 = "ACGTTG"
_FLANK_3 = "CGTACA"
_EXPECTED_SHA256 = "8e0e474a503df6747754560baa0cb9fb92b18ab09ece8b89f226cec781a2d7ff"


def _run(optimize: Any, constraints: Any) -> tuple[Any, str]:
    result = optimize(
        _PROTEIN,
        "ecoli",
        input_type="protein",
        constraints=constraints,
    )
    payload = json.dumps(
        asdict(result), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return result, hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=15)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--print-checksum", action="store_true")
    parser.add_argument(
        "--source-root",
        type=Path,
        help="repository root whose src/ should be imported (default: installed package)",
    )
    args = parser.parse_args()
    if args.samples < 11:
        parser.error("--samples must be at least 11")
    if args.warmups < 1:
        parser.error("--warmups must be at least 1")
    if args.source_root is not None:
        sys.path.insert(0, str(args.source_root.resolve() / "src"))

    from codon_constrain import Constraints, optimize

    constraints = Constraints(
        gc_min=_GC_RANGE[0],
        gc_max=_GC_RANGE[1],
        forbidden_motifs=_MOTIFS,
        homopolymer_max=_HOMOPOLYMER_MAX,
        flank_5=_FLANK_5,
        flank_3=_FLANK_3,
    )
    first, checksum = _run(optimize, constraints)
    if args.print_checksum:
        print(checksum)
        return
    if checksum != _EXPECTED_SHA256:
        raise RuntimeError(f"expected result digest {_EXPECTED_SHA256}, got {checksum}")
    if (
        first.exact is not True
        or first.report["exact_optimum"] is not True
        or first.report["translation_preserved"] is not True
    ):
        raise RuntimeError("optimizer no longer reports an exact translation-preserving result")
    expected_constraints = {
        "gc_in_range": True,
        "motif_free_both_strands": True,
        "homopolymer_ok": True,
        "forbidden_motifs_checked": [
            "GAATTC",
            "GAGACC",
            "GGATCC",
            "GGTCTC",
        ],
    }
    if first.report["constraints"] != expected_constraints:
        raise RuntimeError(f"optimizer constraint report changed: {first.report['constraints']!r}")

    for _ in range(args.warmups - 1):
        _, warmup_checksum = _run(optimize, constraints)
        if warmup_checksum != _EXPECTED_SHA256:
            raise RuntimeError("optimizer result changed during warmup")

    timings: list[float] = []
    for _ in range(args.samples):
        start = time.perf_counter_ns()
        _, sample_checksum = _run(optimize, constraints)
        elapsed = (time.perf_counter_ns() - start) / 1_000_000_000
        if sample_checksum != _EXPECTED_SHA256:
            raise RuntimeError("optimizer result changed between timed samples")
        timings.append(elapsed)

    output = {
        "benchmark": "public optimize() exact constrained optimization",
        "checksum_sha256": checksum,
        "dimensions": {
            "protein_residues": len(_PROTEIN),
            "coding_bases": len(first.coding_sequence),
            "construct_bases": len(first.sequence),
            "forbidden_motifs_input": len(_MOTIFS),
            "forbidden_motifs_expanded": len(
                first.report["constraints"]["forbidden_motifs_checked"]
            ),
            "gc_range": list(_GC_RANGE),
            "homopolymer_max": _HOMOPOLYMER_MAX,
            "flank_bases": len(_FLANK_5) + len(_FLANK_3),
        },
        "samples": len(timings),
        "warmups": args.warmups,
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "max_seconds": max(timings),
    }
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    main()

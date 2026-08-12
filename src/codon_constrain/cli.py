"""Command-line interface for codon-constrain."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .optimizer import Constraints, OptimizationError, optimize


def _read_fasta(path: Path) -> tuple[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise OptimizationError(f"cannot read input FASTA: {error}") from error
    if not lines or not lines[0].startswith(">"):
        raise OptimizationError("input must be a single-record FASTA beginning with '>'")
    if any(line.startswith(">") for line in lines[1:]):
        raise OptimizationError("input FASTA must contain exactly one record")
    identifier_parts = lines[0][1:].strip().split(maxsplit=1)
    identifier = identifier_parts[0] if identifier_parts else ""
    if not identifier:
        raise OptimizationError("FASTA record identifier is empty")
    sequence = "".join(line.strip() for line in lines[1:])
    return identifier, sequence


def _write_fasta(path: Path, identifier: str, sequence: str) -> None:
    body = "\n".join(sequence[index : index + 70] for index in range(0, len(sequence), 70))
    path.write_text(f">{identifier}\n{body}\n", encoding="utf-8")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="codon-constrain",
        description="Exact synonymous codon optimization under finite-state DNA constraints.",
    )
    parser.add_argument("input", type=Path, help="single-record protein or coding-DNA FASTA")
    parser.add_argument("--output", type=Path, required=True, help="optimized construct FASTA")
    parser.add_argument("--report", type=Path, required=True, help="JSON evidence report")
    parser.add_argument("--host", choices=("ecoli", "human"), default="ecoli")
    parser.add_argument("--input-type", choices=("auto", "protein", "dna"), default="auto")
    parser.add_argument("--gc-min", type=float, default=0.0)
    parser.add_argument("--gc-max", type=float, default=1.0)
    parser.add_argument("--forbid", action="append", default=[], metavar="MOTIF")
    parser.add_argument("--homopolymer-max", type=int, default=99)
    parser.add_argument("--flank-5", default="")
    parser.add_argument("--flank-3", default="")
    parser.add_argument("--beam-width", type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the CLI, returning a process-compatible status code."""

    parser = _parser()
    args = parser.parse_args(argv)
    try:
        identifier, sequence = _read_fasta(args.input)
        constraints = Constraints(
            args.gc_min,
            args.gc_max,
            tuple(args.forbid),
            args.homopolymer_max,
            args.flank_5,
            args.flank_3,
        )
        result = optimize(
            sequence,
            args.host,
            input_type=args.input_type,
            constraints=constraints,
            beam_width=args.beam_width,
        )
        _write_fasta(args.output, f"{identifier}|codon-constrain|{args.host}", result.sequence)
        serialized = json.dumps(result.report, indent=2, sort_keys=True) + "\n"
        args.report.write_text(serialized, encoding="utf-8")
    except (OptimizationError, OSError) as error:
        print(f"codon-constrain: {error}", file=sys.stderr)
        return 2
    return 0

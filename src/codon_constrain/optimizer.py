"""Exact finite-state synonymous codon optimization."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal
from itertools import pairwise
from typing import Literal, TypeAlias

Host: TypeAlias = Literal["ecoli", "human"]
InputType: TypeAlias = Literal["auto", "protein", "dna"]


class OptimizationError(ValueError):
    """Raised when input is invalid or no feasible synonymous sequence exists."""


CODONS_BY_AA: dict[str, tuple[str, ...]] = {
    "A": ("GCA", "GCC", "GCG", "GCT"),
    "C": ("TGC", "TGT"),
    "D": ("GAC", "GAT"),
    "E": ("GAA", "GAG"),
    "F": ("TTC", "TTT"),
    "G": ("GGA", "GGC", "GGG", "GGT"),
    "H": ("CAC", "CAT"),
    "I": ("ATA", "ATC", "ATT"),
    "K": ("AAA", "AAG"),
    "L": ("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"),
    "M": ("ATG",),
    "N": ("AAC", "AAT"),
    "P": ("CCA", "CCC", "CCG", "CCT"),
    "Q": ("CAA", "CAG"),
    "R": ("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"),
    "S": ("AGC", "AGT", "TCA", "TCC", "TCG", "TCT"),
    "T": ("ACA", "ACC", "ACG", "ACT"),
    "V": ("GTA", "GTC", "GTG", "GTT"),
    "W": ("TGG",),
    "Y": ("TAC", "TAT"),
}

# Relative adaptiveness within each synonymous family. Values are normalized,
# bundled host profiles; only ratios within an amino-acid family affect scores.
HOST_WEIGHTS: dict[Host, dict[str, float]] = {
    "ecoli": {
        "GCA": 0.52,
        "GCC": 0.70,
        "GCG": 1.0,
        "GCT": 0.62,
        "TGC": 1.0,
        "TGT": 0.86,
        "GAC": 1.0,
        "GAT": 0.79,
        "GAA": 1.0,
        "GAG": 0.56,
        "TTC": 1.0,
        "TTT": 0.78,
        "GGA": 0.30,
        "GGC": 1.0,
        "GGG": 0.34,
        "GGT": 0.89,
        "CAC": 1.0,
        "CAT": 0.78,
        "ATA": 0.15,
        "ATC": 1.0,
        "ATT": 0.82,
        "AAA": 1.0,
        "AAG": 0.32,
        "CTA": 0.15,
        "CTC": 0.20,
        "CTG": 1.0,
        "CTT": 0.19,
        "TTA": 0.16,
        "TTG": 0.18,
        "ATG": 1.0,
        "AAC": 1.0,
        "AAT": 0.82,
        "CCA": 0.27,
        "CCC": 0.16,
        "CCG": 1.0,
        "CCT": 0.23,
        "CAA": 0.92,
        "CAG": 1.0,
        "AGA": 0.04,
        "AGG": 0.02,
        "CGA": 0.36,
        "CGC": 1.0,
        "CGG": 0.54,
        "CGT": 0.95,
        "AGC": 0.46,
        "AGT": 0.28,
        "TCA": 0.22,
        "TCC": 0.40,
        "TCG": 0.41,
        "TCT": 1.0,
        "ACA": 0.25,
        "ACC": 1.0,
        "ACG": 0.57,
        "ACT": 0.45,
        "GTA": 0.46,
        "GTC": 0.51,
        "GTG": 1.0,
        "GTT": 0.72,
        "TGG": 1.0,
        "TAC": 1.0,
        "TAT": 0.78,
    },
    "human": {
        "GCA": 0.61,
        "GCC": 1.0,
        "GCG": 0.28,
        "GCT": 0.67,
        "TGC": 1.0,
        "TGT": 0.85,
        "GAC": 1.0,
        "GAT": 0.89,
        "GAA": 1.0,
        "GAG": 0.96,
        "TTC": 1.0,
        "TTT": 0.86,
        "GGA": 0.59,
        "GGC": 1.0,
        "GGG": 0.73,
        "GGT": 0.52,
        "CAC": 1.0,
        "CAT": 0.72,
        "ATA": 0.49,
        "ATC": 1.0,
        "ATT": 0.80,
        "AAA": 1.0,
        "AAG": 0.96,
        "CTA": 0.40,
        "CTC": 0.81,
        "CTG": 1.0,
        "CTT": 0.52,
        "TTA": 0.20,
        "TTG": 0.52,
        "ATG": 1.0,
        "AAC": 1.0,
        "AAT": 0.85,
        "CCA": 0.66,
        "CCC": 1.0,
        "CCG": 0.36,
        "CCT": 0.70,
        "CAA": 0.36,
        "CAG": 1.0,
        "AGA": 1.0,
        "AGG": 0.99,
        "CGA": 0.39,
        "CGC": 0.86,
        "CGG": 0.89,
        "CGT": 0.38,
        "AGC": 1.0,
        "AGT": 0.60,
        "TCA": 0.58,
        "TCC": 0.82,
        "TCG": 0.24,
        "TCT": 0.74,
        "ACA": 0.71,
        "ACC": 1.0,
        "ACG": 0.33,
        "ACT": 0.66,
        "GTA": 0.39,
        "GTC": 0.72,
        "GTG": 1.0,
        "GTT": 0.58,
        "TGG": 1.0,
        "TAC": 1.0,
        "TAT": 0.82,
    },
}

CODON_TO_AA: dict[str, str] = {codon: aa for aa, codons in CODONS_BY_AA.items() for codon in codons}
_DNA: frozenset[str] = frozenset("ACGT")
_COMPLEMENT: Mapping[int, int | str | None] = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Constraints:
    """Finite-state constraints applied to the complete flanked construct."""

    gc_min: float = 0.0
    gc_max: float = 1.0
    forbidden_motifs: tuple[str, ...] = ()
    homopolymer_max: int = 99
    flank_5: str = ""
    flank_3: str = ""


@dataclass(frozen=True)
class OptimizationResult:
    """A feasible optimized construct and its machine-readable evidence."""

    protein: str
    coding_sequence: str
    sequence: str
    host: str
    log_score: float
    geometric_score: float
    gc_fraction: float
    exact: bool
    report: dict[str, object]


@dataclass(frozen=True)
class GreedyResult:
    """Unconstrained per-position greedy baseline and constraint assessment."""

    coding_sequence: str
    sequence: str
    feasible: bool
    violations: tuple[str, ...]
    log_score: float


# state: GC count, bounded motif suffix, last base, homopolymer run
_State: TypeAlias = tuple[int, str, str, int]
_Path: TypeAlias = tuple[float, str]
_CompletePath: TypeAlias = tuple[float, str, _State]


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of validated or unvalidated DNA text."""

    return sequence.translate(_COMPLEMENT)[::-1]


def translate_dna(sequence: str) -> str:
    """Translate a stop-free coding DNA sequence."""

    dna = sequence.upper()
    if not dna or len(dna) % 3:
        raise OptimizationError("coding DNA must be non-empty and have a length divisible by 3")
    invalid = sorted(set(dna) - _DNA)
    if invalid:
        raise OptimizationError(f"coding DNA contains invalid base(s): {''.join(invalid)}")
    amino_acids: list[str] = []
    for index in range(0, len(dna), 3):
        codon = dna[index : index + 3]
        amino_acid = CODON_TO_AA.get(codon)
        if amino_acid is None:
            codon_number = index // 3 + 1
            raise OptimizationError(
                f"coding DNA contains stop codon {codon} at codon {codon_number}"
            )
        amino_acids.append(amino_acid)
    return "".join(amino_acids)


def _protein_from_input(sequence: str, input_type: InputType) -> tuple[str, str | None]:
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise OptimizationError("input sequence is empty")
    if input_type not in {"auto", "protein", "dna"}:
        raise OptimizationError("input_type must be 'auto', 'protein', or 'dna'")
    use_dna = input_type == "dna" or (input_type == "auto" and set(cleaned) <= _DNA)
    if use_dna:
        return translate_dna(cleaned), cleaned
    invalid = sorted(set(cleaned) - set(CODONS_BY_AA))
    if invalid:
        raise OptimizationError(f"protein contains unsupported residue(s): {''.join(invalid)}")
    return cleaned, None


def _validate_dna(label: str, sequence: str, *, allow_empty: bool = True) -> str:
    dna = sequence.upper()
    if not dna and not allow_empty:
        raise OptimizationError(f"{label} must not be empty")
    invalid = sorted(set(dna) - _DNA)
    if invalid:
        raise OptimizationError(f"{label} contains invalid base(s): {''.join(invalid)}")
    return dna


def _normalized_constraints(constraints: Constraints) -> Constraints:
    if not 0.0 <= constraints.gc_min <= constraints.gc_max <= 1.0:
        raise OptimizationError("GC bounds must satisfy 0 <= gc_min <= gc_max <= 1")
    if constraints.homopolymer_max < 1:
        raise OptimizationError("homopolymer_max must be at least 1")
    motifs = tuple(
        _validate_dna("forbidden motif", motif, allow_empty=False)
        for motif in constraints.forbidden_motifs
    )
    return Constraints(
        constraints.gc_min,
        constraints.gc_max,
        motifs,
        constraints.homopolymer_max,
        _validate_dna("5' flank", constraints.flank_5),
        _validate_dna("3' flank", constraints.flank_3),
    )


def _expanded_motifs(motifs: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({item for motif in motifs for item in (motif, reverse_complement(motif))}))


def _advance(
    state: _State, segment: str, motifs: tuple[str, ...], suffix_limit: int, homopolymer_max: int
) -> _State | None:
    gc_count, suffix, last, run = state
    for base in segment:
        run = run + 1 if base == last else 1
        if run > homopolymer_max:
            return None
        candidate = suffix + base
        if any(candidate.endswith(motif) for motif in motifs):
            return None
        gc_count += base in "GC"
        suffix = candidate[-suffix_limit:] if suffix_limit else ""
        last = base
    return gc_count, suffix, last, run


def _gc_count_bounds(length: int, minimum: float, maximum: float) -> tuple[int, int]:
    low = int((Decimal(str(minimum)) * length).to_integral_value(rounding=ROUND_CEILING))
    high = int((Decimal(str(maximum)) * length).to_integral_value(rounding=ROUND_FLOOR))
    return low, high


def _violations(
    sequence: str, constraints: Constraints, motifs: tuple[str, ...]
) -> tuple[str, ...]:
    failures: list[str] = []
    gc = sum(base in "GC" for base in sequence) / len(sequence)
    if not constraints.gc_min <= gc <= constraints.gc_max:
        failures.append("GC fraction outside requested range")
    if any(motif in sequence for motif in motifs):
        failures.append("forbidden motif present on at least one strand")
    longest = 1
    run = 1
    for previous, current in pairwise(sequence):
        run = run + 1 if current == previous else 1
        longest = max(longest, run)
    if longest > constraints.homopolymer_max:
        failures.append("homopolymer limit exceeded")
    return tuple(failures)


def greedy_optimize(
    protein: str, host: Host = "ecoli", constraints: Constraints | None = None
) -> GreedyResult:
    """Build the unconstrained best-codon baseline, then assess constraints."""

    normalized = _normalized_constraints(constraints or Constraints())
    parsed_protein, _ = _protein_from_input(protein, "protein")
    if host not in HOST_WEIGHTS:
        raise OptimizationError(f"unknown host {host!r}; choose from: {', '.join(HOST_WEIGHTS)}")
    weights = HOST_WEIGHTS[host]
    chosen = [
        min(CODONS_BY_AA[aa], key=lambda codon: (-weights[codon], codon)) for aa in parsed_protein
    ]
    coding = "".join(chosen)
    sequence = normalized.flank_5 + coding + normalized.flank_3
    motifs = _expanded_motifs(normalized.forbidden_motifs)
    failures = _violations(sequence, normalized, motifs)
    score = sum(math.log(weights[codon]) for codon in chosen)
    return GreedyResult(coding, sequence, not failures, failures, score)


def optimize(
    sequence: str,
    host: Host = "ecoli",
    *,
    input_type: InputType = "auto",
    constraints: Constraints | None = None,
    beam_width: int | None = None,
) -> OptimizationResult:
    """Optimize synonymous codons with an exact finite-state dynamic program.

    With ``beam_width=None``, all reachable states are retained and the returned
    sequence is provably optimal for the represented GC, motif-history, and
    homopolymer states. A beam is deterministic but forfeits that proof.
    """

    protein, source_dna = _protein_from_input(sequence, input_type)
    normalized = _normalized_constraints(constraints or Constraints())
    if host not in HOST_WEIGHTS:
        raise OptimizationError(f"unknown host {host!r}; choose from: {', '.join(HOST_WEIGHTS)}")
    if beam_width is not None and beam_width < 1:
        raise OptimizationError("beam_width must be at least 1")

    motifs = _expanded_motifs(normalized.forbidden_motifs)
    suffix_limit = max((len(motif) - 1 for motif in motifs), default=0)
    initial: _State = (0, "", "", 0)
    flank_state = _advance(
        initial,
        normalized.flank_5,
        motifs,
        suffix_limit,
        normalized.homopolymer_max,
    )
    if flank_state is None:
        raise OptimizationError(
            "infeasible: 5' flank itself violates a motif or homopolymer constraint"
        )

    # Each state retains only the highest-scoring path; future feasibility and
    # score depend exclusively on this finite state, which establishes optimality.
    paths: dict[_State, _Path] = {flank_state: (0.0, "")}
    weights = HOST_WEIGHTS[host]
    for position, amino_acid in enumerate(protein, start=1):
        next_paths: dict[_State, _Path] = {}
        for state, (score, coding) in paths.items():
            for codon in CODONS_BY_AA[amino_acid]:
                advanced = _advance(state, codon, motifs, suffix_limit, normalized.homopolymer_max)
                if advanced is None:
                    continue
                candidate = (score + math.log(weights[codon]), coding + codon)
                previous = next_paths.get(advanced)
                if (
                    previous is None
                    or candidate[0] > previous[0]
                    or (candidate[0] == previous[0] and candidate[1] < previous[1])
                ):
                    next_paths[advanced] = candidate
        if not next_paths:
            raise OptimizationError(
                f"infeasible after residue {position} ({amino_acid}): all synonymous "
                "paths violate motif or homopolymer constraints"
            )
        if beam_width is not None and len(next_paths) > beam_width:
            ranked = sorted(next_paths.items(), key=lambda item: (-item[1][0], item[1][1], item[0]))
            next_paths = dict(ranked[:beam_width])
        paths = next_paths

    complete: list[_CompletePath] = []
    total_length = len(normalized.flank_5) + len(protein) * 3 + len(normalized.flank_3)
    gc_low, gc_high = _gc_count_bounds(total_length, normalized.gc_min, normalized.gc_max)
    terminal_constraint_failure = False
    for state, (score, coding) in paths.items():
        final_state = _advance(
            state,
            normalized.flank_3,
            motifs,
            suffix_limit,
            normalized.homopolymer_max,
        )
        if final_state is None:
            terminal_constraint_failure = True
            continue
        if gc_low <= final_state[0] <= gc_high:
            complete.append((score, coding, final_state))
    if not complete:
        reason = (
            "3' flank boundary violates motif or homopolymer constraints"
            if terminal_constraint_failure
            else "GC range is unreachable"
        )
        raise OptimizationError(
            f"infeasible: no synonymous sequence satisfies all constraints ({reason})"
        )

    score, coding, final_state = min(complete, key=lambda item: (-item[0], item[1]))
    construct = normalized.flank_5 + coding + normalized.flank_3
    gc_fraction = final_state[0] / len(construct)
    geometric = math.exp(score / len(protein))
    baseline = greedy_optimize(protein, host, normalized)
    constraint_report = {
        "gc_in_range": normalized.gc_min <= gc_fraction <= normalized.gc_max,
        "motif_free_both_strands": not any(motif in construct for motif in motifs),
        "homopolymer_ok": not any(
            base * (normalized.homopolymer_max + 1) in construct for base in "ACGT"
        ),
        "forbidden_motifs_checked": list(motifs),
    }
    report: dict[str, object] = {
        "host": host,
        "input_type": "dna" if source_dna is not None else "protein",
        "translation_preserved": translate_dna(coding) == protein,
        "protein_length": len(protein),
        "coding_length": len(coding),
        "construct_length": len(construct),
        "gc_fraction": gc_fraction,
        "geometric_relative_adaptiveness": geometric,
        "log_relative_adaptiveness": score,
        "exact_optimum": beam_width is None,
        "beam_width": beam_width,
        "constraints": constraint_report,
        "greedy_baseline": {
            "coding_sequence": baseline.coding_sequence,
            "feasible": baseline.feasible,
            "violations": list(baseline.violations),
            "log_relative_adaptiveness": baseline.log_score,
        },
    }
    return OptimizationResult(
        protein, coding, construct, host, score, geometric, gc_fraction, beam_width is None, report
    )

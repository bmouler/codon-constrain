from __future__ import annotations

import itertools
import math

import pytest

from codon_constrain import (
    Constraints,
    OptimizationError,
    greedy_optimize,
    optimize,
    reverse_complement,
    translate_dna,
)
from codon_constrain.optimizer import CODONS_BY_AA, HOST_WEIGHTS


def test_exact_matches_tiny_brute_force_oracle() -> None:
    protein = "KTL"
    constraints = Constraints(
        gc_min=0.44,
        gc_max=0.67,
        forbidden_motifs=("AAAA",),
        homopolymer_max=3,
    )
    result = optimize(protein, "ecoli", input_type="protein", constraints=constraints)

    feasible: list[tuple[float, str]] = []
    motifs = {"AAAA", reverse_complement("AAAA")}
    for codons in itertools.product(*(CODONS_BY_AA[aa] for aa in protein)):
        sequence = "".join(codons)
        gc = sum(base in "GC" for base in sequence) / len(sequence)
        motif_free = all(motif not in sequence for motif in motifs)
        homopolymer_ok = all(base * 4 not in sequence for base in "ACGT")
        if 0.44 <= gc <= 0.67 and motif_free and homopolymer_ok:
            score = sum(math.log(HOST_WEIGHTS["ecoli"][codon]) for codon in codons)
            feasible.append((score, sequence))
    expected_score, expected_sequence = min(feasible, key=lambda item: (-item[0], item[1]))
    assert result.coding_sequence == expected_sequence
    assert result.log_score == pytest.approx(expected_score)
    assert result.exact is True
    assert result.report["exact_optimum"] is True
    assert result.report["translation_preserved"] is True


def test_greedy_failure_dp_success_is_deterministic_capability_example() -> None:
    constraints = Constraints(homopolymer_max=3)
    greedy = greedy_optimize("KK", "ecoli", constraints)
    optimized = optimize("KK", "ecoli", input_type="protein", constraints=constraints)

    assert greedy.coding_sequence == "AAAAAA"
    assert not greedy.feasible
    assert greedy.violations == ("homopolymer limit exceeded",)
    assert optimized.coding_sequence == "AAGAAA"
    assert optimized.report["greedy_baseline"]["feasible"] is False
    assert optimized.report["constraints"] == {
        "gc_in_range": True,
        "motif_free_both_strands": True,
        "homopolymer_ok": True,
        "forbidden_motifs_checked": [],
    }


def test_greedy_reports_gc_and_both_strand_motif_violations() -> None:
    greedy = greedy_optimize(
        "K",
        constraints=Constraints(gc_min=1.0, gc_max=1.0, forbidden_motifs=("TTT",)),
    )
    assert greedy.violations == (
        "GC fraction outside requested range",
        "forbidden motif present on at least one strand",
    )


def test_both_hosts_have_distinct_preferences_and_beam_marks_approximate() -> None:
    ecoli = optimize("A", "ecoli", input_type="protein")
    human = optimize("A", "human", input_type="protein", beam_width=1)
    assert ecoli.coding_sequence == "GCG"
    assert human.coding_sequence == "GCC"
    assert human.exact is False
    assert human.report["beam_width"] == 1


def test_reverse_complement_motif_is_forbidden() -> None:
    # Forbidding TTT also forbids its reverse complement AAA, changing Lys AAA to AAG.
    result = optimize("K", constraints=Constraints(forbidden_motifs=("TTT",)), input_type="protein")
    assert result.coding_sequence == "AAG"
    assert result.report["constraints"]["forbidden_motifs_checked"] == ["AAA", "TTT"]


def test_coding_dna_translation_is_preserved_with_flanks() -> None:
    result = optimize(
        "AAGACTCTG",
        "human",
        input_type="dna",
        constraints=Constraints(flank_5="AC", flank_3="GT", gc_min=0.3, gc_max=0.8),
    )
    assert result.protein == "KTL"
    assert translate_dna(result.coding_sequence) == "KTL"
    assert result.sequence == "AC" + result.coding_sequence + "GT"
    assert result.report["input_type"] == "dna"
    assert result.gc_fraction == pytest.approx(
        sum(base in "GC" for base in result.sequence) / len(result.sequence)
    )
    assert result.geometric_score == pytest.approx(math.exp(result.log_score / len(result.protein)))


def test_auto_detects_dna_and_protein() -> None:
    assert optimize("ATG", input_type="auto").report["input_type"] == "dna"
    assert optimize("MW", input_type="auto").report["input_type"] == "protein"


def test_gc_and_boundary_motif_constraints() -> None:
    result = optimize(
        "K",
        constraints=Constraints(gc_min=0.4, gc_max=0.4, flank_5="C", flank_3="G"),
        input_type="protein",
    )
    assert result.sequence == "CAAAG"
    with pytest.raises(OptimizationError, match="3' flank boundary"):
        optimize(
            "M",
            constraints=Constraints(flank_3="G", forbidden_motifs=("TGG",)),
            input_type="protein",
        )


def test_infeasible_constraints_explain_stage_and_gc() -> None:
    with pytest.raises(OptimizationError, match=r"residue 1 .*all synonymous paths"):
        optimize("M", constraints=Constraints(forbidden_motifs=("ATG",)), input_type="protein")
    with pytest.raises(OptimizationError, match="GC range is unreachable"):
        optimize("M", constraints=Constraints(gc_min=1.0, gc_max=1.0), input_type="protein")
    with pytest.raises(OptimizationError, match="5' flank itself"):
        optimize(
            "M",
            constraints=Constraints(flank_5="AAAA", homopolymer_max=3),
            input_type="protein",
        )


@pytest.mark.parametrize(
    ("sequence", "input_type", "message"),
    [
        ("", "protein", "empty"),
        ("MX", "protein", "unsupported residue"),
        ("AT", "dna", "divisible by 3"),
        ("ATX", "dna", "invalid base"),
        ("TAA", "dna", "stop codon TAA"),
        ("M", "other", "input_type must"),
    ],
)
def test_invalid_sequence_input(sequence: str, input_type: str, message: str) -> None:
    with pytest.raises(OptimizationError, match=message):
        optimize(sequence, input_type=input_type)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"constraints": Constraints(gc_min=-0.1)}, "GC bounds"),
        ({"constraints": Constraints(gc_min=0.8, gc_max=0.2)}, "GC bounds"),
        ({"constraints": Constraints(homopolymer_max=0)}, "at least 1"),
        ({"constraints": Constraints(forbidden_motifs=("AX",))}, "forbidden motif contains"),
        ({"constraints": Constraints(forbidden_motifs=("",))}, "must not be empty"),
        ({"constraints": Constraints(flank_5="N")}, "5' flank contains"),
        ({"constraints": Constraints(flank_3="N")}, "3' flank contains"),
        ({"host": "mouse"}, "unknown host"),
        ({"beam_width": 0}, "beam_width"),
    ],
)
def test_invalid_configuration(kwargs: dict[str, object], message: str) -> None:
    with pytest.raises(OptimizationError, match=message):
        optimize("M", input_type="protein", **kwargs)


def test_translate_and_reverse_complement_helpers() -> None:
    assert reverse_complement("ACGT") == "ACGT"
    assert translate_dna("atgggc") == "MG"
    with pytest.raises(OptimizationError, match="non-empty"):
        translate_dna("")


def test_greedy_validates_host_and_protein() -> None:
    with pytest.raises(OptimizationError, match="unknown host"):
        greedy_optimize("M", "mouse")
    with pytest.raises(OptimizationError, match="unsupported residue"):
        greedy_optimize("*")

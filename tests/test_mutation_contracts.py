from __future__ import annotations

import itertools
import math
from pathlib import Path

import pytest

from codon_constrain import Constraints, OptimizationError, greedy_optimize, optimize, translate_dna
from codon_constrain.cli import _parser, _read_fasta, _write_fasta
from codon_constrain.optimizer import (
    _advance,
    _gc_count_bounds,
    _normalized_constraints,
    _protein_from_input,
    _violations,
)


def test_parser_contract_and_defaults(tmp_path: Path) -> None:
    parser = _parser()
    assert parser.prog == "codon-constrain"
    assert parser.description == (
        "Exact synonymous codon optimization under finite-state DNA constraints."
    )
    args = parser.parse_args(["input.fa", "--output", "out.fa", "--report", "report.json"])
    assert vars(args) == {
        "input": Path("input.fa"),
        "output": Path("out.fa"),
        "report": Path("report.json"),
        "host": "ecoli",
        "input_type": "auto",
        "gc_min": 0.0,
        "gc_max": 1.0,
        "forbid": [],
        "homopolymer_max": 99,
        "flank_5": "",
        "flank_3": "",
        "beam_width": None,
    }
    assert parser.format_help() == (
        "usage: codon-constrain [-h] --output OUTPUT --report REPORT\n"
        "                       [--host {ecoli,human}]\n"
        "                       [--input-type {auto,protein,dna}] [--gc-min GC_MIN]\n"
        "                       [--gc-max GC_MAX] [--forbid MOTIF]\n"
        "                       [--homopolymer-max HOMOPOLYMER_MAX] [--flank-5 FLANK_5]\n"
        "                       [--flank-3 FLANK_3] [--beam-width BEAM_WIDTH]\n"
        "                       input\n"
        "\n"
        "Exact synonymous codon optimization under finite-state DNA constraints.\n"
        "\n"
        "positional arguments:\n"
        "  input                 single-record protein or coding-DNA FASTA\n"
        "\n"
        "options:\n"
        "  -h, --help            show this help message and exit\n"
        "  --output OUTPUT       optimized construct FASTA\n"
        "  --report REPORT       JSON evidence report\n"
        "  --host {ecoli,human}\n"
        "  --input-type {auto,protein,dna}\n"
        "  --gc-min GC_MIN\n"
        "  --gc-max GC_MAX\n"
        "  --forbid MOTIF\n"
        "  --homopolymer-max HOMOPOLYMER_MAX\n"
        "  --flank-5 FLANK_5\n"
        "  --flank-3 FLANK_3\n"
        "  --beam-width BEAM_WIDTH\n"
    )


def test_fasta_contract_exact_io(tmp_path: Path) -> None:
    source = tmp_path / "source.fa"
    source.write_text(">record description\n acg \n tta\n", encoding="utf-8")
    assert _read_fasta(source) == ("record", "acgtta")
    source.write_text(">record two descriptions\nACG\n", encoding="utf-8")
    assert _read_fasta(source) == ("record", "ACG")

    output = tmp_path / "output.fa"
    _write_fasta(output, "id", "A" * 141)
    assert output.read_text(encoding="utf-8") == f">id\n{'A' * 70}\n{'A' * 70}\nA\n"

    errors = [
        ("", "input must be a single-record FASTA beginning with '>'"),
        ("M\n", "input must be a single-record FASTA beginning with '>'"),
        (">one\n>two\nM\n", "input FASTA must contain exactly one record"),
        (">   \nM\n", "FASTA record identifier is empty"),
    ]
    for index, (content, message) in enumerate(errors):
        invalid = tmp_path / f"invalid-{index}.fa"
        invalid.write_text(content, encoding="utf-8")
        with pytest.raises(OptimizationError) as error:
            _read_fasta(invalid)
        assert str(error.value) == message


def test_translation_and_input_contract_exact_errors() -> None:
    assert translate_dna("atggcttgg") == "MAW"
    cases = [
        ("", "coding DNA must be non-empty and have a length divisible by 3"),
        ("AT", "coding DNA must be non-empty and have a length divisible by 3"),
        ("AXZ", "coding DNA contains invalid base(s): XZ"),
        ("ATGTAG", "coding DNA contains stop codon TAG at codon 2"),
    ]
    for sequence, message in cases:
        with pytest.raises(OptimizationError) as error:
            translate_dna(sequence)
        assert str(error.value) == message

    assert _protein_from_input(" a t g ", "auto") == ("M", "ATG")
    assert _protein_from_input(" mw ", "auto") == ("MW", None)
    assert _protein_from_input("acg", "protein") == ("ACG", None)
    with pytest.raises(OptimizationError) as error:
        _protein_from_input("", "auto")
    assert str(error.value) == "input sequence is empty"
    with pytest.raises(OptimizationError) as error:
        _protein_from_input("M", "wrong")  # type: ignore[arg-type]
    assert str(error.value) == "input_type must be 'auto', 'protein', or 'dna'"
    with pytest.raises(OptimizationError) as error:
        _protein_from_input("MX*", "protein")
    assert str(error.value) == "protein contains unsupported residue(s): *X"


def test_constraint_helper_contracts() -> None:
    normalized = _normalized_constraints(Constraints(0.25, 0.75, ("acg",), 2, "at", "gc"))
    assert normalized == Constraints(0.25, 0.75, ("ACG",), 2, "AT", "GC")
    assert _gc_count_bounds(5, 0.21, 0.79) == (2, 3)

    assert _advance((0, "", "", 0), "AGC", ("TAG",), 2, 2) == (2, "GC", "C", 1)
    assert _advance((0, "TA", "A", 1), "G", ("TAG",), 2, 3) is None
    assert _advance((0, "", "A", 2), "A", (), 0, 2) is None

    constraints = Constraints(0.5, 0.5, ("TAG",), 2)
    assert _violations("AAATAG", constraints, ("TAG",)) == (
        "GC fraction outside requested range",
        "forbidden motif present on at least one strand",
        "homopolymer limit exceeded",
    )
    assert _violations("ATGC", Constraints(), ()) == ()


def test_optimization_result_and_report_are_complete() -> None:
    result = optimize(
        "KK",
        "ecoli",
        input_type="protein",
        constraints=Constraints(homopolymer_max=3, flank_5="C", flank_3="G"),
    )
    assert result.protein == "KK"
    assert result.coding_sequence == "AAGAAA"
    assert result.sequence == "CAAGAAAG"
    assert result.host == "ecoli"
    assert result.log_score == pytest.approx(math.log(0.32))
    assert result.geometric_score == pytest.approx(math.sqrt(0.32))
    assert result.gc_fraction == pytest.approx(3 / 8)
    assert result.exact is True
    assert result.report == {
        "host": "ecoli",
        "input_type": "protein",
        "translation_preserved": True,
        "protein_length": 2,
        "coding_length": 6,
        "construct_length": 8,
        "gc_fraction": pytest.approx(3 / 8),
        "geometric_relative_adaptiveness": pytest.approx(math.sqrt(0.32)),
        "log_relative_adaptiveness": pytest.approx(math.log(0.32)),
        "exact_optimum": True,
        "beam_width": None,
        "constraints": {
            "gc_in_range": True,
            "motif_free_both_strands": True,
            "homopolymer_ok": True,
            "forbidden_motifs_checked": [],
        },
        "greedy_baseline": {
            "coding_sequence": "AAAAAA",
            "feasible": False,
            "violations": ["homopolymer limit exceeded"],
            "log_relative_adaptiveness": 0.0,
        },
    }


def test_greedy_result_contract() -> None:
    result = greedy_optimize("AR", "human", Constraints(flank_5="T", flank_3="A"))
    assert result.coding_sequence == "GCCAGA"
    assert result.sequence == "TGCCAGAA"
    assert result.feasible is True
    assert result.violations == ()
    assert result.log_score == 0.0


def test_optimizer_default_auto_and_constraint_report_boundaries() -> None:
    auto = optimize("ATG")
    assert auto.report["input_type"] == "dna"
    exact_gc = optimize(
        "M",
        input_type="protein",
        constraints=Constraints(gc_min=0.333, gc_max=0.334),
    )
    assert exact_gc.report["constraints"] == {
        "gc_in_range": True,
        "motif_free_both_strands": True,
        "homopolymer_ok": True,
        "forbidden_motifs_checked": [],
    }
    motif = optimize(
        "K",
        input_type="protein",
        constraints=Constraints(forbidden_motifs=("TTT",), homopolymer_max=3),
    )
    assert motif.report["constraints"] == {
        "gc_in_range": True,
        "motif_free_both_strands": True,
        "homopolymer_ok": True,
        "forbidden_motifs_checked": ["AAA", "TTT"],
    }


def test_cli_serialization_and_argument_forwarding(tmp_path: Path) -> None:
    from codon_constrain.cli import main

    source = tmp_path / "source.fa"
    output = tmp_path / "output.fa"
    report = tmp_path / "report.json"
    source.write_text(">sample\nGCT\n", encoding="utf-8")
    status = main(
        [
            str(source),
            "--output",
            str(output),
            "--report",
            str(report),
            "--host",
            "human",
            "--input-type",
            "dna",
            "--gc-min",
            "0.2",
            "--gc-max",
            "1.0",
            "--forbid",
            "AAAA",
            "--homopolymer-max",
            "3",
            "--flank-5",
            "C",
            "--flank-3",
            "G",
            "--beam-width",
            "1",
        ]
    )
    assert status == 0
    assert output.read_text(encoding="utf-8") == ">sample|codon-constrain|human\nCGCCG\n"
    parsed = optimize(
        "GCT",
        "human",
        input_type="dna",
        constraints=Constraints(0.2, 1.0, ("AAAA",), 3, "C", "G"),
        beam_width=1,
    ).report
    assert (
        report.read_text(encoding="utf-8")
        == __import__("json").dumps(parsed, indent=2, sort_keys=True) + "\n"
    )


def test_exact_optimizer_matches_exhaustive_oracle_across_state_collisions() -> None:
    from codon_constrain.optimizer import CODONS_BY_AA, HOST_WEIGHTS

    for amino_acids in itertools.product("AKLRS", repeat=3):
        protein = "".join(amino_acids)
        candidates = []
        for codons in itertools.product(*(CODONS_BY_AA[amino_acid] for amino_acid in protein)):
            sequence = "".join(codons)
            score = sum(math.log(HOST_WEIGHTS["ecoli"][codon]) for codon in codons)
            candidates.append((score, sequence))
        expected_score, expected_sequence = min(
            candidates, key=lambda candidate: (-candidate[0], candidate[1])
        )
        actual = optimize(protein, input_type="protein")
        assert actual.coding_sequence == expected_sequence
        assert actual.log_score == pytest.approx(expected_score)


def test_boundaries_and_exact_diagnostics() -> None:
    assert _protein_from_input("ACGTACGTACGT", "auto") == ("TYVR", "ACGTACGTACGT")
    for constraints, message in (
        (Constraints(gc_max=1.1), "GC bounds must satisfy 0 <= gc_min <= gc_max <= 1"),
        (Constraints(homopolymer_max=0), "homopolymer_max must be at least 1"),
        (Constraints(forbidden_motifs=("",)), "forbidden motif must not be empty"),
    ):
        with pytest.raises(OptimizationError) as error:
            _normalized_constraints(constraints)
        assert str(error.value) == message
    assert _normalized_constraints(Constraints(homopolymer_max=1)).homopolymer_max == 1

    with pytest.raises(OptimizationError) as error:
        optimize("M", "mouse", input_type="protein")  # type: ignore[arg-type]
    assert str(error.value) == "unknown host 'mouse'; choose from: ecoli, human"
    with pytest.raises(OptimizationError) as error:
        optimize("M", input_type="protein", beam_width=0)
    assert str(error.value) == "beam_width must be at least 1"
    with pytest.raises(OptimizationError) as error:
        optimize(
            "M",
            input_type="protein",
            constraints=Constraints(flank_5="AAAA", homopolymer_max=3),
        )
    assert str(error.value) == (
        "infeasible: 5' flank itself violates a motif or homopolymer constraint"
    )
    with pytest.raises(OptimizationError) as error:
        optimize("M", input_type="protein", constraints=Constraints(gc_min=1.0))
    assert str(error.value) == (
        "infeasible: no synonymous sequence satisfies all constraints (GC range is unreachable)"
    )


def test_violation_boundaries_and_report_homopolymers() -> None:
    assert _violations("GC", Constraints(gc_min=1.0, gc_max=1.0), ()) == ()
    assert _violations("AG", Constraints(homopolymer_max=1), ()) == ()
    assert _violations("AA", Constraints(homopolymer_max=1), ()) == ("homopolymer limit exceeded",)
    assert _violations("AAT", Constraints(homopolymer_max=2), ()) == ()
    result = optimize(
        "KK",
        input_type="protein",
        constraints=Constraints(homopolymer_max=3),
    )
    assert result.report["constraints"]["homopolymer_ok"] is True


def test_required_cli_outputs_are_enforced() -> None:
    parser = _parser()
    for argv, missing in (
        (["in.fa", "--report", "r.json"], "--output"),
        (["in.fa", "--output", "out.fa"], "--report"),
    ):
        with pytest.raises(SystemExit) as error:
            parser.parse_args(argv)
        assert error.value.code == 2
        assert missing in parser.format_usage()


def test_cli_auto_type_is_forwarded(tmp_path: Path) -> None:
    from codon_constrain.cli import main

    source = tmp_path / "source.fa"
    output = tmp_path / "out.fa"
    report = tmp_path / "report.json"
    source.write_text(">protein\nACG\n", encoding="utf-8")
    assert (
        main(
            [
                str(source),
                "--output",
                str(output),
                "--report",
                str(report),
                "--input-type",
                "protein",
            ]
        )
        == 0
    )
    assert '"input_type": "protein"' in report.read_text(encoding="utf-8")


def test_suffix_history_and_beam_boundaries_are_behavioral() -> None:
    boundary = optimize("A", input_type="protein", beam_width=4)
    unbounded = optimize("A", input_type="protein")
    assert boundary.coding_sequence == unbounded.coding_sequence
    assert boundary.log_score == unbounded.log_score

    motif = optimize(
        "KK",
        input_type="protein",
        constraints=Constraints(forbidden_motifs=("AAGA",)),
    )
    assert "AAGA" not in motif.sequence
    assert "TCTT" not in motif.sequence


@pytest.mark.parametrize(
    "base",
    ["A", "C", "G", "T"],
)
def test_homopolymer_report_checks_every_base(base: str) -> None:
    result = optimize(
        {"A": "K", "C": "P", "G": "G", "T": "F"}[base],
        input_type="protein",
        constraints=Constraints(homopolymer_max=3),
    )
    assert result.report["constraints"]["homopolymer_ok"] is True


def test_gc_report_includes_both_endpoints() -> None:
    lower = optimize("M", input_type="protein", constraints=Constraints(gc_min=1 / 3 - 1e-9))
    upper = optimize("M", input_type="protein", constraints=Constraints(gc_max=1 / 3 + 1e-9))
    assert lower.report["constraints"]["gc_in_range"] is True
    assert upper.report["constraints"]["gc_in_range"] is True

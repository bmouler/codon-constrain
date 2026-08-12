from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from codon_constrain.cli import _read_fasta, main
from codon_constrain.optimizer import OptimizationError


def test_installed_module_cli_end_to_end_from_clean_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path = tmp_path / "input.faa"
    output_path = tmp_path / "optimized.fasta"
    report_path = tmp_path / "report.json"
    input_path.write_text(">enzyme representative\nKK\n", encoding="utf-8")
    command = Path(sys.executable).with_name("codon-constrain")

    completed = subprocess.run(
        [
            str(command),
            str(input_path),
            "--output",
            str(output_path),
            "--report",
            str(report_path),
            "--input-type",
            "protein",
            "--host",
            "ecoli",
            "--homopolymer-max",
            "3",
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert output_path.read_text(encoding="utf-8") == ">enzyme|codon-constrain|ecoli\nAAGAAA\n"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    # The subprocess covers __main__.py externally; repeat its two-line dispatch
    # in-process so coverage also observes the packaged entry point.
    monkeypatch.setattr("codon_constrain.cli.main", lambda: 0)
    namespace: dict[str, object] = {
        "__name__": "codon_constrain.__main__",
        "__package__": "codon_constrain",
    }
    entrypoint = Path(__file__).parents[1] / "src" / "codon_constrain" / "__main__.py"
    with pytest.raises(SystemExit) as exit_info:
        exec(compile(entrypoint.read_text(encoding="utf-8"), str(entrypoint), "exec"), namespace)
    assert exit_info.value.code == 0
    assert report["translation_preserved"] is True
    assert report["greedy_baseline"]["violations"] == ["homopolymer limit exceeded"]
    assert report["exact_optimum"] is True


def test_cli_writes_wrapped_fasta_and_accepts_every_option(tmp_path: Path) -> None:
    source = tmp_path / "source.faa"
    output = tmp_path / "out.fa"
    report = tmp_path / "out.json"
    source.write_text(">long\n" + "M" * 25 + "\n", encoding="utf-8")
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
            "protein",
            "--gc-min",
            "0.3",
            "--gc-max",
            "0.4",
            "--forbid",
            "AAAA",
            "--homopolymer-max",
            "9",
            "--flank-5",
            "C",
            "--flank-3",
            "G",
            "--beam-width",
            "2",
        ]
    )
    assert status == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert len(lines[1]) == 70
    assert len(lines[2]) == 7


def test_cli_reports_domain_and_output_errors(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    invalid = tmp_path / "invalid.fa"
    invalid.write_text(">bad\nX\n", encoding="utf-8")
    directory_output = tmp_path / "directory"
    directory_output.mkdir()
    status = main(
        [
            str(invalid),
            "--output",
            str(directory_output),
            "--report",
            str(tmp_path / "r.json"),
            "--input-type",
            "protein",
        ]
    )
    assert status == 2
    assert "unsupported residue" in capsys.readouterr().err

    valid = tmp_path / "valid.fa"
    valid.write_text(">ok\nM\n", encoding="utf-8")
    status = main(
        [
            str(valid),
            "--output",
            str(directory_output),
            "--report",
            str(tmp_path / "r.json"),
            "--input-type",
            "protein",
        ]
    )
    assert status == 2
    assert "Is a directory" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("M\n", "beginning"),
        (">one\nM\n>two\nM\n", "exactly one"),
        (">   \nM\n", "identifier is empty"),
    ],
)
def test_fasta_validation(tmp_path: Path, content: str, message: str) -> None:
    path = tmp_path / "input.fa"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(OptimizationError, match=message):
        _read_fasta(path)


def test_missing_fasta_explains_io_error(tmp_path: Path) -> None:
    with pytest.raises(OptimizationError, match="cannot read input FASTA"):
        _read_fasta(tmp_path / "missing.fa")

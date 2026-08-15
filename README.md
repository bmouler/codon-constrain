# codon-constrain

[![CI](https://github.com/bmouler/codon-constrain/actions/workflows/ci.yml/badge.svg)](https://github.com/bmouler/codon-constrain/actions/workflows/ci.yml)
![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)
![Types](https://img.shields.io/badge/types-mypy%20strict-blue)
![Mutation](https://img.shields.io/badge/mutation-95%25%20killed-brightgreen)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

`codon-constrain` is a dependency-free-runtime synonymous codon optimizer. It maximizes the sum of log relative codon adaptiveness for bundled *E. coli* or human profiles while preserving translation and jointly enforcing whole-construct GC, forbidden motifs on either strand, homopolymer length, and fixed DNA flanks.

## Install

```bash
python -m pip install codon-constrain
```

For development:

```bash
python -m pip install -e '.[dev]'
pytest --cov=codon_constrain --cov-branch --cov-fail-under=100
ruff check .
```

## Quickstart

Given `enzyme.faa`:

```text
>enzyme
KK
```

run:

```bash
codon-constrain enzyme.faa \
  --output optimized.fasta \
  --report report.json \
  --host ecoli \
  --input-type protein \
  --gc-min 0.30 --gc-max 0.80 \
  --forbid GAATTC \
  --homopolymer-max 3
```

The FASTA contains the complete flanked construct. The JSON records translation preservation, the CAI-like geometric mean of relative adaptiveness, GC fraction, constraint checks, whether the solution is exact, and the unconstrained greedy baseline.

Python API:

```python
from codon_constrain import Constraints, optimize

result = optimize(
    "KK",
    host="ecoli",
    input_type="protein",
    constraints=Constraints(homopolymer_max=3),
)
assert result.coding_sequence == "AAGAAA"
assert result.report["translation_preserved"]
```

Coding-DNA FASTA is also accepted with `--input-type dna` (or `auto` when unambiguous). Stop codons and unsupported residues are rejected with position-aware messages.

## Algorithm and exactness

```mermaid
flowchart LR; A[protein or coding DNA] --> S[DP states: GC count × suffix × homopolymer run]; S --> T[extend by synonymous codons]; T --> C{constraints: GC, motifs, flanks}; C -->|feasible| B[max log adaptiveness]; B --> R[traceback: optimal construct + report]
```

At each amino-acid position, dynamic programming extends every reachable finite state by one synonymous codon. A state consists of:

- total GC count;
- the suffix needed to recognize a forbidden motif at the next base;
- the final base and current homopolymer run length.

Motifs are expanded to include their reverse complements. Fixed 5' and 3' flanks pass through the same state machine, so violations spanning a flank/coding boundary are detected. For identical states, only the highest-scoring prefix is retained because all future feasibility and score contributions depend solely on that state. Therefore the default unpruned mode proves the maximum score over these represented finite-state constraints; this is dynamic programming, not full-sequence enumeration. `--beam-width N` retains only the best `N` states per residue to cap memory, and the report honestly marks that result non-exact.

The objective is

```text
sum(log(relative_adaptiveness(codon)))
```

and its reported CAI-like value is the geometric mean, `exp(log_score / codon_count)`. Bundled weights are normalized within each synonymous family; they are useful deterministic host profiles rather than an expression guarantee.

## Reproducible capability evidence

The test suite contains an exhaustive tiny oracle for protein `KTL`: it enumerates only that test's 36 synonymous sequences and verifies that exact DP returns the same feasible optimum under GC, motif, and homopolymer constraints. This independently checks the optimality recurrence.

A deterministic material improvement over greedy selection is exercised end to end:

```text
protein:             KK
host:                E. coli
greedy best codons:  AAAAAA   (violates homopolymer_max=3)
exact DP result:     AAGAAA   (translation KK; feasible)
```

Run the evidence:

```bash
pytest tests/test_optimizer.py::test_exact_matches_tiny_brute_force_oracle \
       tests/test_optimizer.py::test_greedy_failure_dp_success_is_deterministic_capability_example \
       tests/test_cli.py::test_installed_module_cli_end_to_end_from_clean_directory
```

### End-to-end performance

`python benchmarks/benchmark_optimize.py --source-root . --samples 15 --warmups 3` runs the
public exact optimizer for a 104-residue protein under a global GC range, four forbidden motifs,
a homopolymer bound, and fixed flanks, then materializes `OptimizationResult` and serializes its
complete dataclass payload for the checksum.

On an Apple M3 Max with CPython 3.11.12 on 2026-08-15, frozen baseline `07a9b99fb702`
measured **297.357 ms** median and the cached finite-state implementation **47.964 ms**, a
**6.200x speedup** over 15 samples after three warmups. Both runs produced SHA-256
`8e0e474a503df6747754560baa0cb9fb92b18ab09ece8b89f226cec781a2d7ff`. Fixture
construction and interpreter startup are outside the timed region; exact search, constraint
reporting, traceback, result materialization, and conservative checksum serialization are
included. Use `--source-root` to point the
same script at another worktree for a direct comparison.

### Mutation testing

From the repository root, reproduce the mutation run with:

```bash
source .venv/bin/activate
mutmut run
mutmut results
```

The run generated 807 mutants and killed 768 (95.17%). The remaining 39 were individually reviewed as behavior-equivalent under the supported contract, not missed mutants. There were zero suspicious results and zero timeouts. Observable CLI help changes were killed rather than classified as equivalent.

| Reviewed equivalent rationale | Count |
|---|---:|
| Explicit UTF-8 versus the platform default, and type-cast identities | 9 |
| FASTA header tokenization with an equivalent `maxsplit` | 1 |
| Sentinels made redundant by validated DNA and motif alphabets | 6 |
| Surplus suffix/initial state retained by the DP | 7 |
| Deterministic tie-breaking and beam/pruning boundaries | 8 |
| Reporting, default-value, and formatting identities | 8 |
| **Total reviewed equivalents** | **39** |

## Limitations

- Optimization models synonymous codon choice only; it does not predict expression, mRNA folding, ribosome traffic, splicing, synthesis success, or biological safety.
- Inputs must use the 20 canonical amino acids or stop-free `ACGT` coding DNA. Ambiguous bases, stop codons, and selenocysteine are unsupported.
- GC is enforced globally across the complete flanked construct, not in sliding windows.
- Long motifs and broad GC ranges can produce many states. Beam search bounds memory but removes the optimality proof.
- Host tables are built-in relative profiles and are not tissue-, strain-, condition-, or gene-specific. Validate designs experimentally and against current domain-specific requirements.

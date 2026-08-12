from hypothesis import given, settings
from hypothesis import strategies as st

from codon_constrain import greedy_optimize, optimize, reverse_complement, translate_dna

proteins = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=12)


@settings(max_examples=50)
@given(protein=proteins)
def test_optimized_coding_sequence_translates_to_input_protein(protein: str) -> None:
    result = optimize(protein, host="ecoli", input_type="protein")

    assert translate_dna(result.coding_sequence) == protein


@settings(max_examples=50)
@given(protein=proteins)
def test_constrained_optimum_never_beats_unconstrained_greedy(protein: str) -> None:
    result = optimize(protein, host="ecoli", input_type="protein")
    greedy = greedy_optimize(protein, host="ecoli")

    assert result.log_score <= greedy.log_score + 1e-9


@given(dna=st.text(alphabet="ACGT"))
def test_reverse_complement_is_an_involution(dna: str) -> None:
    assert reverse_complement(reverse_complement(dna)) == dna

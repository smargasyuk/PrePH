from dataclasses import dataclass
import functools
import pytest
import numpy as np
from preph.fold import Index_seq, FindMinEnLocAlkmer
from preph.asset_loader import load_kmer_table


def test_e2e_readme():
    k = 3
    GT_threshold = 2
    seq = 'AAAGGGC'
    seq_compl = 'AAAGCCCAAAAAACCTTT'
    energy_threshold = -1
    handle_length_threshold = 0
    need_suboptimal = True

    kmers_stacking_matrix = load_kmer_table(k, GT_threshold)
    seq_indxd = Index_seq(seq.encode("ascii"), k)
    seq_compl_indxd = Index_seq(seq_compl.encode("ascii"), k)
    res = FindMinEnLocAlkmer(seq.encode("ascii"), seq_compl.encode("ascii"), 
                            seq_indxd, seq_compl_indxd, k, energy_threshold, handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
    
    assert len(res) == 2
    assert res[0][0] == -10
    assert res[1][0] == -7.2


def test_bench_FindMinEnLocAlkmer(benchmark):
    k = 5
    GT_threshold = 2
    kmers_stacking_matrix = load_kmer_table(k, GT_threshold)
    
    @benchmark
    def e2e_bench():
        seq = 'GTAGAAAAGAAAAATGACAGAGACCAACAGGAACTGAATTGTTTAGAGTGTAGTTTGAAGCTTTCAAAGGTTGTTCTCAGCTAAACTTCAGAACTGACAAAAAGTATGAGTGTCTCTTTTATTCCATAATGTTTTATATTATCCTGAAAAAAACTACCCTTTGGCCTTCAATGAAGCCTAGAATATTATTGCCATCATATTAGTCTTGCTAGACAATTTATAGTTTTTTATTATTTTATCTTTTAG'
        seq_compl = 'GGTAGAGTAGAAGAAAAAGATATAAACAGGAAGGAAGTACCCAGGTTTTATAAATCCCAACAACTGGTAATTTAGAATGAGGGATTTTGGAGCTAACCTAAGAATATAGTGGCTTTTTTCTGATGGAGTCTTTCTCTGTCGCCCAGGCTGGAGTGCAGTGGCACAATCTCGACTCATTGCAACCTCTGCCTCCTGGGTTCAAACGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACAGGC'
        energy_threshold = -1
        handle_length_threshold = 0
        need_suboptimal = True

        seq_indxd = Index_seq(seq.encode("ascii"), k)
        seq_compl_indxd = Index_seq(seq_compl.encode("ascii"), k)
        res = FindMinEnLocAlkmer(seq.encode("ascii"), seq_compl.encode("ascii"), 
                                seq_indxd, seq_compl_indxd, k, energy_threshold, handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
 
        assert len(res) == 156
        assert res[0][0] == -18.3
        assert res[1][0] == -17.1 


@functools.lru_cache(maxsize=None)
def _stacking_matrix(k: int, GT_threshold: int) -> np.ndarray:
    return load_kmer_table(k, GT_threshold)


def find_min_energy(seq: str, seq_compl: str, k: int = 3, GT_threshold: int = 2,
                    energy_threshold: float = -1, handle_length_threshold: float = 0,
                    need_suboptimal: bool = True):
    """Run the full preparation pipeline (stacking matrix, kmer indexing) and
    return the list of PCCRs found. FindMinEnLocAlkmer returns 0 when nothing
    is found; this helper normalizes that to an empty list."""
    kmers_stacking_matrix = _stacking_matrix(k, GT_threshold)
    seq_indxd = Index_seq(seq.encode("ascii"), k)
    seq_compl_indxd = Index_seq(seq_compl.encode("ascii"), k)
    res = FindMinEnLocAlkmer(seq.encode("ascii"), seq_compl.encode("ascii"),
                            seq_indxd, seq_compl_indxd, k, energy_threshold,
                            handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
    return res if isinstance(res, list) else []


@dataclass(frozen=True)
class E2eCase:
    name: str
    seq: str
    seq_compl: str
    k: int
    energy_threshold: float
    handle_length_threshold: float
    need_suboptimal: bool
    expected_n: int
    expected_energies: tuple  # energies of all returned PCCRs, best first
    expected_best: tuple  # (start_j, end_j, end_i-k-2, start_i-k-2, align1, align2, structure) of the best PCCR


def _case(name, seq, seq_compl, k=3, energy_threshold=-1, handle_length_threshold=0,
          need_suboptimal=True, expected_n=None, expected_energies=None, expected_best=None):
    return E2eCase(name, seq, seq_compl, k, energy_threshold, handle_length_threshold,
                   need_suboptimal, expected_n, expected_energies, expected_best)


CASES = [
    # The original example from the README.
    _case("baseline", "AAAGGGC", "AAAGCCCAAAAAACCTTT",
          expected_n=2,
          expected_energies=(-10.0, -7.2),
          expected_best=(3, 6, 3, 6, b"GGGC", b"GCCC", b"(((())))")),
    # need_suboptimal=False must return only the optimal structure.
    _case("no_suboptimal", "AAAGGGC", "AAAGCCCAAAAAACCTTT", need_suboptimal=False,
          expected_n=1,
          expected_energies=(-10.0,),
          expected_best=(3, 6, 3, 6, b"GGGC", b"GCCC", b"(((())))")),
    # energy_threshold=-9.9 filters out the -7.2 suboptimal structure.
    _case("energy_cutoff", "AAAGGGC", "AAAGCCCAAAAAACCTTT", energy_threshold=-9.9,
          expected_n=1,
          expected_energies=(-10.0,),
          expected_best=(3, 6, 3, 6, b"GGGC", b"GCCC", b"(((())))")),
    # handle_length_threshold=5 excludes the best structure (its handle is only 4 nt),
    # so FindMinEnLocAlkmer finds nothing at all.
    _case("handle_cutoff", "AAAGGGC", "AAAGCCCAAAAAACCTTT", handle_length_threshold=5,
          expected_n=0,
          expected_energies=()),
    # No complementary region: all A vs all C.
    _case("no_match", "AAAAAAAAAAAA", "CCCCCCCCCCCC",
          expected_n=0,
          expected_energies=()),
    # A perfectly complementary insert in the middle of the sequences.
    _case("perfect_complement", "AAAAAACCCGGGTTTTTT", "TTTTTTGGGCCCAAAAAA",
          expected_n=4,
          expected_energies=(-10.1, -6.6, -4.5, -4.5),
          expected_best=(9, 14, 6, 11, b"GGGTTT", b"GGGCCC", b"(((((())))))")),
    # Antisense (shifted) match: best structure is a long 11 bp stem.
    _case("antisense_match", "AAGGCTTTAATG", "CATTAAGGCTTAAA",
          expected_n=2,
          expected_energies=(-15.1, -1.8),
          expected_best=(1, 11, 0, 10, b"AGGCTTTAATG", b"CATTAAGGCTT", b"((((((((((()))))))))))")),
    # Same sequences with k=5: only the 5-mer kmer-based structure is found.
    _case("k5_short", "AAAGCCCAAAAAACCTTT", "AAAGGGC", k=5,
          expected_n=1,
          expected_energies=(-7.2,),
          expected_best=(13, 17, 0, 4, b"CCTTT", b"AAAGG", b"((((()))))")),
]


@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_e2e(case: E2eCase):
    res = find_min_energy(case.seq, case.seq_compl, k=case.k,
                          energy_threshold=case.energy_threshold,
                          handle_length_threshold=case.handle_length_threshold,
                          need_suboptimal=case.need_suboptimal)

    assert len(res) == case.expected_n
    assert [r[0] for r in res] == pytest.approx(case.expected_energies)
    if case.expected_n:
        start_j, end_j, end_i, start_i, align1, align2, structure = res[0][1:]
        assert (int(start_j), int(end_j), int(end_i), int(start_i),
                bytes(align1), bytes(align2), bytes(structure)) == case.expected_best
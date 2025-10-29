from ..src.fold import Index_seq, FindMinEnLocAlkmer
from pathlib import Path
import numpy as np


def test_e2e():
    k = 3
    GT_threshold = 2
    seq = 'AAAGGGC'
    seq_compl = 'AAAGCCCAAAAAACCTTT'
    energy_threshold = -1
    handle_length_threshold = 0
    need_suboptimal = True

    SCRIPT_PARENT_FOLDER = Path(__file__).resolve().parent
    kmers_stacking_matrix = np.load(SCRIPT_PARENT_FOLDER / ("../data/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy"))
    seq_indxd = Index_seq(seq.encode("ascii"), k)
    seq_compl_indxd = Index_seq(seq_compl.encode("ascii"), k)
    res = FindMinEnLocAlkmer(seq.encode("ascii"), seq_compl.encode("ascii"), 
                            seq_indxd, seq_compl_indxd, k, energy_threshold, handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
    
    assert len(res) == 2
    assert res[0][0] == -10
    assert res[1][0] == -7.2


def test_bench_FindMinEnLocAlkmer(benchmark):
    SCRIPT_PARENT_FOLDER = Path(__file__).resolve().parent
    k = 5
    GT_threshold = 2
    kmers_stacking_matrix = np.load(SCRIPT_PARENT_FOLDER / ("../data/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy"))
    
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
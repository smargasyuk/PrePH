#!/usr/bin/env python

from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort
import numpy as np
import numpy.typing as npt
import typing as tp
from dataclasses import dataclass
from math import ceil
import binascii, itertools, sys, getopt, os
from functools import partial
from sys import getsizeof
inf = float('inf')
# os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pathlib import Path

# Dictionary for nts (used in 1x1, 2x1, 2x2 loops in last 2 dims)
Dic_nt = {b'@': 0, b'A': 1, b'C': 2, b'G': 3, b'T': 4}
# Dictionary of basepairs (used in 1x1, 2x1, 2x2 loops in first 2 dims)
Dic_bp = {b'CG': 0, b'GC': 1, b'GT': 2, b'TG': 3, b'AT': 4, b'TA': 5}

SCRIPT_PARENT_FOLDER = Path(__file__).resolve().parent
stacking_matrix = load(SCRIPT_PARENT_FOLDER / "../lib/stacking_matrix.npy")
bulge_list = load(SCRIPT_PARENT_FOLDER / "../lib/bulge_list.npy")
intl11_matrix = load(SCRIPT_PARENT_FOLDER / "../lib/intl11_matrix.npy")
intl12_matrix = load(SCRIPT_PARENT_FOLDER / "../lib/intl12_matrix.npy")
intl22_matrix = load(SCRIPT_PARENT_FOLDER / "../lib/intl22_matrix.npy")


# Adding for long bulges
TerminalAU = 50


def bsget(bstr: bytes, loc: int) -> bytes:
    return bstr[loc : loc + 1]
    # return bstr[loc]


def Seq_to_bin(seq: bytes) -> int: # Uses bitwise shift to make bin from sequence
    Dict = {b'A': 0b0, b'T': 0b10, b'G': 0b11, b'C': 0b1}
    s = 0
    # print(seq)
    for char in [b'%c' % i for i in seq]:
        # print(char)
        try:
            s = s << 2 | Dict[char]
        except KeyError:
            s = False
    return (s)


def Index_seq(seq: bytes, k: int) -> tp.List[int]: #Uses bitwise shift to divide bin into kmers
    seq_bin = Seq_to_bin(seq)
    if seq_bin is False:
        return (False)
    else:
        seq_indxd_tmp = []
        mask = 2 ** (k * 2) - 1
        for i in range(len(seq) - (k - 1)):
            seq_indxd_tmp.append(mask & seq_bin)
            seq_bin >>= 2
        return (seq_indxd_tmp)

# creates the initial state
# first two return types: npt.NDArray[np.float64], np.ndarray[np.object_], npt.NDArray[np.str_]
def Initiate_with_kmers(seq, seq_compl, seq_indxd_tmp, seq_compl_indxd_tmp, kmers_stacking_matrix, k) -> tp.Tuple[np.ndarray, np.ndarray, float, int, int, int, int, int, int, int, int, bytes, bytes]:
    min_energy = end_pos_i = end_pos_j = start_pos_i = start_pos_j = 0
    #seq_indxd_tmp = Index_seq(seq, k=k)
    #seq_compl_indxd_tmp = Index_seq(seq_compl[::-1], k=k)
    if (seq_indxd_tmp == False) | (seq_compl_indxd_tmp == False):
        return (False, False, False, False, False, False, False, False, False, False, False, False, False)
    
    seq_length = len(seq)
    seq = seq + b'$' * (k + 2)  # (horizontally) _
    seq_compl = b'$' * (k + 2) + seq_compl  # (vertically) |
    seq_compl_length = len(seq_compl)
    D = full((len(seq_compl), len(seq)), inf)  # distance matrix
    # zero_coords = empty((), dtype=object)
    # zero_coords[()] = (0, 0)
    # B = full((len(seq_compl), len(seq)), zero_coords, dtype=object)  # backtracker matrix
    B = np.zeros((len(seq_compl), len(seq), 2), 'i')
    S = empty([len(seq_compl), len(seq)],
                dtype="S" + str(len(seq) + len(seq_compl)))  # dot bracket structure matrix
    for I, kmer_i in enumerate(seq_compl_indxd_tmp):
        i = seq_compl_length - I - k
        for J, kmer_j in enumerate(seq_indxd_tmp):
            j = seq_length - J - 1
            D[i][j] = kmers_stacking_matrix[kmer_j, kmer_i]
            if D[i][j] != inf:
                B[i][j] = (i + k - 1, j - k + 1)
                S[i][j] = b'(' * k + b')' * k
                if D[i][j] < min_energy:
                    min_energy = D[i][j]
                    end_pos_i = i
                    end_pos_j = j
                    start_pos_i = i + k - 1
                    start_pos_j = j - k + 1
    return (
        D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp,
        seq_length, seq_compl_length, seq, seq_compl)


def End_coords(argmin_, i, j, old_end_pos_i, old_end_pos_j, k=3):
    Dict_end_coords = {1: (- 1, + 1),
                       2: (- k, + 1 + k), 3: (- 1 - k, k),
                       4: (- k, + 2 + k), 5: (- 2 - k, k),
                       6: (- 1 - k, + 1 + k),
                       7: (- 1 - k, + 2 + k), 8: (- 2 - k, + 1 + k),
                       9: (- 2 - k, + 2 + k)}
    if argmin_ == 0:
        return (old_end_pos_i, old_end_pos_j)
    else:
        add_i, add_j = Dict_end_coords[argmin_]
        return (i + add_i, j + add_j)


def Start_coords(argmin_, backtrack, old_start_pos_i, old_start_pos_j):
    if argmin_ == 0:
        return (old_start_pos_i, old_start_pos_j)
    else:
        return (backtrack[0], backtrack[1])


def Backtrack(argmin_, old_coords, new_coords):
    if argmin_ == 0:
        return ((0, 0))
    elif argmin_ == 1:
        return (old_coords)
    else:
        return (new_coords)


def Check_ranges_overlap(x1, x2, y1, y2):
    return (((x1 < y2) & (y1 < x2)))


# segment1[start[i,j],end[i,j]], segment2 - part of the square [start[i,j],end[i,j]]
def Check_segments_intersection(segment1, segment2, slope):
    a1 = float(-segment1[0][0] + segment1[1][0]) / (segment1[0][1] - segment1[1][1])
    b1 = -segment1[0][0] - a1 * segment1[0][1]
    if slope == 'v':
        # check vertical
        if segment1[0][1] <= segment2[0][1] & segment2[1][1] <= segment1[1][1]:
            intersection_i = -(a1 * segment2[0][1] + b1)
            if (intersection_i > min(segment1[0][0], segment2[0][0])) | (
                    intersection_i < max(segment1[1][0], segment2[1][0])):
                return False  # intersection is out of bound
            return (int(intersection_i), segment2[0][1])
        else:
            return (False)
    # check horizontal
    elif segment1[0][0] >= segment2[0][0] & segment2[1][0] >= segment1[1][0]:
        intersection_j = (-segment2[0][0] - b1) / a1
        if (intersection_j < max(segment1[0][1], segment2[0][1])) | (
                intersection_j > min(segment1[1][1], segment2[1][1])):
            return False  # intersection is out of bound
        return (segment2[0][0], int(intersection_j))
    else:
        return (False)


def updateStep(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp, seq, seq_compl, kmers_stacking_matrix, k, i, j, I, J):
    # Dictionary for nts (used in 1x1, 2x1, 2x2 loops in last 2 dims)
    Dic_nt = {b'@': 0, b'A': 1, b'C': 2, b'G': 3, b'T': 4}
    # Dictionary of basepairs (used in 1x1, 2x1, 2x2 loops in first 2 dims)
    Dic_bp = {b'CG': 0, b'GC': 1, b'GT': 2, b'TG': 3, b'AT': 4, b'TA': 5}
    cut_idx = S[i][j].find(b')')
    S_head = S[i][j][:cut_idx]
    S_tail = S[i][j][cut_idx:]
    # stem
    new_en = D[i][j] + stacking_matrix[Dic_bp.get(bsget(seq_compl, i - 1) + bsget(seq, j + 1), 6)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 6)]
    
    # print(seq, seq_compl)
    # print("Debug new_en parts: ", D[i][j], stacking_matrix[Dic_bp.get(bsget(seq_compl, i - 1) + bsget(seq, j + 1), 6)], [Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 6)])
    # print("Debug new_en parts: ", Dic_bp,  bsget(seq, j) + bsget(seq_compl, i))

    argmin_ = argmin([0, D[i - 1][j + 1], new_en])
    B[i - 1][j + 1] = Backtrack(argmin_, B[i - 1][j + 1], B[i][j])
    D[i - 1][j + 1] = [0, D[i - 1][j + 1], new_en][argmin_]
    S[i - 1][j + 1] = [b'*', S[i - 1][j + 1], S_head + b'()' + S_tail][argmin_]
    # bulge01 (seq has 1 more nt)
    new_en = D[i][j] + bulge_list[1] + stacking_matrix[Dic_bp.get(bsget(seq_compl, i - 1) + bsget(seq, j + 2), 6)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 6)] + kmers_stacking_matrix[
                    seq_indxd_tmp[J - 1 - k], seq_compl_indxd_tmp[I + k]]
    argmin_ = argmin([0, D[i - k][j + 1 + k], new_en])
    B[i - k][j + 1 + k] = Backtrack(argmin_, B[i - k][j + 1 + k], B[i][j])
    D[i - k][j + 1 + k] = [0, D[i - k][j + 1 + k], new_en][argmin_]
    S[i - k][j + 1 + k] = [b'*', S[i - k][j + 1 + k], S_head + b'.' + b'(' * k + b')' * k + S_tail][argmin_]
    # bulge10 (seq_compl has 1 more nt)
    new_en = D[i][j] + bulge_list[1] + stacking_matrix[Dic_bp.get(bsget(seq_compl, i - 2) + bsget(seq, j + 1), 6)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 6)] + kmers_stacking_matrix[
                    seq_indxd_tmp[J - k], seq_compl_indxd_tmp[I + 1 + k]]
    argmin_ = argmin([0, D[i - 1 - k][j + k], new_en])
    B[i - 1 - k][j + k] = Backtrack(argmin_, B[i - 1 - k][j + k], B[i][j])
    D[i - 1 - k][j + k] = [0, D[i - 1 - k][j + k], new_en][argmin_]
    S[i - 1 - k][j + k] = [b'*', S[i - 1 - k][j + k], S_head + b'(' * k + b')' * k + b'.' + S_tail][argmin_]
    # bulge 02
    new_en = D[i][j] + bulge_list[2] + (
        TerminalAU if Dic_bp.get(bsget(seq_compl, i - 1) + bsget(seq, j + 3), 6) > 1 or Dic_bp.get(
            bsget(seq, j) + bsget(seq_compl, i),
                6) > 1 else 0) + kmers_stacking_matrix[seq_indxd_tmp[J - 2 - k], seq_compl_indxd_tmp[I + k]]
    argmin_ = argmin([0, D[i - k][j + 2 + k], new_en])
    B[i - k][j + 2 + k] = Backtrack(argmin_, B[i - k][j + 2 + k], B[i][j])
    D[i - k][j + 2 + k] = [0, D[i - k][j + 2 + k], new_en][argmin_]
    S[i - k][j + 2 + k] = [b'*', S[i - k][j + 2 + k], S_head + b'..' + b'(' * k + b')' * k + S_tail][
        argmin_]
    # bulge 20
    new_en = D[i][j] + bulge_list[2] + (
        TerminalAU if Dic_bp.get(bsget(seq_compl, i - 3) + bsget(seq, j + 1), 6) > 1 or Dic_bp.get(
            bsget(seq, j) + bsget(seq_compl, i),
                6) > 1 else 0) + kmers_stacking_matrix[seq_indxd_tmp[J - k], seq_compl_indxd_tmp[I + 2 + k]]
    argmin_ = argmin([0, D[i - 2 - k][j + k], new_en])
    B[i - 2 - k][j + k] = Backtrack(argmin_, B[i - 2 - k][j + k], B[i][j])
    D[i - 2 - k][j + k] = [0, D[i - 2 - k][j + k], new_en][argmin_]
    S[i - 2 - k][j + k] = [b'*', S[i - 2 - k][j + k], S_head + b'(' * k + b')' * k + b'..' + S_tail][
        argmin_]
    # loop11
    new_en = D[i][j] + intl11_matrix[Dic_bp.get(bsget(seq_compl, i - 2) + bsget(seq, j + 2), 7)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 7)][Dic_nt.get(bsget(seq_compl, i - 1), 5)][
        Dic_nt.get(bsget(seq, j + 1), 5)] + kmers_stacking_matrix[
                    seq_indxd_tmp[J - 1 - k], seq_compl_indxd_tmp[I + 1 + k]]
    argmin_ = argmin([0, D[i - 1 - k][j + 1 + k], new_en])
    B[i - 1 - k][j + 1 + k] = Backtrack(argmin_, B[i - 1 - k][j + 1 + k], B[i][j])
    D[i - 1 - k][j + 1 + k] = [0, D[i - 1 - k][j + 1 + k], new_en][argmin_]
    S[i - 1 - k][j + 1 + k] = \
        [b'*', S[i - 1 - k][j + 1 + k], S_head + b'.' + b'(' * k + b')' * k + b'.' + S_tail][argmin_]
    # loop12
    new_en = D[i][j] + intl12_matrix[Dic_bp.get(bsget(seq_compl, i - 2) + bsget(seq, j + 3), 7)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 7)][Dic_nt.get(bsget(seq_compl, i - 1), 5)][
        Dic_nt.get(bsget(seq, j + 1), 5)][Dic_nt.get(bsget(seq, j + 2), 5)] + kmers_stacking_matrix[
                    seq_indxd_tmp[J - 2 - k], seq_compl_indxd_tmp[I + 1 + k]]
    argmin_ = argmin([0, D[i - 1 - k][j + 2 + k], new_en])
    B[i - 1 - k][j + 2 + k] = Backtrack(argmin_, B[i - 1 - k][j + 2 + k], B[i][j])
    D[i - 1 - k][j + 2 + k] = [0, D[i - 1 - k][j + 2 + k], new_en][argmin_]
    S[i - 1 - k][j + 2 + k] = \
        [b'*', S[i - 1 - k][j + 2 + k], S_head + b'..' + b'(' * k + b')' * k + b'.' + S_tail][argmin_]
    # loop21
    new_en = D[i][j] + intl12_matrix[Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 7)][
        Dic_bp.get(bsget(seq_compl, i - 3) + bsget(seq, j + 2), 7)][Dic_nt.get(bsget(seq, j + 1), 5)][
        Dic_nt.get(bsget(seq_compl, i - 2), 5)][Dic_nt.get(bsget(seq_compl, i - 1), 5)] + kmers_stacking_matrix[
                    seq_indxd_tmp[J - 1 - k], seq_compl_indxd_tmp[I + 2 + k]]
    argmin_ = argmin([0, D[i - 2 - k][j + 1 + k], new_en])
    B[i - 2 - k][j + 1 + k] = Backtrack(argmin_, B[i - 2 - k][j + 1 + k], B[i][j])
    D[i - 2 - k][j + 1 + k] = [0, D[i - 2 - k][j + 1 + k], new_en][argmin_]
    S[i - 2 - k][j + 1 + k] = \
        [b'*', S[i - 2 - k][j + 1 + k], S_head + b'.' + b'(' * k + b')' * k + b'..' + S_tail][argmin_]
    # loop22
    new_en = D[i][j] + intl22_matrix[Dic_bp.get(bsget(seq_compl, i - 3) + bsget(seq, j + 3), 7)][
        Dic_bp.get(bsget(seq, j) + bsget(seq_compl, i), 7)][Dic_nt.get(bsget(seq_compl, i - 2), 5)][
        Dic_nt.get(bsget(seq_compl, i - 1), 5)][Dic_nt.get(bsget(seq, j + 1), 5)][Dic_nt.get(bsget(seq, j + 2), 5)] + \
                kmers_stacking_matrix[seq_indxd_tmp[J - 2 - k], seq_compl_indxd_tmp[I + 2 + k]]
    argmin_ = argmin([0, D[i - 2 - k][j + 2 + k], new_en])
    B[i - 2 - k][j + 2 + k] = Backtrack(argmin_, B[i - 2 - k][j + 2 + k], B[i][j])
    D[i - 2 - k][j + 2 + k] = [0, D[i - 2 - k][j + 2 + k], new_en][argmin_]
    S[i - 2 - k][j + 2 + k] = \
        [b'*', S[i - 2 - k][j + 2 + k], S_head + b'..' + b'(' * k + b')' * k + b'..' + S_tail][argmin_]
    # check if found min Energy
    argmin_ = argmin([min_energy, D[i - 1][j + 1],
                        D[i - k][j + 1 + k], D[i - 1 - k][j + k],
                        D[i - k][j + 2 + k], D[i - 2 - k][j + k],
                        D[i - 1 - k][j + 1 + k],
                        D[i - 1 - k][j + 2 + k], D[i - 2 - k][j + 1 + k],
                        D[i - 2 - k][j + 2 + k]])
    end_pos_i, end_pos_j = End_coords(argmin_, i, j, end_pos_i, end_pos_j, k=k)
    start_pos_i, start_pos_j = Start_coords(argmin_, B[i][j], start_pos_i, start_pos_j)
    min_energy = [min_energy, D[i - 1][j + 1],
                    D[i - k][j + 1 + k], D[i - 1 - k][j + k],
                    D[i - k][j + 2 + k], D[i - 2 - k][j + k],
                    D[i - 1 - k][j + 1 + k],
                    D[i - 1 - k][j + 2 + k], D[i - 2 - k][j + 1 + k],
                    D[i - 2 - k][j + 2 + k]][argmin_]
    
    return D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j


def updateLoop(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp, seq_length, seq_compl_length, seq, seq_compl, kmers_stacking_matrix, k):
    for i in range(len(seq_compl) - k, k + 2, -1):
        I = seq_compl_length - k - i
        for j in range(k - 1, len(seq) - k - 3):
            J = seq_length - j - 1
            if (D[i][j] != inf and D[i][j] != 0):  # found kmer stacking
                D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j = updateStep(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp, seq, seq_compl, kmers_stacking_matrix, k, i, j, I, J)
    return D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j


def getAlignments(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq, seq_compl, k, need_suboptimal: bool, energy_threshold: float, handle_length_threshold: float):
        alignments = []
        # save optimal structure
        alignments.append((min_energy / 100, start_pos_j, end_pos_j, end_pos_i - k - 2, start_pos_i - k - 2,
                            seq[start_pos_j:end_pos_j + 1], seq_compl[end_pos_i:start_pos_i + 1],
                            S[end_pos_i, end_pos_j]))
        if need_suboptimal:
            min_energy_potential = 0
            while True:
                # LOOK FOR SUBOPTIMAL STRUCTURES
                ## put zeros in square from start to end of the best alignment
                len_zero_i = start_pos_i - end_pos_i + 1
                len_zero_j = end_pos_j - start_pos_j + 1
                zero_matrix_D = full((len_zero_i, len_zero_j), 0)
                D[end_pos_i:start_pos_i + 1, start_pos_j: end_pos_j + 1] = zero_matrix_D
                ## find next suboptimal structure
                argmin_ = D.argmin()
                end_pos_i_potential, end_pos_j_potential = unravel_index(argmin_, D.shape)
                start_pos_i_potential, start_pos_j_potential = Start_coords(1, B[end_pos_i_potential][
                    end_pos_j_potential], 0, 0)
                min_energy_potential = D[end_pos_i_potential][end_pos_j_potential]
                ## check it's energy is higher than threshold
                if min_energy_potential / 100.0 > energy_threshold:
                    break
                    ## check it doesn't overlap any previous alignments
                overlap = False
                for alignment in alignments:
                    start_pos_j = alignment[1]
                    end_pos_j = alignment[2] + k + 2
                    end_pos_i = alignment[3]
                    start_pos_i = alignment[4] + k + 2
                    if (Check_ranges_overlap(start_pos_j, end_pos_j, start_pos_j_potential,
                                                end_pos_j_potential) & Check_ranges_overlap(end_pos_i_potential,
                                                                                            start_pos_i_potential,
                                                                                            end_pos_i, start_pos_i)):
                        overlap = True
                        break
                        ## if it overlaps, delete it from D matrix
                if overlap:
                    D[end_pos_i_potential][end_pos_j_potential] = 0
                ## if it doesn't overlap, save the suboptimal structure
                else:
                    min_energy = min_energy_potential
                    start_pos_i = start_pos_i_potential
                    end_pos_i = end_pos_i_potential
                    start_pos_j = start_pos_j_potential
                    end_pos_j = end_pos_j_potential
                    if (end_pos_j - start_pos_j + 1 >= handle_length_threshold) & (
                            start_pos_i - end_pos_i + 1 >= handle_length_threshold):  # check the al is long enough
                        alignments.append(
                            (min_energy / 100, start_pos_j, end_pos_j, end_pos_i - k - 2, start_pos_i - k - 2,
                                seq[start_pos_j:end_pos_j + 1], seq_compl[end_pos_i:start_pos_i + 1],
                                S[end_pos_i, end_pos_j]))
        return (alignments)


def FindMinEnLocAlkmer(seq, seq_compl, seq_indxd, seq_compl_indxd, k, energy_threshold, handle_length_threshold,
                       need_suboptimal, kmers_stacking_matrix):
    
    # print(seq, seq_compl, energy_threshold, handle_length_threshold, need_suboptimal)
    
    if (seq_indxd == False) | (seq_compl_indxd == False):
        return 0
    
    seq_indxd_tmp = seq_indxd[:]
    seq_compl_indxd_tmp = seq_compl_indxd[:]
    D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp, \
    seq_length, seq_compl_length, seq, seq_compl = \
        Initiate_with_kmers(seq, seq_compl, seq_indxd_tmp, seq_compl_indxd_tmp, kmers_stacking_matrix, k)

    # exit if input is empty or 
    if (min_energy == 0) | (seq == False):
        return 0

    seq_indxd_tmp.extend([kmers_stacking_matrix.shape[0] - 1] * (k + 2))
    seq_compl_indxd_tmp.extend([kmers_stacking_matrix.shape[0] - 1] * (k + 2))

    # main loop: go through matrices and fill them in
    D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j = updateLoop(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq_indxd_tmp, seq_compl_indxd_tmp, seq_length, seq_compl_length, seq, seq_compl, kmers_stacking_matrix, k)

    
    # exit if the optimal alignment is short / low energy
    if (min_energy / 100.0 > energy_threshold) | (end_pos_j - start_pos_j + 1 < handle_length_threshold) | (
            start_pos_i - end_pos_i + 1 < handle_length_threshold):
        return 0
    
    # get alignments from state
    alignments = getAlignments(D, B, S, min_energy, end_pos_i, end_pos_j, start_pos_i, start_pos_j, seq, seq_compl, k, need_suboptimal, energy_threshold, handle_length_threshold)

    return (alignments)

def main(argv):
    seq = b''
    seq_compl = b''
    k = 5
    handle_length_threshold = 10
    energy_threshold = -15
    need_suboptimal = True
    GT_threshold = 2
    try:
        opts, args = getopt.getopt(argv, "h:f:s:k:a:e:u:d:",
                                   ["help=", "first_seq=", "second_seq=", "k=", "handle_len_min", "energy_max",
                                    "need_subopt", "gt_threshold"])
    except getopt.GetoptError:
        print(
            b'fold.py -f <first_seq> -s <second_seq> -k <kmer_length> -a <handle_len_min> -e <energy_max, kcal/mol> -u <need_subopt> -d <gt_threshold>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', b'--help'):
            print(
                b'fold.py -f <first_seq> -s <second_seq> -k <kmer_length> -a <handle_len_min> -e <energy_max, kcal/mol> -u <need_subopt> -d <gt_threshold>')
            sys.exit()
        elif opt in ("-f", "--first_seq"):
            seq = arg
        elif opt in ("-s", "--second_seq"):
            seq_compl = arg
        elif opt in ("-k", "--k"):
            k = int(arg)
        elif opt in ("-a", "--handle_len_min"):
            handle_length_threshold = int(arg)
        elif opt in ("-e", "--energy_max"):
            energy_threshold = float(arg)
        elif opt in ("-u", "--need_subopt"):
            need_suboptimal = arg
        elif opt in ("-d", "--gt_threshold"):
            GT_threshold = int(arg)

    if need_suboptimal == 'False':
        need_suboptimal = False
    elif need_suboptimal == 'True':
        need_suboptimal = True

    kmers_stacking_matrix = load(SCRIPT_PARENT_FOLDER / ("../data/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy"))
    seq_indxd = Index_seq(seq.encode("ascii"), k)
    seq_compl_indxd = Index_seq(seq_compl.encode("ascii"), k)
    res = FindMinEnLocAlkmer(seq.encode("ascii"), seq_compl.encode("ascii"), seq_indxd, seq_compl_indxd, k, energy_threshold, handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
    print(res)
    return(res)

if __name__ == '__main__':
    main(sys.argv[1:])

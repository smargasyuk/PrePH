#!/usr/bin/env python2
import itertools
import pandas as pd
import numpy as np
from numpy import full
import sys, getopt, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

inf = float('inf')

Dic_bp = {'CG': 0, 'GC': 1, 'GT': 2, 'TG': 3, 'AT': 4, 'TA': 5}
stacking_matrix = np.load("../lib/stacking_matrix.npy")
Dict_nts = {'A': 0b0, 'T': 0b10, 'G': 0b11, 'C': 0b1}
List_pairs = ['AT', 'TA', 'GC', 'CG', 'TG', 'GT']
bases = ['A', 'T', 'G', 'C']


def CalculateStackingEnergy(seq, seq_compl):
    energy = 0
    i = 0
    j = 0
    while i < len(seq) - 1 and j < len(seq) - 1:
        energy_add = stacking_matrix[Dic_bp.get(seq_compl[i + 1] + seq[j + 1], 6)][Dic_bp.get(seq[j] + seq_compl[i], 6)]
        i += 1
        j += 1
        energy += energy_add
    return (energy)


def Seq_to_bin(seq):
    seq_bin = 1
    for char in seq:
        seq_bin = seq_bin << 2 | Dict_nts[char]
    return (seq_bin)


def Precalculatekmers(k, GT_threshold, to_remove):
    try: 
        os.makedirs('../data/')
    except OSError:
        if not os.path.isdir('../data/'):
            raise
    kmers_for_bin = [''.join(p) for p in itertools.product(bases, repeat=k)]
    with open('../data/kmers_list_' + str(k) + 'mers_no$.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % kmer for kmer in kmers_for_bin)
    kmers_bin = []
    for kmer in kmers_for_bin:
        kmer_bin = Seq_to_bin(kmer)
        kmers_bin.append(kmer_bin)
    kmers_bin[:] = [x - len(bases) ** k for x in kmers_bin]
    Dict_kmers = dict(zip(kmers_bin, kmers_for_bin))
    kmers_array = np.full((len(bases) ** k, len(bases) ** k), inf)
    for i in range(len(bases) ** k):
        for j in range(len(bases) ** k):
            #pairs = zip(Dict_kmers[i], Dict_kmers[j])
            pairs = zip(Dict_kmers[i], Dict_kmers[j][::-1])
            GT_count = sum(1 for pair in pairs if (pair == ('G', 'T')) | (pair == ('T', 'G')))
            stem_count = sum(1 for pair in pairs if pair[0] + pair[1] in List_pairs)
            if GT_count <= GT_threshold and stem_count == k and not (Dict_kmers[i] in to_remove) and not (Dict_kmers[j][::-1] in to_remove):
                energy = CalculateStackingEnergy(Dict_kmers[i], Dict_kmers[j][::-1])
                kmers_array[i][j] = energy
    np.save("../data/" + str(k) + str(GT_threshold) + "mers_stacking_energy_no$", kmers_array)
    kmers_array_for_dollar = full((kmers_array.shape[0] + 1, kmers_array.shape[0] + 1), inf)
    kmers_array_for_dollar[:-1, :-1] = kmers_array
    np.save("../data/" + str(k) + str(GT_threshold) +  "mers_stacking_energy_binary", kmers_array_for_dollar)
    return (0)

def main(argv):
    k = 5
    GT_threshold = 2
    to_remove = ''
    try:
        opts, args = getopt.getopt(argv, "h:k:g:r:",
                                   ["help=", "k=", "gt_threshold"])
    except getopt.GetoptError:
        print('PrecalculateStrackingEnergies.py -k <kmer_length> -g <gt_amount_in_kmer_max> -r <kmers_to_remove>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('PrecalculateStrackingEnergies.py -k <kmer_length> -g <gt_amount_in_kmer_max> -r <kmers_to_r_move>')
            sys.exit()
        elif opt in ("-k", "--k"):
            k = int(arg)
        elif opt in ("-g", "--gt_threshold"):
            GT_threshold = int(arg)
        elif opt in ("-r", "--kmers_to_r_move"):
            to_remove = arg
    if to_remove != '':
        text_file = open(to_remove, "r")
        to_remove = text_file.read().split('\n')
        print('I would remove ' + str(len(to_remove)) + ' kmers')
    Precalculatekmers(k, GT_threshold, to_remove)


if __name__ == '__main__':
    main(sys.argv[1:])

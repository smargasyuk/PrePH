#!/usr/bin/env python3
import itertools
import numpy as np
import sys
import os
import click

from .asset_loader import load_static_data, get_kmer_table_path

os.chdir(os.path.dirname(os.path.abspath(__file__)))

inf = float('inf')

Dic_bp = {'CG': 0, 'GC': 1, 'GT': 2, 'TG': 3, 'AT': 4, 'TA': 5}
stacking_matrix = load_static_data()[0]
Dict_nts = {'A': 0b0, 'T': 0b10, 'G': 0b11, 'C': 0b1}
List_pairs = ['AT', 'TA', 'GC', 'CG', 'TG', 'GT']
bases = ['A', 'T', 'G', 'C']


def CalculateStackingEnergy(seq, seq_compl):
    energy = 0
    i = 0
    j = 0
    while i < len(seq) - 1 and j < len(seq) - 1:
        energy_add = stacking_matrix[
            Dic_bp.get(seq_compl[i + 1] + seq[j + 1], 6)
        ][
            Dic_bp.get(seq[j] + seq_compl[i], 6)
        ]
        i += 1
        j += 1
        energy += energy_add
    return energy


def Seq_to_bin(seq):
    seq_bin = 1
    for char in seq:
        seq_bin = seq_bin << 2 | Dict_nts[char]
    return seq_bin


def Precalculatekmers(k, GT_threshold, to_remove):
    target_dir = get_kmer_table_path(k, GT_threshold).parent
    os.makedirs(target_dir, exist_ok=True)

    kmers_for_bin = [''.join(p) for p in itertools.product(bases, repeat=k)]
    with open(target_dir / ('kmers_list_' + str(k) + 'mers_no$.txt'), 'w') as filehandle:
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
            pairs = list(zip(Dict_kmers[i], Dict_kmers[j][::-1]))

            GT_count = sum(
                1 for pair in pairs
                if pair == ('G', 'T') or pair == ('T', 'G')
            )
            stem_count = sum(
                1 for pair in pairs
                if pair[0] + pair[1] in List_pairs
            )

            if (
                GT_count <= GT_threshold
                and stem_count == k
                and Dict_kmers[i] not in to_remove
                and Dict_kmers[j][::-1] not in to_remove
            ):
                energy = CalculateStackingEnergy(Dict_kmers[i], Dict_kmers[j][::-1])
                kmers_array[i][j] = energy

    np.save(get_kmer_table_path(k, GT_threshold, suffix="mers_stacking_energy_no$.npy"), kmers_array)

    kmers_array_for_dollar = np.full(
        (kmers_array.shape[0] + 1, kmers_array.shape[0] + 1), inf
    )
    kmers_array_for_dollar[:-1, :-1] = kmers_array
    np.save(
        get_kmer_table_path(k, GT_threshold),
        kmers_array_for_dollar,
    )
    return 0

@click.command()
@click.option("-k", type=int, default=5)
@click.option("-g", '--gt-threshold', type=int, default=2)
@click.option("-r", '--kmers-to-remove', type=str, default='')
def main(k, gt_threshold, kmers_to_remove):
    if kmers_to_remove != '':
        with open(kmers_to_remove, "r") as text_file:
            to_remove = text_file.read().split('\n')
        print('I would remove ' + str(len(to_remove)) + ' kmers')
    else:
        to_remove = []

    Precalculatekmers(k, gt_threshold, to_remove)


if __name__ == '__main__':
    main(sys.argv[1:])
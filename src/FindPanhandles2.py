from typing import Any
from .fold import FindMinEnLocAlkmer, Index_seq
import polars as pl
import click
import numpy as np
import itertools
from pathlib import Path
import enum
from multiprocessing import Pool
from dataclasses import dataclass
from collections.abc import Iterator
from tqdm import tqdm


class PairingMode(enum.Enum):
    any = enum.auto()
    cis = enum.auto()
    trans = enum.auto()


def format_preph_results(entries_in):
    results_table = pl.DataFrame([{'row': result['r1']['group_name'], 'energy': alignment[0], 'interval1': result['r1']['interval_id'], 'interval2': result['r2']['interval_id'],
                    'start_al1': alignment[1], 'end_al1': alignment[2],
                    'start_al2': alignment[3], 'end_al2': alignment[4],
                    'alignment1': alignment[5], 'alignment2': alignment[6], 'structure': alignment[7]} for result in entries_in for alignment in result['alignment']])\
            .with_columns(pl.col(x).map_elements(lambda y: y.decode('ascii'), return_dtype=str) for x in ['alignment1', 'alignment2', 'structure'])\
            .filter((pl.col("start_al1") < pl.col("start_al2")) | (pl.col('interval1') != pl.col('interval2')))
    return results_table


def attach_index_columns(df, k):
    df2 = df\
        .with_columns((pl.col('chrom') + "_" + pl.col("chromStart").cast(str) + "_" + pl.col("chromEnd").cast(str)+ "_" + pl.col("strand")).alias("interval_id"))\
        .with_columns(pl.col('sequence').map_elements(lambda x: x.encode('ascii'), return_dtype=pl.Binary).alias('sequence_bs'))\
        .with_columns(pl.col('sequence_bs').map_elements(lambda x: Index_seq(x, k), return_dtype=pl.List(pl.Int64)).alias('sequence_idx'))
    return df2


def generate_pairs(df, pairing: PairingMode, match_column:str = "name") -> Iterator[tuple[dict[str, Any], dict[str, Any]]]:
    iter0 = df.sort(by='chromStart').iter_rows(named=True)
    iter1 = itertools.combinations_with_replacement(iter0, 2)
    if pairing == PairingMode.any:
        return iter1
    if pairing == PairingMode.cis:
        return filter(lambda pair: pair[0][match_column] == pair[1][match_column], iter1)
    if pairing == PairingMode.trans:
        return filter(lambda pair: pair[0][match_column] != pair[1][match_column], iter1)


@dataclass
class PrephParameters:
    k: int
    gt_threshold: int
    energy_max: float
    panhandle_len_min: int
    need_subopt: bool
    kmers_stacking_matrix: np.ndarray


def pack_preph_params(k, gt_threshold, energy_max, panhandle_len_min, need_subopt, kmers_stacking_matrix):
    return PrephParameters(k, gt_threshold, energy_max, panhandle_len_min, need_subopt, kmers_stacking_matrix)


def apply_preph_to_row_pair(args) -> dict[str, list[Any] | int]:
    r1, r2, preph_params = args
    pp = preph_params
    # print(r1['name'], r2['name'], flush=True)
    pw_results = FindMinEnLocAlkmer(r1['sequence_bs'], r2['sequence_bs'], r1['sequence_idx'], r2['sequence_idx'], pp.k, pp.energy_max, pp.panhandle_len_min, pp.need_subopt, pp.kmers_stacking_matrix)
    return {'r1': r1, 'r2': r2, 'alignment': pw_results}


@click.command()
@click.option("--input", required=True)
@click.option("--output", required=True)
@click.option("-k", type=int, default=5)
@click.option("-d", '--gt-threshold', type=int, default=2)
@click.option("-e", '--energy-max', type=float, default=-15)
@click.option("-a", '--panhandle-len-min', type=int, default=10)
@click.option('--need-subopt/--no-need-subopt', default=True)
@click.option('--pairing', type=click.Choice(PairingMode), default=PairingMode.any)
@click.option('-j', '--jobs', type=int, default=8)
def main(input, output, k, gt_threshold, energy_max, panhandle_len_min, need_subopt, pairing, jobs):
    df1 = pl.read_csv(input, separator='\t')
    df2 = attach_index_columns(df1, k)
    
    SCRIPT_PARENT_FOLDER = Path(__file__).resolve().parent
    kmers_stacking_matrix = np.load(SCRIPT_PARENT_FOLDER /  ("../data/" + str(k) + str(gt_threshold) + "mers_stacking_energy_binary.npy"))
    preph_params = pack_preph_params(k, gt_threshold, energy_max, panhandle_len_min, need_subopt, kmers_stacking_matrix)

    # iterate over tuples (group_id, grouped_df)
    group_it = df2.group_by('group_name')
    # for each group, generate an iterator over pairs of rows, then chain them together
    pair_it0 = itertools.chain.from_iterable(generate_pairs(v, pairing) for group_id, v in group_it)
    pair_it1 = ((r1, r2, preph_params) for r1, r2 in pair_it0)

    with Pool(jobs) as p:
        results0 = list(tqdm(p.imap_unordered(apply_preph_to_row_pair, pair_it1)))

    results: list[dict[str, list[Any]]] = [pwr for pwr in results0 if pwr['alignment'] != 0]  # pyright: ignore[reportAssignmentType]

    results_table = format_preph_results(results)
    results_table.write_csv(output, separator='\t')

if __name__ == '__main__':
    main()

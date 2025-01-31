from .fold import FindMinEnLocAlkmer, Index_seq
import polars as pl
import click
import numpy as np
from itertools import combinations_with_replacement
from pathlib import Path


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


@click.command()
@click.option("--input", required=True)
@click.option("-k", type=int, default=5)
@click.option("-d", '--gt-threshold', type=int, default=2)
@click.option("-e", '--energy-max', type=float, default=-15)
@click.option("-a", '--panhandle-len-min', type=int, default=10)
@click.option('--need-subopt/--no-need-subopt', default=True)
@click.option("--output", required=True)
def main(input, output, k, gt_threshold, energy_max, panhandle_len_min, need_subopt):
    df1 = pl.read_csv(input, separator='\t')
    df2 = attach_index_columns(df1, k)
    
    SCRIPT_PARENT_FOLDER = Path(__file__).resolve().parent
    kmers_stacking_matrix = np.load(SCRIPT_PARENT_FOLDER /  ("../data/" + str(k) + str(gt_threshold) + "mers_stacking_energy_binary.npy"))

    results = []
    for group_id, v in df2.group_by('group_name'):
        print(group_id[0])
        v = v.sort('chromStart')
        for r1, r2 in combinations_with_replacement(v.iter_rows(named=True), 2):
            print(r1['name'], r2['name'])
            pw_results = FindMinEnLocAlkmer(r1['sequence_bs'], r2['sequence_bs'], r1['sequence_idx'], r2['sequence_idx'], k, energy_max, panhandle_len_min, need_subopt, kmers_stacking_matrix)
            if pw_results != 0:
                results.append({'r1': r1, 'r2': r2, 'alignment': pw_results})

    results_table = format_preph_results(results)
    results_table.write_csv(output, separator='\t')

if __name__ == '__main__':
    main()

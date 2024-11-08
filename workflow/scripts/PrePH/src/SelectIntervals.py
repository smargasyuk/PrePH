#!/usr/bin/env python2
import pandas as pd
import sys, getopt
from subprocess import call
from functools import partial
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import argparse

def select_genes(path_to_anno, gene_type, output_dir):
    headers = 0
    with open(path_to_anno) as f:
        while 1:
            line = f.readline()
            if line.startswith("#"):
                headers += 1
            else:
                break
    df = pd.read_csv(path_to_anno, sep="\t", skiprows=headers, header=None, names=['V1','V2','V3','V4','V5','V6','V7','V8', 'V9'])
    df_genes = df.loc[df.V3 == 'gene']
    if gene_type == 'coding':
        print('Selected coding genes')
        df_genes = df_genes.loc[df_genes.V9.str.contains('protein_coding')]
    elif gene_type == 'noncoding':
        print('Selected non-coding genes')
        df_genes = df_genes.loc[~df_genes.V9.str.contains('protein_coding')]
    else:
        print('Selected all genes')
    gene_id_place = [i for i, s in enumerate(df_genes.iloc[0, 8].split(";")) if 'gene_id' in s][0]
    df_genes['gene_id'] = df_genes.V9.str.split("; ", n=2, expand=True)[0]
    df_genes['gene_id'] = df_genes["gene_id"].str.split(' "', n=3, expand=True)[1].str[:-1]
    gene_name_place = [i for i, s in enumerate(df_genes.iloc[0, 8].split(";")) if 'gene_name' in s][0]
    df_genes['gene_name'] = df_genes.V9.str.split("; ", n=5, expand=True)[gene_name_place]
    df_genes['gene_name'] = df_genes["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    df_genes['gene_id_name'] = df_genes['gene_id'] + '_' + df_genes['gene_name']
    df_coding_bed = df_genes[['V1', 'V4', 'V5', 'gene_id_name', 'V6', 'V7']]
    df_coding_bed.to_csv(output_dir + "/genes.bed", sep="\t", index=False, header=False)
    return(0)



def select_intronic_regions(path_to_anno, flank_length, gene_type, output_dir):
    headers = 0
    with open(path_to_anno) as f:
        while 1:
            line = f.readline()
            if line.startswith("#"):
                headers += 1
            else:
                break
    df = pd.read_csv(path_to_anno, sep="\t", skiprows=headers, header=None)
    # select CDSs
    if gene_type == 'coding':
        CDS = df[df[2] == 'CDS']
    else:
        CDS = df[df[2] == 'exon']
    CDS['gene_id'] = CDS[8].str.split("; ", n=2, expand=True)[0]
    CDS.gene_id = CDS.gene_id.replace("gene_id ", "", regex=True)
    CDS.gene_id = CDS.gene_id.replace("\"", "", regex=True)  
    gene_name_place = [i for i, s in enumerate(CDS.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    CDS['gene_name'] = CDS[8].str.split("; ", n=5, expand=True)[gene_name_place]
    CDS.gene_name = CDS["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    CDS['gene_id_name'] = CDS['gene_id'] + '_' + CDS['gene_name']
    CDS = CDS[[0, 3, 4, 'gene_id_name', 5, 6]]
    CDS.to_csv("{0}/CDS.bed".format(output_dir), sep="\t", index=False, header=False)
    # subtract CDS from genes
    x = 'bedtools subtract -a {0}/genes.bed -b {0}/CDS.bed -s > {0}/genes_no_CDS.bed'.format(output_dir)
    call([x], shell = True)
    dt = pd.read_csv('{0}/genes_no_CDS.bed'.format(output_dir),  sep="\t", header=None)
    dt = dt.sort_values(by =[3, 1, 2])
    # add flanks
    dt[1] = dt[1] - flank_length
    dt[2] = dt[2] + flank_length
    dt[4] = 1
    dt.to_csv('{0}/introns_with_flanks_python.bed'.format(output_dir), sep="\t", index=False, header=False)
    call(['rm {0}/genes_no_CDS.bed'.format(output_dir)], shell = True)
    call(['rm {0}/CDS.bed'.format(output_dir)], shell = True)
    return(0)

def filter(filter_table, length_min, output_dir):
    ## sort
    call(['bedtools sort -i {0}/conin_python_long.tsv > {0}/conin_python_long_sorted.tsv'.format(output_dir)], shell = True)
    ## sort repeats and subtract and merge
    x = 'bedtools sort -i ' + filter_table + \
    '|bedtools subtract -a {0}/conin_python_long_sorted.tsv -b stdin'.format(output_dir) + \
    '|bedtools sort -i stdin | bedtools merge -s -c 4 -o distinct -i stdin > {0}/conin_python_long_sorted_filtered.tsv'.format(output_dir)
    call([x], shell = True)
    ## select only long enough
    df = pd.read_csv('{0}/conin_python_long_sorted_filtered.tsv'.format(output_dir), sep="\t", header=None)
    df.columns = ['chrom', 'chromStart', 'chromEnd', 'gene']
    df = df[df.chromEnd - df.chromStart + 1 >= length_min]
    df.to_csv('{0}/conin_python_long_sorted_filtered.tsv'.format(output_dir), sep="\t", index=False, header=False)
    print('Filtered')

def make_output_table(output_dir):
    # separate intersected genes
    df = pd.read_csv('{0}/conin_python_long_sorted_filtered.tsv'.format(output_dir), sep="\t", header=None)
    df.columns = ['chrom', 'chromStart', 'chromEnd', 'gene']
    df1 = df[df['gene'].str.contains(',')]
    df2 = df[~df['gene'].str.contains(',')]
    df1_gene_names = df1["gene"].str.split(',', n=2, expand=True)[0].append(df1["gene"].str.split(',', n=2, expand=True)[1])
    df1_coords = df1[['chrom', 'chromStart', 'chromEnd']].append(df1[['chrom', 'chromStart', 'chromEnd']])
    df1_new = pd.concat([df1_coords, df1_gene_names], axis = 1)
    df1_new.columns = ['chrom', 'chromStart', 'chromEnd', 'gene']
    df = df1_new.append(df2)

    genes = pd.read_csv('{0}/genes.bed'.format(output_dir), sep="\t", header=None)
    genes.columns = ['chrom', 'start_gene', 'end_gene', 'gene', 'score', 'strand']
    # add genes coords
    df = df.merge(genes[['start_gene', 'end_gene', 'gene', 'strand']], on = 'gene', how = 'inner')
    # header
    df.drop_duplicates(subset = ['chrom', 'chromStart', 'chromEnd', 'strand', 'start_gene', 'end_gene'], keep = 'first', inplace = True) 
    df['name'] = range(1,df.shape[0]+1)
    df['score'] = 1
    df = df.sort_values(by =['chrom', 'chromStart', 'chromEnd'])
    df = df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'start_gene', 'end_gene']]
    df.to_csv('{0}/conin_python_long_filtered_final.tsv'.format(output_dir), sep="\t", index=False)
    print('selected intervals')


def main(argv):

    parser = argparse.ArgumentParser(description='Extract conserved intervals')
    parser.add_argument('-a', '--annotation', dest="path_to_anno", required=True)
    parser.add_argument('-c', '--cons_regions', dest="path_to_cons", required=True)
    parser.add_argument('-l', '--handle_len_min', dest="length_min", default="10")
    parser.add_argument('-f', '--flanks', dest="flank_length", default=10, type=int)
    parser.add_argument('-t', '--coding', dest="gene_type", default="coding")
    parser.add_argument('-r', '--repeats_table', dest="filter_table")
    parser.add_argument('-o', '--output-dir', default='../data')
    args = parser.parse_args()


    select_genes(args.path_to_anno, args.gene_type, args.output_dir)
    select_intronic_regions(args.path_to_anno, args.flank_length, args.gene_type, args.output_dir)
    call(['./Select_conins_new.sh', args.path_to_cons, str(args.length_min), args.output_dir])
    if args.filter_table:
        filter(args.filter_table, int(args.length_min), args.output_dir)
    make_output_table(args.output_dir)


if __name__ == '__main__':
    main(sys.argv[1:])

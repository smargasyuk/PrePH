#!/usr/bin/env python2
import pandas as pd
import sys, getopt
from subprocess import call
from functools import partial
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def select_coding_genes(path_to_anno):
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
    df_genes = df_genes.loc[df_genes.V9.str.contains('protein_coding')]

    
    gene_id_place = [i for i, s in enumerate(df_genes.iloc[1, 8].split(";")) if 'gene_id' in s][0]
    df_genes['gene_id'] = df_genes.V9.str.split("; ", n=2, expand=True)[0]
    df_genes['gene_id'] = df_genes["gene_id"].str.split(' "', n=3, expand=True)[1].str[:-1]

    gene_name_place = [i for i, s in enumerate(df_genes.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    df_genes['gene_name'] = df_genes.V9.str.split("; ", n=5, expand=True)[gene_name_place]
    df_genes['gene_name'] = df_genes["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    
    df_genes['gene_id_name'] = df_genes['gene_id'] + '_' + df_genes['gene_name']
    df_coding_bed = df_genes[['V1', 'V4', 'V5', 'gene_id_name', 'V6', 'V7']]
    df_coding_bed.to_csv("../data/coding_genes.bed", sep="\t", index=False, header=False)
    return(0)


def select_intronic_regions(path_to_anno, flank_length):
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
    CDS = df[df[2] == 'CDS']
    CDS['gene_id'] = CDS[8].str.split("; ", n=2, expand=True)[0]
    CDS.gene_id = CDS.gene_id.replace("gene_id ", "", regex=True)
    CDS.gene_id = CDS.gene_id.replace("\"", "", regex=True)  
    gene_name_place = [i for i, s in enumerate(CDS.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    CDS['gene_name'] = CDS[8].str.split("; ", n=5, expand=True)[gene_name_place]
    CDS.gene_name = CDS["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    CDS['gene_id_name'] = CDS['gene_id'] + '_' + CDS['gene_name']
    CDS = CDS[[0, 3, 4, 'gene_id_name', 5, 6]]
    CDS.to_csv("../data/CDS.bed", sep="\t", index=False, header=False)

    # subtract CDS from genes
    x = 'bedtools subtract -a ../data/coding_genes.bed -b ../data/CDS.bed -s > ../data/coding_genes_no_CDS.bed' 
    call([x], shell = True)
    dt = pd.read_csv('../data/coding_genes_no_CDS.bed',  sep="\t", header=None)
    dt = dt.sort_values(by =[3, 1, 2])

    # add flanks
    dt[1] = dt[1] - flank_length
    dt[2] = dt[2] + flank_length
    dt[4] = 1
    dt.to_csv('../data/introns_with_flanks_python.bed', sep="\t", index=False, header=False)
    call(['rm ../data/coding_genes_no_CDS.bed'], shell = True)
    call(['rm ../data/CDS.bed'], shell = True)
    return(0)


def main(argv):
    path_to_anno="../data/dm6/Drosophila_melanogaster.BDGP6.22.96_chr.gtf"
    path_to_cons="../data/dm6/phastConsElements124way.txt"
    length_min=10
    flank_length=10

    try:
        opts, args = getopt.getopt(argv, "h:a:c:l:i:r:f:",
                                   ["help=", "annotation=", "cons_regions=", "handle_len_min=", "always_intronic=", "chromosomes=", "flanks="])
    except getopt.GetoptError:
        print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks>')
            sys.exit()
        elif opt in ("-a", "--annotation"):
            path_to_anno = arg
        elif opt in ("-c", "--cons_regions"):
            path_to_cons = arg
        elif opt in ("-l", "--handle_len_min"):
            length_min = arg
        elif opt in ("-f", "--flanks"):
            flank_length = int(arg)

    select_coding_genes(path_to_anno)
    print('selected coding genes')
    
    select_intronic_regions(path_to_anno, flank_length)
    print('selected introns with flanks')
    
    call(['./Select_conins_new.sh', path_to_cons, length_min])

    df = pd.read_csv('../data/conin_python_long_coding_final.tsv', sep="\t")
    df['name'] = range(1,df.shape[0]+1)
    df['score'] = 1
    df = df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'start_gene', 'end_gene']]
    df.drop_duplicates(subset = ['chrom', 'chromStart', 'chromEnd', 'strand', 'start_gene', 'end_gene'], keep = 'first', inplace = True) 
    df.to_csv('../data/conin_python_long_coding_final.tsv', sep="\t", index=False)
    print('selected intervals')


if __name__ == '__main__':
    main(sys.argv[1:])

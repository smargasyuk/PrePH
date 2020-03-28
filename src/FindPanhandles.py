#!/usr/bin/env python2
from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort
from math import ceil
import re, sys, getopt, itertools, binascii, time, subprocess
from functools import partial
#sys.path.append('../../../tools/')
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
from Bio.Seq import Seq
sys.path.insert(0, './')
from fold import * 


inf = float('inf')


def GetSequencesForDF(genome, row):
    return (str(genome[row['chrom']][row['chromStart'] - 1:row['chromEnd']]))

def MakeComplement(row):
    if row['strand'] == '-':
        return(str(Seq(row['sequences']).complement()))
    else:
        return(row['sequences'])

def Find_panhandles_one_gene(lock, df, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                             need_suboptimal, kmers_stacking_matrix, path_to_ph, gene):
    with open('../data/genes_done.txt', 'a') as done_f:
        with lock:
            done_f.write(str(gene) + '\n')
    print('Working with gene ' + gene)
    results_one_gene = []
    results_one_gene_table = pd.DataFrame(
        {'gene': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    interval_indexes = df.index[df['gene_chr_start_end_strand'] == gene].tolist()
    for left_interval_index in range(interval_indexes[0], interval_indexes[len(interval_indexes) - 1]):
        seq1 = df.loc[left_interval_index, 'sequences'][:]
        seq1_indxd = df.loc[left_interval_index, 'sequences_indxd'][:]
        for right_interval_index in range(left_interval_index, interval_indexes[len(interval_indexes) - 1] + 1):
        #for right_interval_index in range(left_interval_index + 1, interval_indexes[len(interval_indexes) - 1] + 1):
            seq2 = df.loc[right_interval_index, 'sequences'][:]
            seq2_indxd = df.loc[right_interval_index, 'sequences_indxd'][:]
            if abs(df.loc[right_interval_index, 'chromStart'] - df.loc[left_interval_index, 'chromEnd']) \
                    <= panhandle_length_threshold:
                align = FindMinEnLocAlkmer(seq1, seq2, seq1_indxd, seq2_indxd, k, energy_threshold,
                                                     handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
                if align != 0:
                    results_one_gene.append([align, df.loc[left_interval_index, 'interval_chr_start_end_strand'],
                                             df.loc[right_interval_index, 'interval_chr_start_end_strand']])
    for result in results_one_gene:
        for alignment in result[0]:
            results_one_gene_table = results_one_gene_table.append(
                {'gene': gene, 'energy': alignment[0], 'interval1': result[1], 'interval2': result[2],
                 'start_al1': alignment[1], 'end_al1': alignment[2],
                 'start_al2': alignment[3], 'end_al2': alignment[4],
                 'alignment1': alignment[5], 'alignment2': alignment[6], 'structure': alignment[7]}, ignore_index=True)
    with open(path_to_ph, 'a') as f:
        with lock:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=False)



def Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix, strandness, annotation_file, path_to_ph):
    start_time = time.time()
    df = pd.read_csv(path_to_intervals, sep='\t')
    if 'sequences' in list(df.columns.values):
        seq = True
    else:
        seq = False
   

    # Check genes 
    if not('start_gene' in list(df.columns.values)):
        if(annotation_file == ''):
            print('I assume there is only one gene..')
            df.loc[:, 'start_gene'] = 1
            df.loc[:, 'end_gene'] = max(df.chromEnd) + 10 
        else:
            print('Attaching genes from annotation..')
            ncol = df.shape[1]
            cmd = 'bedtools intersect -a '+ path_to_intervals + ' -b '  + annotation_file + ' -s -wa -wb -header | awk -F"\t" \'{ if ($' + str(ncol + 3) + ' == "gene") { print } }\'| cut -f1-' + str(ncol) + ',' + str(ncol+4) + ',' + str(ncol+5)
            a = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            if sys.version_info[0] < 3: 
                from StringIO import StringIO
            else:
                from io import StringIO
            b = StringIO(a.communicate()[0].decode('utf-8'))
            df = pd.read_csv(b, sep="\t", header=None)
            if seq:
                df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'sequences', 'start_gene', 'end_gene']
            else:
                df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'start_gene', 'end_gene']

    # check chrom
    df.chrom.replace("_", "", regex=True, inplace=True) 

    # Make index for gene and interval
    df["gene_chr_start_end_strand"] = df.chrom + "_" + df.start_gene.map(str) + "_" + df.end_gene.map(str) + "_" + df.strand
    df["interval_chr_start_end_strand"] = df["chrom"] + "_" + df["chromStart"].map(str) + "_" + df[
        "chromEnd"].map(str) + "_" + df["strand"]
    df['chromStart'] = df['chromStart'].astype(int)
    df['chromEnd'] = df['chromEnd'].astype(int)
    df['start_gene'] = df['start_gene'].astype(int)
    df['end_gene'] = df['end_gene'].astype(int)
    
    # Order intervals
    df.sort_values(by=['gene_chr_start_end_strand', 'chromStart', 'chromEnd'], inplace=True)
    df.reset_index(inplace = True)
    
    # Attach sequences from the genome
    if not ('sequences' in list(df.columns.values)):
        print('Attaching sequences..')
        genome = Fasta(genome_file)
        GetSequencesForDF2 = partial(GetSequencesForDF, genome)
        df.loc[:, 'sequences'] = df.apply(GetSequencesForDF2, axis=1)
        if strandness:
            print("Making complement of minus strand..")
            df.loc[:, 'sequences'] = df.apply(MakeComplement, axis=1) 
        df.to_csv("../out/intervals_with_seqs.tsv", sep='\t', index=False)
    df.sequences = map(lambda x: x.upper(), df['sequences'])
    
    # Index sequences
    df["sequences_indxd"] = df['sequences'].apply(lambda x: Index_seq(x, k))
    df = df.loc[df.sequences_indxd != False]
    
    # Create headers for out files
    print("Creating files..")
    with open('../data/genes_done.txt', 'w') as done_f:
        done_f.write('Started alignment: \n')
    results_one_gene_table = pd.DataFrame(
        {'gene': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    with open(path_to_ph, 'w') as f:
        results_one_gene_table.to_csv(f, sep='\t', index=False, header=True)
    
    # Look for phs in parallel threads
    print('Start to align..')
    p = mp.Pool(processes=threads)
    m = mp.Manager()
    lock = m.Lock() # Allows multiple threads write into one file
    Find_panhandles_one_gene2 = partial(Find_panhandles_one_gene, lock, df, energy_threshold, handle_length_threshold,
                                        panhandle_length_threshold, k, need_suboptimal, kmers_stacking_matrix, path_to_ph)
    genes = df["gene_chr_start_end_strand"].unique()
    p.map(Find_panhandles_one_gene2, genes)
    p.close()
    p.join()
    print("all done!")
    print(time.time() - start_time)
    return (0)

def MakePretty(path_to_ph, annotation_file):
    df = pd.read_csv(path_to_ph, sep="\t")
    # Divide columns
    df[["start_gene", "end_gene"]] = df.gene.str.split('_', expand=True)[[1,2]].astype(int)
    df[["strand"]] = df.gene.str.split('_', expand=True)[[3]]
    df[["chr"]] = df.gene.str.split('_', expand=True)[[0]]
    df[["interval1_start", "interval1_end"]] = df.interval1.str.split('_', expand=True)[[1,2]].astype(int)
    df[["interval2_start", "interval2_end"]] = df.interval2.str.split('_', expand=True)[[1,2]].astype(int)

    # Make absolute coordinates
    df["panhandle_start"] = df.interval1_start + df.start_al1
    df["panhandle_left_hand"] = df.interval1_start + df.end_al1
    df["panhandle_right_hand"] = df.interval2_start + df.start_al2
    df["panhandle_end"] = df.interval2_start + df.end_al2

    # Calculate handle length
    df["al1_length"] = df.end_al1 - df.start_al1 + 1
    df["al2_length"] = df.end_al2 - df.start_al2 + 1

    # Remove broken phs (ph in one conserved interval that have end lefter than start)
    df = df.loc[df.panhandle_start < df.panhandle_right_hand]
    df = df.loc[df.panhandle_left_hand < df.panhandle_right_hand]
    df.drop(['interval1', 'interval2','start_al1','end_al1','start_al2','end_al2',
             'interval1_start','interval2_start', 'interval1_end','interval2_end', 'gene'], axis=1, inplace=True)

    # Add gene ids and names from annotation
    if annotation_file != '':
        print('Adding gene ids and names from annotation..')
        cmd = 'awk -F"\t" \'{ if ($3 == "gene") { print } }\' ' + annotation_file
        a = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        if sys.version_info[0] < 3: 
            from StringIO import StringIO
        else:
            from io import StringIO
        b = StringIO(a.communicate()[0].decode('utf-8'))
        anno = pd.read_csv(b, sep="\t", header=None)
        anno[["gene_id", "gene_name"]] = anno[8].str.split(';', expand=True)[[0,4]]
        anno.gene_id.replace("gene_id", "", regex=True, inplace=True)
        anno.gene_id.replace("\"", "", regex=True, inplace=True)
        anno.gene_name.replace("gene_name", "", regex=True, inplace=True)
        anno.gene_name.replace("\"", "", regex=True, inplace=True)

        anno = anno[[0,3,4,6,'gene_id','gene_name']]
        anno.columns = ['chr', 'start_gene', 'end_gene', 'strand','gene_id', 'gene_name']
        df = pd.merge(df, anno, on=['chr','start_gene','end_gene','strand'], how="inner")
        df.drop_duplicates(subset=["chr", "panhandle_start", "panhandle_left_hand",
                                   "panhandle_right_hand", "panhandle_end"], inplace=True)

    # Make alignments 5'-3' for all strands
    df['rev1'] = df.alignment1.apply(lambda x: x[::-1])
    df['rev2'] = df.alignment2.apply(lambda x: x[::-1])
    df.loc[df.strand == '-', 'alignment1'] = df.loc[df.strand == '-'].rev1
    df.loc[df.strand == '-', 'alignment2'] = df.loc[df.strand == '-'].rev2
    df.drop(['rev1', 'rev2'], axis=1)

    # Make structures 5'-3' for all strands
    df['rev1'] = df.alignment1.apply(lambda x: x[::-1])
    df['rev2'] = df.alignment2.apply(lambda x: x[::-1])
    df.loc[df.strand == '-', 'alignment1'] = df.loc[df.strand == '-'].rev1
    df.loc[df.strand == '-', 'alignment2'] = df.loc[df.strand == '-'].rev2
    df.drop(['rev1', 'rev2'], axis=1)

    # Add unique ids to phs
    df.sort_values(by=["chr", "panhandle_start", "panhandle_left_hand",
                               "panhandle_right_hand", "panhandle_end"], inplace=True)
    df["id"] = range(1, df.shape[0] + 1)
    df.to_csv("../out/panhandles_preprocessed.tsv", sep="\t", index=False)

def main(argv):
    #path_to_intervals = '../out/intervals_with_seqs.tsv'
    #genome_file = './../../../genomes/GRCh37.p13.genome.fa'
    path_to_intervals = ''
    genome_file = ''
    k = 5
    panhandle_length_threshold = 10000
    handle_length_threshold = 10
    threads = 5
    energy_threshold = -15
    need_suboptimal = True
    GT_threshold = 2
    strandness = True
    annotation_file = ''
    path_to_ph = "../out/panhandles.tsv"

    try:
        opts, args = getopt.getopt(argv, "h:i:g:k:p:a:t:e:u:d:s:n:o:",
                                   ["help=", "intervals=", "genome=", "k=", "panh_len_max", "handle_len_min", "threads",
                                    "energy_max", "need_subopt", "gt_threshold", "strand", "annotation", "out"])
    except getopt.GetoptError:
        print(
            'FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_length> -p <panhandle_len_max> ' +
            '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -s <strand>' +
            '-n <annotation> -o <out>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_length> -p <panhandle_len_max> ' +
                '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -s <strand>' +
                '-n <annotation> -o <out>')
            sys.exit()
        elif opt in ("-i", "--intervals"):
            path_to_intervals = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-k", "--kmer"):
            k = int(arg)
        elif opt in ("-p", "--panh_len_max"):
            panhandle_length_threshold = int(arg)
        elif opt in ("-a", "--handle_len_min"):
            handle_length_threshold = int(arg)
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-e", "--energy_max"):
            energy_threshold = float(arg)
        elif opt in ("-u", "--need_subopt"):
            need_suboptimal = bool(arg)
        elif opt in ("-d", "--gt_threshold"):
            GT_threshold = int(arg)
        elif opt in ("-s", "--strand"):
            strandness = bool(arg)
        elif opt in ("-n", "--annotation"):
            annotation_file = arg
        elif opt in ("-0", "--out"):
            path_to_ph = arg
    kmers_stacking_matrix = load("../lib/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy")
    Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold,
                    panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix, strandness, annotation_file, path_to_ph)
    MakePretty(path_to_ph, annotation_file)


if __name__ == '__main__':
    main(sys.argv[1:])

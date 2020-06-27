#!/usr/bin/env python2
from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort
from math import ceil
import re, sys, getopt, itertools, binascii, time, subprocess, os
from functools import partial
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append('../../tools/')
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
        return(str(Seq(row['sequences']).reverse_complement()))
    else:
        return(row['sequences'])

def Find_panhandles_one_gene(lock, df, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                             need_suboptimal, kmers_stacking_matrix, path_to_ph, gene):
    with open(path_to_ph + '/progress.txt', 'a') as done_f:
        with lock:
            done_f.write(str(gene) + '\n')
    #print('Working with gene ' + gene)
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
    with open(path_to_ph + '/panhandles.tsv', 'a') as f:
        with lock:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=False)

def Find_panhandles_one_row(lock, df, energy_threshold, handle_length_threshold, k,
                             need_suboptimal, kmers_stacking_matrix, path_to_ph, row):
    with open(path_to_ph + '/progress.txt', 'a') as done_f:
        with lock:
            done_f.write(str(row) + '\n')
    #print('Working with row ' + str(row))
    results_one_row = []
    results_one_row_table = pd.DataFrame(
        {'row': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    
    seq1 = df.loc[0, 'sequences'][:]
    seq1_indxd = df.loc[0, 'sequences_indxd'][:]
    seq2 = df.loc[row, 'sequences'][:]
    seq2_indxd = df.loc[row, 'sequences_indxd'][:]
    align = FindMinEnLocAlkmer(seq1, seq2, seq1_indxd, seq2_indxd, k, energy_threshold,
        handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
    if align != 0:
        results_one_row.append([align, df.loc[0, 'interval_chr_start_end_strand'],
            df.loc[row, 'interval_chr_start_end_strand']])
    for result in results_one_row:
        for alignment in result[0]:
            results_one_row_table = results_one_row_table.append(
                {'row': row, 'energy': alignment[0], 'interval1': result[1], 'interval2': result[2],
                 'start_al1': alignment[1], 'end_al1': alignment[2],
                 'start_al2': alignment[3], 'end_al2': alignment[4],
                 'alignment1': alignment[5], 'alignment2': alignment[6], 'structure': alignment[7]}, ignore_index=True)
    with open(path_to_ph + '/panhandles.tsv', 'a') as f:
        with lock:
            results_one_row_table.to_csv(f, sep='\t', index=False, header=False)
    

def Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix, strandness, annotation_file, 
                    path_to_ph, first_to_all):
    start_time = time.time()
    df = pd.read_csv(path_to_intervals, sep='\t')
    if 'sequences' in list(df.columns.values):
        seq = True
    else:
        seq = False
   

    one_gene = False
    # Check genes 
    if not('start_gene' in list(df.columns.values)):
        if(annotation_file == ''):
            print('I assume there is only one gene..')
            df.loc[:, 'start_gene'] = 1
            df.loc[:, 'end_gene'] = max(df.chromEnd) + 10 
            one_gene = True
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

    # Make index for gene 
    if one_gene:
        df["gene_chr_start_end_strand"] = "chrSmth" + "_" + df.start_gene.map(str) + "_" + df.end_gene.map(str) + "_" + df.strand
    else:
        df["gene_chr_start_end_strand"] = df.chrom + "_" + df.start_gene.map(str) + "_" + df.end_gene.map(str) + "_" + df.strand

    # Make index for interval
    df["interval_chr_start_end_strand"] = df["chrom"] + "_" + df["chromStart"].map(str) + "_" + df[
        "chromEnd"].map(str) + "_" + df["strand"]
    df['chromStart'] = df['chromStart'].astype(int)
    df['chromEnd'] = df['chromEnd'].astype(int)
    df['start_gene'] = df['start_gene'].astype(int)
    df['end_gene'] = df['end_gene'].astype(int)
    
    # Attach sequences from the genome
    if not ('sequences' in list(df.columns.values)):
        print('Attaching sequences..')
        genome = Fasta(genome_file)
        GetSequencesForDF2 = partial(GetSequencesForDF, genome)
        df.loc[:, 'sequences'] = df.apply(GetSequencesForDF2, axis=1)
        if strandness:
            print("Making reverse complement of minus strand..")
            df.loc[:, 'sequences'] = df.apply(MakeComplement, axis=1) 
        df.to_csv("../data/intervals_with_seqs.tsv", sep='\t', index=False)
    
    # Check sequences are all upper case    
    df.sequences = map(lambda x: x.upper(), df['sequences'])
    
    # Check sequences are all from DNA, not RNA
    df.sequences = df.sequences.replace("U", "T", regex=True) 

    # Check if sequences have non ATGC symbols
    df["wrong_symbol"] = df['sequences'].apply(lambda x: any(c not in 'ATGC' for c in x))
    if any(df.wrong_symbol == True):
        print('I have removed ' + str(df.loc[df.wrong_symbol == True].shape[0]) + ' rows with non ATGC symbols')
        df = df.loc[df.wrong_symbol != True]

    # Index sequences
    df["sequences_indxd"] = df['sequences'].apply(lambda x: Index_seq(x, k))
    df = df.loc[df.sequences_indxd != False]
    
    # Look for phs in parallel threads
    if first_to_all == False:
        # Create headers for out files
        print("Creating files..")
        with open(path_to_ph + '/progress.txt', 'w') as done_f:
            done_f.write('Started alignment: \n')
        results_one_gene_table = pd.DataFrame(
            {'gene': [], 'energy': [],
             'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
             'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
        with open(path_to_ph + '/panhandles.tsv', 'w') as f:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=True)
        # Order intervals
        df.sort_values(by=['gene_chr_start_end_strand', 'chromStart', 'chromEnd'], inplace=True)
        df.reset_index(inplace = True)
        # Work
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
    elif first_to_all == True:
        print('I will compare only the first sequence to all the others!')
        print("Creating files..")
        with open(path_to_ph + '/progress.txt', 'w') as done_f:
            done_f.write('Started alignment: \n')
        results_one_row_table = pd.DataFrame(
            {'row': [], 'energy': [],
             'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
             'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
        with open(path_to_ph + '/panhandles.tsv', 'w') as f:
            results_one_row_table.to_csv(f, sep='\t', index=False, header=True)
        # Work
        print('Start to align..') 
        p = mp.Pool(processes=threads)
        m = mp.Manager()
        lock = m.Lock() # Allows multiple threads write into one file
        Find_panhandles_one_row2 = partial(Find_panhandles_one_row, lock, df, energy_threshold, handle_length_threshold, k,
            need_suboptimal, kmers_stacking_matrix, path_to_ph)
        rows = range(1, df.shape[0])
        p.map(Find_panhandles_one_row2, rows)
        p.close()
        p.join()
    print("all done!")
    print(time.time() - start_time)
    return (0)

def MakePretty(path_to_ph, annotation_file, RNA_RNA_interaction):
    if RNA_RNA_interaction == True:
        print('I treat the results as RNA-RNA interaction')
    df = pd.read_csv(path_to_ph + '/panhandles.tsv', sep="\t")
    if df.shape[0] == 0:
        print('No structures found!')
    else:
        # Divide columns
        if  'gene' in list(df.columns.values):
            df[["start_gene", "end_gene"]] = df.gene.str.split('_', expand=True)[[1,2]].astype(int)

        if(RNA_RNA_interaction == False): 
            df[["strand"]] = df.interval1.str.split('_', expand=True)[[3]]
            df[["chr"]] = df.interval1.str.split('_', expand=True)[[0]]
        else:
            df[["strand_1"]] = df.interval1.str.split('_', expand=True)[[3]]
            df[["chr_1"]] = df.interval1.str.split('_', expand=True)[[0]]
            df[["strand_2"]] = df.interval2.str.split('_', expand=True)[[3]]
            df[["chr_2"]] = df.interval2.str.split('_', expand=True)[[0]]

        df[["interval1_start", "interval1_end"]] = df.interval1.str.split('_', expand=True)[[1,2]].astype(int)
        df[["interval2_start", "interval2_end"]] = df.interval2.str.split('_', expand=True)[[1,2]].astype(int)

        # Make absolute coordinates
        if(RNA_RNA_interaction == False): 
            # + strand
            df["panhandle_start"] = df.interval1_start + df.start_al1
            df["panhandle_left_hand"] = df.interval1_start + df.end_al1
            df["panhandle_right_hand"] = df.interval2_start + df.start_al2
            df["panhandle_end"] = df.interval2_start + df.end_al2
            # - strand (reverse complement)
            df.loc[df.strand == '-', 'panhandle_left_hand'] = df.loc[df.strand == '-'].interval1_end - df.loc[df.strand == '-'].start_al1
            df.loc[df.strand == '-', 'panhandle_start'] = df.loc[df.strand == '-'].interval1_end - df.loc[df.strand == '-'].end_al1
            df.loc[df.strand == '-', 'panhandle_end'] = df.loc[df.strand == '-'].interval2_end - df.loc[df.strand == '-'].start_al2
            df.loc[df.strand == '-', 'panhandle_right_hand'] = df.loc[df.strand == '-'].interval2_end - df.loc[df.strand == '-'].end_al2

            df = df.loc[~((df.panhandle_left_hand > df.panhandle_right_hand) & (df.panhandle_start < df.panhandle_end))]

            x = df.loc[df.panhandle_start > df.panhandle_end]
            y = df.loc[~(df.panhandle_start > df.panhandle_end)]
            if x.shape[0] != 0:
                x['al1'] = x[['alignment2']]
                x['alignment2'] = x[['alignment1']]
                x['alignment1'] = x[['al1']]
                x['lh'] = x[['panhandle_left_hand']]
                x['panhandle_left_hand'] = x[['panhandle_end']]
                x['panhandle_end'] = x[['lh']]
                x['st'] = x[['panhandle_start']]
                x['panhandle_start'] = x[['panhandle_right_hand']]
                x['panhandle_right_hand'] = x[['st']]
                x.drop(['al1', 'lh', 'st'], inplace = True, axis = 1)
                df = y.append(x)

        else:
            # + strand
            df.loc[df.strand_1 == '+',"panhandle_start"] = df.loc[df.strand_1 == '+'].interval1_start + df.loc[df.strand_1 == '+'].start_al1
            df.loc[df.strand_1 == '+',"panhandle_left_hand"] = df.loc[df.strand_1 == '+'].interval1_start + df.loc[df.strand_1 == '+'].end_al1
            df.loc[df.strand_2 == '+', "panhandle_right_hand"] = df.loc[df.strand_2 == '+'].interval2_start + df.loc[df.strand_2 == '+'].start_al2
            df.loc[df.strand_2 == '+', "panhandle_end"] = df.loc[df.strand_2 == '+'].interval2_start + df.loc[df.strand_2 == '+'].end_al2
            # - strand (reverse complement)
            df.loc[df.strand_1 == '-', 'panhandle_left_hand'] = df.loc[df.strand_1 == '-'].interval1_end - df.loc[df.strand_1 == '-'].start_al1
            df.loc[df.strand_1 == '-', 'panhandle_start'] = df.loc[df.strand_1 == '-'].interval1_end - df.loc[df.strand_1 == '-'].end_al1
            df.loc[df.strand_2 == '-', 'panhandle_end'] = df.loc[df.strand_2 == '-'].interval2_end - df.loc[df.strand_2 == '-'].start_al2
            df.loc[df.strand_2 == '-', 'panhandle_right_hand'] = df.loc[df.strand_2 == '-'].interval2_end - df.loc[df.strand_2 == '-'].end_al2



        # Calculate handle length
        df["al1_length"] = df.panhandle_left_hand - df.panhandle_start + 1
        df["al2_length"] = df.panhandle_end - df.panhandle_right_hand + 1

        # Remove broken phs (ph in one conserved interval that have end lefter than start)
        #if ((first_to_all == False) & ()):
            #df = df.loc[df.panhandle_start < df.panhandle_right_hand]
            #df = df.loc[df.panhandle_left_hand < df.panhandle_right_hand]
        df.drop(['interval1', 'interval2','start_al1','end_al1','start_al2','end_al2',
                 'interval1_start','interval2_start', 'interval1_end','interval2_end'], axis=1, inplace=True)
        if  'gene' in list(df.columns.values):
            df.drop(['gene'], axis=1, inplace=True)

        # Add gene ids and names from annotation
        if ((annotation_file != '') & ('start_gene' in list(df.columns.values))):
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
        #df['rev1'] = df.alignment1.apply(lambda x: x[::-1])
        #df['rev2'] = df.alignment2.apply(lambda x: x[::-1])
        #df.loc[df.strand == '-', 'alignment1'] = df.loc[df.strand == '-'].rev1
        #df.loc[df.strand == '-', 'alignment2'] = df.loc[df.strand == '-'].rev2
        #df.drop(['rev1', 'rev2'], axis=1)

        # Make structures 5'-3' for all strands

        # Add unique ids to phs
        if RNA_RNA_interaction == False:
            df.sort_values(by=["chr", "panhandle_start", "panhandle_left_hand",
                                       "panhandle_right_hand", "panhandle_end"], inplace=True)
        else:
            df.sort_values(by=["chr_1", "panhandle_start", "panhandle_left_hand",
                                       "panhandle_right_hand", "panhandle_end"], inplace=True)
        df["id"] = range(1, df.shape[0] + 1)
        df.to_csv(path_to_ph + '/panhandles_preprocessed.tsv', sep="\t", index=False)

def main(argv):
    path_to_intervals = '../data/conin_python_long_coding_final.tsv'
    genome_file = ''
    k = 5
    panhandle_length_threshold = 10000
    handle_length_threshold = 10
    threads = 1
    energy_threshold = -15
    need_suboptimal = True
    GT_threshold = 2
    strandness = True
    annotation_file = ''
    path_to_ph = "../data/panhandles"
    first_to_all = False
    RNA_RNA_interaction = False 

    try:
        opts, args = getopt.getopt(argv, "h:i:g:k:p:a:t:e:u:d:s:n:o:r:c:",
                                   ["help=", "intervals=", "genome=", "k=", "panh_len_max", "handle_len_min", "threads",
                                    "energy_max", "need_subopt", "gt_threshold", "strand", "annotation", "out", "first_to_all", "RNA_RNA_interaction"])
    except getopt.GetoptError:
        print(
            'FindPanhandles.py -i <intervals.bed> -g <genome.fa> -k <kmer_length> -p <panhandle_len_max> ' +
            '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -s <strand> -r <first_to_all> ' +
            '-c <RNA_RNA_interaction> -n <annotation.gtf> -o <out_folder>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'FindPanhandles.py -i <intervals.bed> -g <genome.fa> -k <kmer_length> -p <panhandle_len_max> ' +
                '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -s <strand True> -r <first_to_all False> ' +
                '-c <RNA_RNA_interaction False> -n <annotation.gtf> -o <out_folder>')
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
            need_suboptimal = arg
        elif opt in ("-d", "--gt_threshold"):
            GT_threshold = int(arg)
        elif opt in ("-s", "--strand"):
            strandness = arg
        elif opt in ("-n", "--annotation"):
            annotation_file = arg
        elif opt in ("-o", "--out"):
            path_to_ph = arg
        elif opt in ("-r", "--first_to_all"):
            first_to_all = arg   
        elif opt in ("-c", "--RNA_RNA_interaction"):
            RNA_RNA_interaction = arg                
    if first_to_all == 'False':
        first_to_all = False
    elif first_to_all == 'True':
        first_to_all = True
    if RNA_RNA_interaction == 'False':
        RNA_RNA_interaction = False
    elif RNA_RNA_interaction == 'True':
        RNA_RNA_interaction = True   
    if need_suboptimal == 'False':
        need_suboptimal = False
    elif need_suboptimal == 'True':
        need_suboptimal = True
    if strandness == 'False':
        strandness = False
    elif strandness == 'True':
        strandness = True

    kmers_stacking_matrix = load("../data/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy")

    Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold,
                    panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, 
                    kmers_stacking_matrix, strandness, 
                    annotation_file, path_to_ph,
                    first_to_all)
    MakePretty(path_to_ph, annotation_file, RNA_RNA_interaction)


if __name__ == '__main__':
    main(sys.argv[1:])

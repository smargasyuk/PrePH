#!/usr/bin/env python2
import re, sys, getopt, itertools, binascii, time, subprocess
from functools import partial
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq

def split_overlap(seq, size, overlap):        
    for i in range(0, len(seq) - overlap, size - overlap):  
        yield(seq[i:i + size])
def split_overlap_i(seq, size, overlap):        
    for i in range(0, len(seq) - overlap, size - overlap):  
        yield(i)
def find_end(row):
    return(row['chromStart'] + len(row['sequences']) - 1)

def MakeDataFrame(size, overlap, record):
    seq = str(record.seq)
    seqs = [x for x in split_overlap(seq, size, overlap)]
    df = pd.DataFrame(seqs, columns=['sequences'])
    df['chrom'] = record.id
    df['chromStart'] = [x+1 for x in split_overlap_i(seq, size, overlap)] 
    df['chromEnd'] = df.apply(find_end, axis = 1)
    df['name'] = range(1,df.shape[0]+1)
    df['score'] = 1
    df['strand'] = '+'
    df = df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'sequences']]
    return(df)

def main(argv):
    path_to_input = ''
    path_to_output = ''
    size = 1000
    overlap = 30

    try:
        opts, args = getopt.getopt(argv, "i:o:s:v:",
                                   ["help=", "input=", "output=", "size=", "overlap="])
    except getopt.GetoptError:
        print('MakeBedForVirusGenome.py -i <input.fa> -o <output.bed> -s <size> -v <overlap>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('MakeBedForVirusGenome.py -i <input.fa> -o <output.bed> -s <size> -v <overlap>')
            sys.exit()
        elif opt in ("-i", "--input"):
            path_to_input = arg
        elif opt in ("-o", "--output"):
            path_to_output = arg
        elif opt in ("-s", "--size"):
            size = int(arg)
        elif opt in ("-v", "--overlap"):
            overlap = int(arg)

    
    records = list(SeqIO.parse(path_to_input, "fasta"))

    MakeDataFrame2 = partial(MakeDataFrame, size,  overlap)
    dfs = map(MakeDataFrame2, records)
    df = pd.concat(dfs)
    df.to_csv(path_to_output, sep='\t', index=False, header=True)





    



if __name__ == '__main__':
    main(sys.argv[1:])

![DOI:250842115](https://zenodo.org/badge/doi/10.5281/zenodo.4601044.svg)

# PrePH
Predict PanHandles
PrePH is a set of Python scripts for finding Pairs of Complementary regions 

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# System Requirements
## OS Requirements
This tool is supported for macOS. The tool has been tested on the following systems:
- Ubuntu 18.04.4
- CentOS Linux 7

## Python Dependencies
- Python version = 2.7
- pandas 
- numpy 
- getopt
- functools
- pyfaidx 
- multiprocessing 
- Bio == 1.76

## External Dependncies
- bedtools

## Installation
`git clone https://github.com/kalmSveta/PrePH.git`

`export PATH=$PATH:path_to_PrePH/src`

Typical install time on a "normal" desktop computer: less than 3s

Now you can run the scripts from any place

Always use GLOBAL paths to input and output files

## Step 0 - Precalculates stacking energies matrix for kmers
Run `PrecalculateStackingEnergeis.py -k <kmer_length> -g <gt_amount_in_kmer_max>` 

This script needs to be run only once before usage. 
Parameters:
- -k <kmer_length> - is a minimal length of consicutive stacked nt pairs. Must be the same as used for `./FindPanhandles.py` Recommended k = 5, default 5
- -g <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is g = 2, default 2

### Example:
`PrecalculateStackingEnergeis.py -k 5 -g 2`

Expected output: in PrePH/data 3 files will be created:
52mers_stacking_energy_binary.npy, 52mers_stacking_energy_no$.npy, kmers_list_5mers_no$.txt


## Find Pairs of Complementary Regions in two sequences
If you need to find complementary regions in just TWO sequences use this script. If you need to find PCCRs in the wholw genome or subset of it, use the pipeline below.

Run `fold.py -f <first_seq> -s <second_seq> -k <kmer_length> -a <handle_len_min> -e <energy_max, kcal/mol> -u <need_subopt> -d <gt_threshold>`

Parameters:

- -f first DNA sequence, default ''
- -s second DNA sequence, default ''
- -k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for `PrecalculateStackingEnergeis.py`. Recommended k = 5, default 5
- -a <handle_len_min> - minimal length of handles. Recommended = 10, default 10
- -e <energy_max> - maximum energy in kcal/mol. Recommended = -15, default -15
- -u <need_suboptimal> - if True, will try to find suboptimal structures. If False, will return only MFE structure. Default True
- -d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2

### Example
`./fold.py -f AAAGGGC -s AAAGCCCAAAAAACCTTT -k 3 -a 3 -e -1 -u True -d 2`

Expected output:
[(-10.0, 3, 6, 3, 6, 'GGGC', 'GCCC', '(((())))'), (-7.2, 0, 4, 13, 17, 'AAAGG', 'CCTTT', '((((()))))')]


## The pipeline to find PCCRs in the whole genome or subset of it
-----------------------------------------------------------------

## Step 1 - Select conserved intronic intervals (CII)

This script can be used to select conserved intronic intervals from the whole genome

Run `SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks> -t <gene_type of coding/noncoding/all>`

- -a <annotation> - GLOBAL path to genome annotation in gtf format, default ''
- -c <cons_regions> - GLOBAL path to phastConsElements file, default ''
- -l <handle_len_min> - minimum length of CCR, default 10
- -f <flanks> - length of flanks that can intersect CDS, default 10 
- -t <gene_type of coding/noncoding/all> - which genes to select, default 'coding'
  
The output files will be stored in PrePH/data/ directory: 
- conin_python_long_coding_final.tsv - CII. Ths file can be used as input for FindPanhandles.py
- introns_with_flanks_python.bed - introns with flanking regions in CDS
- genes.bed - genes selected from annotation in bed format

### Example:
This is a toy example with a subset of gencode hg19 annotation of only one gene 

`SelectIntervals.py -a global_bath_to_PrePH/lib/example_gencode.v19.annotation.gtf  -c global_bath_to_PrePH/lib/example_phastConsElements100way.txt  -l 10 -f 10 -t coding`

Expected output: 3 output files in PrePH/data/
Expected run time: less than 9s

## Step 2 - Predicts panhandles with handles in the intervals 
Run `FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max>  -a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -o <out> -n <annotation> -s <strandness>  -d <gt_amount_in_kmer_max> -r <first_to_all>`


This script can be run in 3 modes:

1) if -i input file has gene coordinates or -n annotation_file is provided, then PrePH will compare seuqences from same genes with each other 
2) if -i input file does not have gene coordinates and -n annotaion_file is not provided, then PrePH will compare each sequence with each
3) if -r first_to_all is True, then regardless of input and annotation_file PrePH will compare the first sequence with all the others

Parameters:
- -i <intervals_df> - GLOBAL path to input file of intervals in bed 6 format (with possible additional 3 columns)  - tab-separated file with header like this:
**The header is compulsory!**

| chrom | chromStart | chromEnd | name | score | strand | sequences | start_gene | end_gene |
| :---: |   :---:    |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
|chr1  |136219| 136233| 1 | 1 | + | GGCTTTGATAAAAA |  135223 |  138932 |
|chr1 | 230 | 243 | 2 | 1 | - | TTTTTATAAAGCC | 105 | 310 | 
|chr1 | 136330 | 136343 | 3 | 1 | + | GGCCAGCAGATGG | 135223 | 138932 |
|chr1 | 1005 | 1028  | 4 | 1 | - | TGACAAACCACAGGACACTACAC | 105 | 310 | 

This file can be taken from step 1.

**if column 'sequences' is absent:**
`-g <genome.fa>` GLOBAL path to genome must be provided in fasta format. It should have sequences which include interval sequences. E.g. a whole human genome. 'sequences' colunm will be extracted automatically

The comparison of seuqences is performed only inside one gene, that is why PrePh needs to know gene coordinates. All sequences in one gene will be compared in a pairwise manner to find complementary regions between them

**if colums gene_start and gene_end are absent:**
but `-n <annotation>` is provided, the genes will be outomatically identified from the genome annotation.
If `-n <annotation>` is also absent, PrePH assumes that all the sequences belong to one gene and will create artificial coordinates for the gene.   
  

- -k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for `PrecalculateStackingEnergeis.py`. Recommended k = 5, default 5
- -p <panhandle_len_max> - maximum length of looped out region of panhandle. Recommened = 10000, default 10000
- -a <handle_len_min> - minimal length of handles. Recommended = 10, default 10
- -t <threads> - number of threads to run in parallel, default 1
- -e <energy_max> - maximum energy in kcal/mol. Recommended = -15, default -15
- -u <need_suboptimal> - if True, will try to find suboptimal structures. If False, will return only MFE structure for each pair of compared CII. Default True
- -d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2, default 2
- -s <strandness> - If True (recommended), the sequences extracted from the genome with strand == '-' will be reverse complemented. If sequences are alredy present in the input_file, this parameter is ignored. Default True
- -n <annotation> - genome annotation in gtf or gff format, default ''
- -r <first_to_all> - PrePH will compare ONLY one sequence (FIRST) to all the others. In this case genes are ignored. Default False
- -o <out> - GLOBAL path to output directory, default NA
- -c <RNA_RNA_interaction> - this parameter is in the test mode. Use it to find interamolecular interactions. Defaut False.   
  
The final output will be stored in file called \<out>/panhandles_preprocessed.tsv


### Test example
To test the scripts run:

`FindPanhandles.py -i ../lib/example_intervals.bed -k 5 -p 10000 -a 10 -t 5 -e -15 -u True -d 2 -s True -o output`

Expected output: panhandles_preprocessed.tsv output files in PrePH/data/
The ouput file has the following columns:
- alignment1, alignment2 - sequences of CCRs
- energy - energy of PCCR in kcal/mol
- structure - PCCR structure in dot-bracket notation
- start_gene, end_gene - coordinates of gene
- strand
- chr     
- panhandle_start, panhandle_left_hand - coordintates of the left CCR
- panhandle_right_hand, panhandle_end - coordinates of the right CCR
- al1_length, al2_length - CCRs length
- id - unique id of PCCR

Expected run time: less than 1s

## Make bed input file from short genome sequence
If you need to find panhandles in a short (e.g virus) genome, you need to divide the sequence into a set of smaller overlapping sequences to reduce memory loading. Recommended -s = 1000, -v = 30

Run `MakeBedForVirusGenome.py -i <input.fa> -o <output.bed> -s <size> -v <overlap>`

- -i <input.fa> - sequence in fasta format. All chomoseomes must be in one sequence
- -o <output.bed> - GLOBAL path to output bed file
- -s <size> - length of output sequences
- -v <overlap> - length of overlaps



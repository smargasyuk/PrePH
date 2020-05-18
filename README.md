# PrePH
Predict PanHandles

Find Pairs of Complementary regions 

Python version = 2.7

## Installation
`git clone https://github.com/kalmSveta/PrePH.git`

`export PATH=$PATH:path_to_PrePH/src`

Now you can run the scripts from any place

Always use GLOBAL paths to input and output files


## Step 0 - Precalculates stacking energies matrix for kmers
Run `./PrecalculateStackingEnergeis.py -k <kmer_length> -g <gt_amount_in_kmer_max>` 

This script needs to be run only once before usage. 
Parameters:
- -k <kmer_length> - is a minimal length of consicutive stacked nt pairs. Must be the same as used for `./FindPanhandles.py` Recommended k = 5
- -g <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is g = 2

## Step 1 - Select conserved intronic intervals (CII)
Run `SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks>`

- -a <annotation> - GLOBAL path to genome annotation in gtf format
- -c <cons_regions> - GLOBAL path to phastConsElements file
- -l <handle_len_min> - minimum length of CCR
- -f <flanks> - length of flanks that can intersect CDS
  
The output files will be stored in ../data/ directory: 
- conin_python_long_coding_final.tsv - CII. Ths file can be used as input for FindPanhandles.py
- introns_with_flanks_python.bed - introns with flanking regions in CDS
- coding_genes.bed - coding genes selected from annotation in bed format

## Step 2 - Predicts panhandles with handles in the intervals 
Run `FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max>  -a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -o <out> -n <annotation> -s <strandness>  -d <gt_amount_in_kmer_max> -r <first_to_all>`

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
  

- -k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for `PrecalculateStackingEnergeis.py`. Recommended k = 5
- -p <panhandle_len_max> - maximum length of looped out region of panhandle. Recommened = 10000
- -a <handle_len_min> - minimal length of handles. Recommended = 10
- -t <threads> - number of threads to run in parallel
- -e <energy_max> - maximum energy in kcal/mol. Recommended = -15
- -u <need_suboptimal> - if True, will try to find suboptimal structures
- -d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2
- -s <strandness> - account for the strand
- -n <annotation> - genome annotation in gtf or gff format
- -r <first_to_all> - PrePH will compare ONLY one sequence (FIRST) to all the others. In this case genes are ignored
- -o <out> - GLOBAL path to output file 
  
The final output will be stored in file called \<out>_preprocessed


### Test example
To test the scripts run:

`PrecalculateStackingEnergeis.py -k 5 -g 2` 

`FindPanhandles.py -i ../lib/exammple.bed -k 5 -p 10000 -a 10 -t 5 -e -15 -u True -d 2 -s True -o output`

## Find Pair of Complementary Regions in two sequences
Run `fold.py -f <first_seq> -s <second_seq> -k <kmer_lentgh> -a <handle_len_min> -e <energy_max> -u <need_subopt> -d <gt_threshold>`

Parameters are the same as above, except:

- -f first DNA sequence
- -s second DNA sequence

### Example
`fold.py -f AAAGGGC -s AAAGCCCAAAAAA -k 3 -a 3 -e -1 -u True -d 2`


## Make bed input file from short genome sequence
If you need to find panhandles in a short (e.g virus) genome, you need to divide the sequence into a set of smaller overlapping sequences to reduce memory loading. Recommended -s = 1000, -v = 30

Run `MakeBedForVirusGenome.py -i <input.fa> -o <output.bed> -s <size> -v <overlap>`

- -i <input.fa> - sequence in fasta format. All chomoseomes must be in one sequence
- -o <output.bed> - GLOBAL path to output bed file
- -s <size> - length of output sequences
- -v <overlap> - length of overlaps



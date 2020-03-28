# PrePH
Predict PanHandles

Find Pairs of Complementary regions 

All scripts should be run from PrePH/src directory

## Step 0 - Precalculates stacking energies matrix for kmers
Run `./PrecalculateStackingEnergeis.py -k <kmer_length> -g <gt_amount_in_kmer_max>` 

This script needs to be run only once before usage. 
Parameters:
- -k <kmer_length> - is a minimal length of consicutive stacked nt pairs. Must be the same as used for `./FindPanhandles.py` Recommended k = 5
- -g <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is g = 2

## Step 1 - Predicts panhandles with handles in the intervals
Run `FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max>  -a <handle_len_min> -t <threads> -e <energy_max> -s <have_seqs_in_input_df> -u <need_suboptimal>`

Parameters:
- -i <intervals_df> - an input file of intervals in bed 6 format (with possible additional 3 columns)  - tab-separated file with header like this:
**The header is compulsory!**

| chrom | chromStart | chromEnd | name | score | strand | sequences | gene_start | gene_end |
| :---: |   :---:    |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
|chr1  |136219| 136233| 1 | 1 | + | GGCTTTGATAAAAA |  135223 |  138932 |
|chr1 | 230 | 243 | 2 | 1 | - | TTTTTATAAAGCC | 105 | 310 | 
|chr1 | 136330 | 136343 | 3 | 1 | + | GGCCAGCAGATGG | 135223 | 138932 |
|chr1 | 1005 | 1028  | 4 | 1 | - | TGACAAACCACAGGACACTACAC | 105 | 310 | 

**if column 'sequences' is absent:**
`-g <genome.fa>` must be provided in fasta format. It should have sequences which include interval sequences. E.g. a whole human genome. 'sequences' colunm will be extracted automatically


The comparison of seuqences is performed only inside one gene, that is why PrePh needs to know gene coordinates
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
- -o <out> - path to output file 
  
The final output will be stored in file called <out>_preprocessed

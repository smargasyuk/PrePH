# PrePH
Predict PanHandles

Find Pairs of Complementary regions 

All scripts should be run from PrePH/src directory

## Step 0 - Precalculates stacking energies matrix for kmers.
Run `./PrecalculateStackingEnergeis.py` 

This script needs to be run only once before usage. 
Parameters:
- -k <kmer_length> - is a minimal length of consicutive stacked nt pairs. Must be the same as used for `./FindPanhandles.py` Recommended k = 5
- -g <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is g = 2

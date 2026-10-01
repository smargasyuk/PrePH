# PrePH

Predict PanHandles

PrePH is a set of Python scripts for finding Pairs of Complementary regions.

Forked from [kalmSveta/PrePH](https://github.com/kalmSveta/PrePH) by Svetlana Kalmykova.

    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Setup

## Installation

Install as a package:

`pip install git+https://github.com/smargasyuk/PrePH.git`

Or install in development mode:

```bash
git clone https://github.com/smargasyuk/PrePH.git
cd PrePH
pip install -e .
```

## k-mer stacking energy precalculation
Run `preph-precalculate-energies -k <kmer_length> -g <gt_amount_in_kmer_max>`

This script needs to be run only once before usage. 
Parameters:
- `-k <kmer_length>` - is a minimal length of consecutive stacked nt pairs. Must be the same as used for other scripts. Recommended `k = 5`, default `5`.
- `-g <gt_amount_in_kmer_max>` - maximal number of GT pairs in kmer. Recommended for `k = 5` is `g = 2`, default `2`.

### Example:
`preph-precalculate-energies -k 5 -g 2`

Expected output: in `$HOME/.local/share/preph/` 3 files will be created:
```
52mers_stacking_energy_binary.npy
52mers_stacking_energy_no$.npy
kmers_list_5mers_no$.txt
```

# Find Pairs of Complementary Regions (PCRs) in two sequences
If you need to find complementary regions in just two sequences use this script. If you need to find PCRs in the whole genome or subset of it, use the pipeline below.

Run `preph-fold-raw -f <first_seq> -s <second_seq> -k <kmer_length> -a <handle_len_min> -e <energy_max, kcal/mol> -u <need_subopt> -d <gt_threshold>`

Parameters:

- `-f` first DNA sequence, default ''.
- `-s` second DNA sequence, default ''.
- `-k <kmer_length>` is a minimal length of consicutive stacked nt pairs. Must be the same as used for `preph-precalculate-energies`. Recommended `k = 5`, default `5`.
- `-a <handle_len_min>` - minimal length of handles. Recommended = `10`, default `10`.
- `-e <energy_max>` - maximum energy in kcal/mol. Recommended = `-15`, default `-15`.
- `-u <need_suboptimal>` - if True, will try to find suboptimal structures. If False, will return only MFE structure. Default True.
- `-d <gt_amount_in_kmer_max>` - maximum number of GT pairs in kmer. Recommended for `k = 5 `is `2.`

### Example:
`preph-fold-raw -f AAAGGGC -s AAAGCCCAAAAAACCTTT -k 3 -a 3 -e -1 -u True -d 2`

Expected output:
```text
[(-10.0, 3, 6, 3, 6, 'GGGC', 'GCCC', '(((())))'),
 (-7.2, 0, 4, 13, 17, 'AAAGG', 'CCTTT', '((((()))))')]
```


# Find PCRs in genomic intervals

Predicts structures and reports them in absolute genomic coordinates. The subworkflow assumes that global genomic coordinates of input intervals are available, allowing predicted structures to be mapped to their genomic locations and displayed.

## 1. Predict panhandles with handles in genomic intervals

Run `preph-genomic-fold --input <intervals_table> --output <output_table>`

### Input

`--input <intervals_table> ` — path to the input file containing intervals in BED6 format with two additional columns. Unlike standard BED6, the coordinates are **1-based and inclusive**. The file must be tab-separated and include a header with the following columns:

| chrom | chromStart | chromEnd | name | score | strand | sequence | group_name |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr10 | 100151391 | 100151417 | 1 | 1 | + | GTCTCGAAACCGAGTCTCGTTTCCAAA | gene1 |
| chr10 | 100152125 | 100152143 | 2 | 1 | + | GAACGTAGTTGGACACGAG | gene1 |
| chr10 | 100159361 | 100159382 | 3 | 1 | + | GTTCTTCCGTTTTTAATCTTCT | gene1 |
| chr10 | 100164084 | 100164098 | 4 | 1 | + | AAGAGTCGGAGGGAC | gene1 |
| chr10 | 100183737 | 100183765 | 5 | 1 | + | AGGGGATCCTTAGTGAGTGGACGTGTCTA | gene1 |

Sequence comparisons are performed only between sequences belonging to the same group. The `name` field is used by the `--pairing` option to determine which sequences are compared:

* `--pairing any` (default) — compare all pairs of sequences within each group, including each sequence with itself.
* `--pairing cis` — compare each sequence only with itself.
* `--pairing trans` — compare each sequence with all other sequences within the same group.

### Output

This step outputs the raw structures. The output includes the input interval coordinates, the relative coordinates of the handles within these intervals, and the structure data: energy, alignment, and dot-bracket structure. The file may contain some non-collinear structures or duplicated *cis*-structures.

| row | energy | interval1 | interval2 | start_al1 | end_al1 | start_al2 | end_al2 | alignment1 | alignment2 | structure |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gene1 | -16.1 | chr10_100151391_100151417_+ | chr10_100152125_100152143_+ | 15 | 25 | 8 | 18 | CTCGTTTCCAA | TTGGACACGAG | (((((.((((())))).))))) |
| gene1 | -17.9 | chr10_100159361_100159382_+ | chr10_100164084_100164098_+ | 0 | 13 | 0 | 14 | GTTCTTCCGTTTTT | AAGAGTCGGAGGGAC | (((((((((((((())))).))))))))) |

### Parameters

- `-k <kmer_length>` -  Minimum length of consecutive stacked nucleotide pairs. Must match the value used by `preph-precalculate-energies`. Recommended k = `5`, default `5`.
* `-d <gt_amount_in_kmer_max>` — Maximum number of GT pairs allowed in a k-mer. For k = `5`, recommended: `2`, default `2`.
- `-a <handle_len_min>` - Minimum length of the handles. Recommended = `10`, default `10`.
- `-e <energy_max>` - Maximum energy in kcal/mol. Recommended = `-15`, default `-15`.
- `--need-subopt / --no-need-subopt` - if enabled, attempts to find suboptimal structures. If disabled, returns only MFE structure for each pair of compared intervals. Enabled by default.
- `-j <jobs>` - Number of threads to run in parallel, default 8.

## 2. Normalize the structures

Run `preph-genomic-normalize --input <preph_df> --output <preph_normalized_df>`

### Input

PrePH output with relative coordinates from `preph-genomic-fold`.

### Output

Structures with absolute genomic coordinates. Duplicated and non-collinear structures are removed, and stable IDs are assigned to each structure.

| row | energy | alignment1 | alignment2 | structure | strand | chr | panhandle_start | panhandle_left_hand | panhandle_right_hand | panhandle_end | al1_length | al2_length | id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gene1 | -16.1 | CTCGTTTCCAA | TTGGACACGAG | (((((.((((())))).))))) | + | chr10 | 100151406 | 100151416 | 100152133 | 100152143 | 11 | 11 | 1 |
| gene1 | -17.9 | GTTCTTCCGTTTTT | AAGAGTCGGAGGGAC | (((((((((((((())))).))))))))) | + | chr10 | 100159361 | 100159374 | 100164084 | 100164098 | 14 | 15 | 2 |

## 3. Build genomic BED for visualization

Run `preph-genomic-to-bed --input <preph_normalized_df> --output <preph_bed_df>`

### Input

Structures with absolute genomic coordinates from `preph-genomic-normalize`.

### Output

BED12 for genome browser.

|||||||||||||
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr10 | 100151405 | 100152143 | id=1,dG=-16.1 | 1 | + | 100151405 | 100152143 | 0,100,0 | 2 | 11,11 | 0,727 |
| chr10 | 100159360 | 100164098 | id=2,dG=-17.9 | 1 | + | 100159360 | 100164098 | 0,100,0 | 2 | 14,15 | 0,4723 |

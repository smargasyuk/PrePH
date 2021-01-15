#!/bin/bash
length=$2

# 0-based
awk -F"\t" '{$2 = $2 - 1; print}' OFS="\t" ../data/introns_with_flanks_python.bed | awk -F"\t" '{$3 = $3 - 1; print}' OFS="\t" > ../data/tmp
mv ../data/tmp ../data/introns_with_flanks_python.bed
# conserved
awk -F"\t" '{print $2,$3,$4}' OFS="\t" $1 | bedtools intersect -a ../data/introns_with_flanks_python.bed -b stdin > ../data/conin_python.bed
# sort & merge
awk -F"\t" '{print $1"_"$4"_"$5"_"$6,$2,$3,$4,$5,$6}' OFS="\t" ../data/conin_python.bed | bedtools sort -i stdin | bedtools merge -s -c 4 -o distinct -i stdin | awk -F"\t" '{n = split($1, a, "_"); print a[1],$2,$3,$4,a[4],a[5]}' OFS="\t" > ../data/tmp
# coding genes
bedtools intersect -a ../data/tmp -b ../data/genes.bed -s -wa -wb | awk -F"\t" '{if ($4 == $10) {print $1,$2,$3,$6,$8,$9}}' OFS="\t" > ../data/conin_python_in_genes.bed

echo "intersected with conserved regions, merge, intersected with genes"

# select only long enough
awk -v len="$length" -F"\t" '{if($3-$2 >= len) { print }}' OFS="\t" ../data/conin_python_in_genes.bed > ../data/conin_python_long_in_genes.tsv
echo "selected only long enough"

# add header
echo -e "chrom\tchromStart\tchromEnd\tstrand\tstart_gene\tend_gene" | cat - ../data/conin_python_long_in_genes.tsv > ../data/conin_python_long_final.tsv
echo "added header"

rm ../data/tmp
rm ../data/conin_python_long_in_genes.tsv
rm ../data/conin_python_in_genes.bed
rm ../data/conin_python.bed
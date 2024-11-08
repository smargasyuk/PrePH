#!/bin/bash
length=$2
output_dir=$3

# 0-based
awk -F"\t" '{$2 = $2 - 1; print}' OFS="\t" $output_dir/introns_with_flanks_python.bed | awk -F"\t" '{$3 = $3 - 1; print}' OFS="\t" > $output_dir/tmp
mv $output_dir/tmp $output_dir/introns_with_flanks_python.bed
# conserved
awk -F"\t" '{print $2,$3,$4}' OFS="\t" $1 | bedtools intersect -a $output_dir/introns_with_flanks_python.bed -b stdin > $output_dir/conin_python.bed
# sort & merge
awk -F"\t" '{print $1"_"$4"_"$5"_"$6,$2,$3,$4,$5,$6}' OFS="\t" $output_dir/conin_python.bed | bedtools sort -i stdin | bedtools merge -s -c 4 -o distinct -i stdin | awk -F"\t" '{n = split($1, a, "_"); print a[1],$2,$3,$4,a[4],a[5]}' OFS="\t" > $output_dir/tmp

echo "intersected with conserved regions, merge, intersected with genes"

# select only long enough
awk -v len="$length" -F"\t" '{if($3-$2 >= len) { print }}' OFS="\t" $output_dir/tmp > $output_dir/conin_python_long.tsv
echo "selected only long enough"

rm $output_dir/tmp
rm $output_dir/conin_python.bed
#rm ../data/conin_python_long_in_genes.tsv
#rm ../data/conin_python_in_genes.bed

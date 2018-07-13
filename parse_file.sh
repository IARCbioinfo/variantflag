#!/bin/sh

args=("$@")

#echo Number of arguments: $#

mut_PATH=${args[0]}
tag=${args[1]}
annot_path=${args[2]}
genome=${args[3]}

#grep "#CHROM" ${mut_PATH}raw_VCF_mutect/${tag}_calls.vcf | cut -f 1-5 > ${mut_PATH}header.txt
#cat <(cut -f 1-2,4-5 ${mut_PATH}${tag}_pass_strelka.vcf) <(cut -f 1-2,4-5 ${mut_PATH}${tag}_pass_mutect.vcf) | sort -u -k1,1V -k2,2n > ${mut_PATH}${tag}_allMutations.vcf
cat <(cut -f 1-5 ${mut_PATH}${tag}_pass_strelka.vcf) <(cut -f 1-5 ${mut_PATH}${tag}_pass_mutect.vcf) | sort -u -k1,1V -k2,2n | sed -e '1d' > ${mut_PATH}tmp
cat ${mut_PATH}header.txt ${mut_PATH}tmp > ${mut_PATH}${tag}_allMutations.vcf
sed -e "1d" ${mut_PATH}${tag}_allMutations.vcf | cut -f 1-2  | awk 'BEGIN { OFS = "\t" } {print $1, $2, $2+1}' > ${mut_PATH}${tag}_allPositions.bed 
#bedtools intersect -a ${mut_PATH}${tag}_allPositions.bed -b ${annot_path}${genome}_repeat_masker.bed ${annot_path}${genome}_tandem_repeat.bed -wa -wb -names rpt str > ${mut_PATH}${tag}_allPositions_intersectRepeat.bed
bedtools intersect -a ${mut_PATH}${tag}_allPositions.bed -b ${annot_path}repeat_masker_${genome}.bed ${annot_path}tandem_repeat_${genome}.bed -wa -wb -names rpt str > ${mut_PATH}${tag}_allPositions_intersectRepeat.bed
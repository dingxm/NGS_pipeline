#!/bin/bash

#merge column from different files in different directory
#for i in $1
#do
#  length=${#i}
###	cut -f 10 $i/genes.fpkm_tracking | paste result/gene_expression2.xlsx - | sed "s/FPKM/${i:0:length-7}/" > result/temp
###	mv result/temp result/gene_expression2.xlsx
#  cp $i/accepted_hits.bam bamfiles/${i:0:length-6}_accepted.bam
#done
#for i in bamfiles/*accepted.bam
#do
#  samtools index $i
#done
#bedtools multicov -split -D -bams bamfiles/*accepted.bam -bed mouse_gene.bed > result/ow_gene_total_counts
#generate genotype from vcf file
awk  'BEGIN {ORS="\t"}{print $1"_"$2; for(i=10;i<=57;i++) print substr($i,1,3); printf "\n" }' test.vcf > genotype
#generage input file for fbat
awk  'BEGIN {FS = "/";ORS="\t"};{for (i=1;i<=NF;i++) print $i; printf "\n"}' genotype
#transpose
awk ' BEGIN {ORS="\t"};{for (i=1;i<=NF;i++) {array[NR,i]=$i; if(big <= NF) big=NF}};END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++)  {print array[j,i]}; printf "\n"} }'
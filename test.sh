#!/bin/bash
#samtools mpileup -C50 -DSuf ../genome.fa -l ~/mm10_coding_exon.bed *remove.bam | bcftools view -bvcg - > ../result/var_raw.bcf
genome_fa="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
samtools mpileup -C50 -DSuf $genome_fa -l ~/hg19_coding_re.bed *recal.bam | bcftools view -bvcg - > var_raw.bcf
awk  'BEGIN {ORS="\t"}{print $1"_"$2; for(i=10;i<=57;i++) print substr($i,1,3); printf "\n" }' test.vcf > genotype
 awk ' BEGIN {ORS="\t"};{for (i=1;i<=NF;i++) {array[NR,i]=$i; if(big <= NF) big=NF}};END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++)  {print array[j,i]}; printf "\n"} }' genotype > genotype_transpose

#!/bin/bash
reference=/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa
dbsnp=/u/home/lixuser/data/dbsnp_138.hg19.vcf
gold_indel=/u/home/lixuser/data/Mills_and_1000G_gold_standard.indels.hg19.vcf
sampleNames=$@
mkdir BWA_Alignment
for i in $sampleNames
do 
echo `date`
RG='@RG\tID:'$i'\tSM:'$i'\tLB:'$i'\tPL:illumina'
bwa mem $reference -t 8 -M -R "$RG" $i/${i}_trimmed_R1.fastq $i/${i}_trimmed_R2.fastq | samtools view -bS - > BWA_Alignment/${i}.bam
done
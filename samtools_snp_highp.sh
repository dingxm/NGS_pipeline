#!/bin/bash
human_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human_coding_bed="/u/home/d/dingxm/hg19_coding_merge.bed"
mouse_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
mouse_coding_bed="/u/home/d/dingxm/mm10_coding_exon.bed"
samtools mpileup -C50 -DSuf $human_genome_fa -l $human_coding_bed *.realigned.recal.bam | bcftools view -bvcg - > result/var_raw.bcf
#samtools mpileup -C50 -DSuf $mouse_genome_fa -l $mouse_coding_bed *.bam | bcftools view -bvcg - > var_raw.bcf
bcftools view result/var_raw.bcf | vcfutils.pl varFilter -D 8000 > result/var_filter.vcf
#!/bin/bash
sampleid=$1
ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.* .
human_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
#Next500
tophat -p 8 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
#Hiseq
cd ${sampleid}_thout
java -Xmx8g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar RIBOSOMAL_INTERVALS=~/hg19_ribsome_interval.txt INPUT=accepted_hits.bam REF_FLAT=~/hg19_refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
ngs.plot.r -G hg19 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
#rm accepted_hits.bam
#rm unmapped.bam
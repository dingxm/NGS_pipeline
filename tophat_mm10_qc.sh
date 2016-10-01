#!/bin/bash
sampleid=$1
ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.* .
mouse_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2014-05-23-16-05-10/Genes/genes.gtf"
tophat -p 8 -G $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
cd ${sampleid}_thout
java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=~/mm10_refFlat.txt RIBOSOMAL_INTERVALS=~/mm10_ribsome_interval.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
ngs.plot.r -G mm10 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
#rm accepted_hits.bam
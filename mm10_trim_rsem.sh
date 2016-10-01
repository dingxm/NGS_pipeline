#!/bin/bash
sampleid=$1
ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.* .
mouse_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2014-05-23-16-05-10/Genes/genes.gtf"
#single read
#tophat -p 6 -G $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
#paired read
/u/local/apps/tophat/2.0.9/tophat --no-coverage-search -p 6 -G $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/*_R1.fastq ${sampleid}/*_R2.fastq
cd ${sampleid}_thout
java -Xmx4g -jar /u/home/d/dingxm/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=/u/home/d/dingxm/mm10_refFlat.txt RIBOSOMAL_INTERVALS=/u/home/d/dingxm/mm10_ribsome_interval.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
#ngs.plot.r -G mm10 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
#!/bin/bash

for i in Sample*thout/accepted_hits.bam
do

java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=$i REF_FLAT=~/hg19_refFlat.txt ASSUME_SORTED=True OUTPUT=read_distribution/${i:0:9}_RNAseq_metrics RIBOSOMAL_INTERVALS=~/hg19_ribsome_interval.txt STRAND_SPECIFICITY=NONE

done


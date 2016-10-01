#!/bin/bash
for i in *sorted.bam
do

java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=$i REF_FLAT=~/hg19_refFlat.txt ASSUME_SORTED=True OUTPUT=../read_distri/${i}_RNAseq_metrics STRAND_SPECIFICITY=NONE

done

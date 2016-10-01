#!/bin/bash
for i in Sample*/*fastq.gz; 
do /u/project/lixuser/data/CAP_miRNA/bin/cutadapt -a AGATCGGAAGAGCACACGTCT -b GTTCAGAGTTCTACAGTCCGACGATC -m 17 -f fastq -q 20 -o ${i%gz}trimmed.gz $i;
done
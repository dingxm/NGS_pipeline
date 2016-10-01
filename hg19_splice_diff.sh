#!/bin/bash
human_genes_gtf="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
mouse_genes_gtf="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2014-05-23-16-05-10/Genes/genes.gtf"
sample1=$1
sample2=$2
python ~/MATS.3.0.8/RNASeq-MATS.py -b1 ${sample1}_thout2/accepted_hits.bam -b2 ${sample2}_thout2/accepted_hits.bam -gtf $human_genes_gtf -o MATS_result123 -t paired -len 100 -a 8 -c 0.0001 -analysis U -expressionChange 10000.0
#!/bin/bash
genome=/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
sampleid=$1
/u/home/d/dingxm/bismark_v0.14.5/bismark $genome -p 4 --bam $sampleid/*R1* 

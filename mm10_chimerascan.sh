#!/bin/bash
sampleid=$1
#mkdir ${sampleid}_chimera
#python ~/local/bin/chimerascan_run.py -v -p 8 --rm-tmp /u/home/lixuser/data/mm10_chimera_index ${sampleid}/*R1* ${sampleid}/*R2* ${sampleid}_chimera
python ~/local/bin/chimerascan_run.py -v -p 8 /u/home/lixuser/data/hg19_chimera_index ${sampleid}/*R1* ${sampleid}/*R2* ${sampleid}_chimera
rm ${sampleid}_chimera/tmp/*.fq
rm ${sampleid}_chimera/tmp/*.bam
find ${sampleid}_chimera/tmp ! -name breakpoints.txt  -type f -delete 

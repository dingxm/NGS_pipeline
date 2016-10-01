#!/bin/bash
mkdir rsem_hg19

for i in $@
do

 ~/kallisto_linux-v0.42.4/kallisto quant -i /u/project/lixuser/data/hg19_rsem_index/hg19_transcrpipt.idx -o hg19_rsem/${i}_count -b 100 ${i}/*R1* ${i}/*R2*

done
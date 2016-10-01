#!/bin/bash
ref_path=/u/project/lixuser/data/Data_analysis/LQ2_analysis/Bowtie2_Alignment
for i in Sample*
do 
cd $i
#perl /u/project/lixuser/data/CAP_miRNA/bin/bwa_sam_converter.pl -i ${i}.sam -c -o ${i}.fa -a ${i}.arf;
/u/project/lixuser/data/CAP_miRNA/bin/quantifier.pl -p ${ref_path}/precusor.dna.fa -m ${ref_path}/mature.dna.fa -r ${i}.fa 
cd ..
done
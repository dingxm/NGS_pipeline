#!/bin/bash
target=$1
control=$2
#ratio=$3 
mkdir ${target}_${control}_MACS
cd ${target}_${control}_MACS
macs -t ../Bowtie_Alignment/${target}.sort.bam -c ../Bowtie_Alignment/${control}.sort.bam --name ${target}_CHIP -g mm -p 0.00001 -BS --call-subpeaks --nomodel --shiftsize 85
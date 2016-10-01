#!/bin/bash
sampleNames=$@
#dir_sic="/u/home/d/dingxm/workshop2_data/software/sickle-master"
for i in $sampleNames
do

cd $i
#sickle pe -g -f ${i}*R1* -r ${i}*R2* -t sanger -o ${i}_trimmed_R1.fastq.gz -p ${i}_trimmed_R2.fastq.gz -s ${i}_single.fastq.gz -q 20 -l 35 
#rm *L00*.fastq.gz
#rm ${i}_single.fastq.gz
java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 ${i}*R1* ${i}*R2* ${i}_trimmed_R1.fastq.gz ${i}_trimmed_R1_unpaired.fastq.gz ${i}_trimmed_R2.fastq.gz ${i}_trimmed_R2_unpaired.fastq.gz ILLUMINACLIP:~/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
rm *unpaired*
cd ..
done

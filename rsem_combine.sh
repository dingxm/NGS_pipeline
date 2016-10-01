#!/bin/bash

#merge column from different files in different directory
#filename = rsem_count
#combine TPM
j=1
for i in *genes.results
do
  if [ $j == 1 ]
  then
      cut -f 1 $i > rsem_tpm.txt
      cut -f 1 $i > rsem_count.txt
      cut -f 1 $i > rsem_fpkm.txt
  fi
  length=${#i}
	cut -f 6 $i | paste rsem_tpm.txt - | sed "s/TPM/${i:0:length-14}/" > temp1
  cut -f 5 $i | paste rsem_count.txt - | sed "s/expected_count/${i:0:length-14}/" > temp2
  cut -f 7 $i | paste rsem_fpkm.txt - | sed "s/FPKM/${i:0:length-14}/" > temp3
	mv temp1 rsem_tpm.txt
  mv temp2 rsem_count.txt
  mv temp3 rsem_fpkm.txt
  j=2
  #cp $i/accepted_hits.bam bamfiles/${i:0:length-6}_accepted.bam
done

#combine count
#j=1
#for i in *genes.results
#do
#  if [ $j == 1 ]
#  then
#      
#  fi
#  length=${#i}
#	cut -f 5 $i | paste rsem_count - | sed "s/expected_count/${i:0:length-14}/" > temp2
#	mv temp rsem_count
#  j=2
#  #cp $i/accepted_hits.bam bamfiles/${i:0:length-6}_accepted.bam
#done
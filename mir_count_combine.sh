#!/bin/bash
j=1
for i in Sample*
do
  if [ $j == 1 ]
  then
      cut -f 1 $i/*csv > mir_count.txt
  fi
  
	cut -f 5 $i/*csv | paste mir_count.txt - | sed "s/seq/${i}/" > temp1
 	mv temp1 mir_count.txt
  j=2
  #cp $i/accepted_hits.bam bamfiles/${i:0:length-6}_accepted.bam
done
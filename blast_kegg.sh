#!/bin/bash
for i in *genes.pep
do
   blastp -db /u/project/lixuser/data/athdb/ath.pep.fasta -query $i -outfmt 6 -out ${i}.blast.tab
done 
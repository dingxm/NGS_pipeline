#!/bin/bash

##########################
## Usage
##########################
usage() {
  echo ""
  echo "   Usage: `basename $0` reads_length"
  echo "   Creates summary of reads distribution in RNAseq"
  echo ""
  echo "   example: `basename $0` 50"
  echo ""
  echo "   The expected input is: arg1=the length of reads"
  exit 2
}

## is the argument supplied?
   if [ "$#" != "1" ]
   then
      usage
      exit 1
   fi
header="Sample\tTotal_Reads\tMapped_Reads\tRibsome_Reads\tCoding_Reads\tUTR_Reads\tIntronic_Reads\tIntergenic_Reads\tPCT_Ribsome_Reads\tPCT_Coding_Reads\tPCT_UTR_Reads\tPCT_Intronic_Reads\tPCT_Intergenic_reads\tPCT_mRNA_Reads"

echo -e $header > summary_of_alignment.xls
read_length=$1
for i in *thout
do
length=${#i}
file1=${i}/align_summary.txt
#cat ${i}/align_summary.txt
sample=${i:0:length-6}
total_reads=$(grep "Input" $file1 | cut -d ":" -f2 | sed -e 's/^[\t]*//' | awk '{if (NR==1) print $1}')
mapped_reads=$(grep "Mapped" $file1 | cut -d ":" -f2 | sed -e 's/^[\t]*//' | awk '{if (NR==1)print $1}')
file2=${i}/${sample}_RNAseq_metrics
ribsome_reads=$(sed -n '8p' $file2 | awk -v len=$read_length '{print int($3/len+0.5)}')
coding_reads=$(sed -n '8p' $file2 | awk -v len=$read_length '{print int($4/len+0.5)}')
utr_reads=$(sed -n '8p' $file2 | awk -v len=$read_length  '{print int($5/len+0.5)}')
intronic_reads=$(sed -n '8p' $file2 | awk -v len=$read_length '{print int($6/len+0.5)}')
intergenic_reads=$(sed -n '8p' $file2 | awk -v len=$read_length '{print int($7/len+0.5)}')
ribsome_pct=$(sed -n '8p' $file2 | awk '{print $11}')
coding_pct=$(sed -n '8p' $file2 | awk '{print $12}')
utr_pct=$(sed -n '8p' $file2 | awk '{print $13}')
intronic_pct=$(sed -n '8p' $file2 | awk '{print $14}')
intergenic_pct=$(sed -n '8p' $file2 | awk '{print $15}')
mrna_pct=$(sed -n '8p' $file2 | awk '{print $16}')
echo -e "$sample\t$total_reads\t$mapped_reads\t$ribsome_reads\t$coding_reads\t$utr_reads\t$intronic_reads\t$intergenic_reads\t$ribsome_pct\t$coding_pct\t$utr_pct\t$intronic_pct\t$intergenic_pct\t$mrna_pct" >> summary_of_alignment.xls
done




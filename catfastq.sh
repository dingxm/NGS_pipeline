#!/bin/bash
##########################
## Usage
##########################
usage() {
  echo ""
  echo "   Usage: `basename $0` PI_name Dir_name"
  echo "   combine reads together"
  echo ""
  echo "   example: `basename $0` XD /u/project/lixuser/data/NextSeq500/151020_NS500551_0105_AHHVV2BGXX "
  echo ""
  exit 2
}

## is the argument supplied?
   if [ "$#" != "3" ]
   then
      usage
      exit 1
   fi
pi=$1
dir=$2
reads=$3
length=${#dir}
mkdir /u/nobackup/lixuser/data/Data_analysis/${pi}_analysis
for i in $dir/*
do
   file_name=`basename "$i"`
#  index=`expr index "$file_name" _S`
#  sample=${file_name:0:($index-1)}#get sample name
   sample=${file_name%_S*gz}
	 mkdir /u/nobackup/lixuser/data/Data_analysis/${pi}_analysis/${sample}  
   #cp $dir/$sample* /u/project/lixuser/data/Data_analysis/${pi}_analysis/${sample}
done
cd /u/nobackup/lixuser/data/Data_analysis/${pi}_analysis
for i in *
do
  if [ $reads == 2 ];  then
    cat $dir/${i}_*R1*.fastq.gz >> $i/${i}_R1.fastq.gz
    cat $dir/${i}_*R2*.fastq.gz >> $i/${i}_R2.fastq.gz
  else 
    cat $dir/${i}_*R1*.fastq.gz >> $i/${i}_R1.fastq.gz
  fi  
done
#for i in $dir/*
#do
#   file_name=${i:length+1}
##  index=`expr index "$file_name" _S`
##  sample=${file_name:0:($index-1)}#get sample name
#  sample=${file_name%_S[1-9]*gz}
#	mkdir /u/project/lixuser/data/Data_analysis/${pi}_analysis/${sample}  
#  cp $dir/${sample}_* /u/project/lixuser/data/Data_analysis/${pi}_analysis/${sample}
#done
#for i in $dir/Sample*
#do
#  sample=${i:length+1}
#  if [ $reads == 2 ]; then
#	   cat $i/*1*.gz >> /u/project/lixuser/data/Data_analysis/${pi}_analysis/$sample/${sample}_R1.fastq.gz
#	   cat $i/*2*.gz >> /u/project/lixuser/data/Data_analysis/${pi}_analysis/$sample/${sample}_R2.fastq.gz
#  else
#      cat $i/*R1*.fastq.gz > /u/project/lixuser/data/Data_analysis/${pi}_analysis/$sample/${sample}_R1.fastq.gz
#  fi
#  #zcat $i/*.gz > /u/home/lixuser/data/${pi}_analysis/$sample/${sample}_R1.fastq
#done
cd /u/nobackup/lixuser/data/Data_analysis/${pi}_analysis
mkdir fastqc
for i in */*.fastq*
do
    ~/FastQC/fastqc -o fastqc $i
done


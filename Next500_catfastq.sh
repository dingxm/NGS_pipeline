#!/bin/bash
##########################
## Usage
##########################
usage() {
  echo ""
  echo "   Usage: `basename $0` PI_name Dir_name reads_type"
  echo "   combine reads together"
  echo ""
  echo "   example: `basename $0` XD /u/project/lixuser/data/NextSeq500/151020_NS500551_0105_AHHVV2BGXX 2"
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
#sampleNames=$@
#dir=$1
#cd $dir
for i in $dir/*
do
   file_name=`basename "$i"`
#  index=`expr index "$file_name" _S`
#  sample=${file_name:0:($index-1)}#get sample name
   sample=${file_name%_S[1-9]*gz}
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
#cd -
#for i in *
#do 
#  qsub -cwd -V -N $i -l h_data=4G,h_rt=24:00:00,arch=intel* -m bea -pe shared 8 ~/scripts/hg19_germline_SNV.sh $i 
#done 
mkdir report_of_raw_fastq_quality
for i in */*.fastq*
do
    ~/FastQC/fastqc -o report_of_raw_fastq_quality $i
done


#!/bin/bash
sampleid=$1
#mkdir rsem_newt
#rsem-calculate-expression --bowtie2 -p 4 --no-bam-output $sampleid/${sampleid}_*.fastq.gz  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/Newt_rsem_index/newt_ref_rsem rsem_newt/$sampleid
#mkdir rsem_pa
#rsem-calculate-expression --bowtie2 -p 4 --no-bam-output $sampleid/*fastq*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/Pseudomonas_aeruginosa_rsem_index/Pseudomonas_aeruginosa_ref_rsem rsem_pa/$sampleid
#mkdir rsem_ss3
#rsem-calculate-expression --bowtie2 -p 4 --no-bam-output $sampleid/*fastq*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/ss3en_rsem_index/ss3_ref_rsem rsem_ss3/$sampleid
cd $sampleid
seqtk trimfq -b 3 *R1* | gzip > ${sampleid}_R1_trim.fastq.gz
rm *R1.fastq.gz
cd ..
mkdir rsem_bat
rsem-calculate-expression --bowtie2 -p 4 --no-bam-output $sampleid/*fastq*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/little_brown_bat/rsem_bat_index/bat rsem_bat/$sampleid
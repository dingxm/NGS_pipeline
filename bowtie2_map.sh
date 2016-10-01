#!/bin/bash
reference=/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Pabies/pabties1
#mkdir Bowtie2_Alignment
for i in Sample_7
do
  sampleId=$i
#reference=$2
#bowtie2 -x bowtie2_index/$reference -1 $sampleId/${sampleId}_R1.fastq -2 $sampleId/${sampleId}_R2.fastq -p 4 | samtools view -bS - > Bowtie2_Alignment/${sampleId}.bam
  /u/home/d/dingxm/bowtie2-2.2.5/bowtie2 -x $reference -U $sampleId/*trimmed* -p 4 > $sampleId/${sampleId}.sam
  samtools view -bS $sampleId/${sampleId}.sam > Bowtie2_Alignment/${sampleId}.bam
  cd Bowtie2_Alignment
  samtools sort ${sampleId}.bam ${sampleId}.sorted
  cd ..
done

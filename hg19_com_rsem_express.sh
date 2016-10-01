#!/bin/bash
##########################
## Usage
##########################
usage() {
  echo ""
  echo "   Usage: `basename $0` Sample_name reads_type sample_number reads_length"
  echo "   combine reads together"
  echo ""
  echo "   example: `basename $0` XD MCF_T1 2 10 50"
  echo ""
  exit 2
}

## is the argument supplied?
   if [ "$#" != "4" ]
   then
      usage
      exit 1
   fi
sampleid=$1
reads=$2
sample_number=$3
length=$4
mkdir report_of_quality_of_raw_data
cd $sampleid
java -jar /u/home/d/dingxm/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 ${sampleid}*R1*.gz ${sampleid}*R2*.gz ${sampleid}_trimmed_R1.fastq.gz ${sampleid}_trimmed_R1_unpaired.fastq.gz ${sampleid}_trimmed_R2.fastq.gz ${sampleid}_trimmed_R2_unpaired.fastq.gz ILLUMINACLIP:/home/dingxm/NGS_software/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
#java -jar ~/NGS_software/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 -threads 8 ${sampleid}*R1*.gz ${sampleid}_trimmed_R1.fastq.gz ILLUMINACLIP:/home/dingxm/NGS_software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
~/NGS_software/FastQC/fastqc -o ../report_of_quality_of_raw_data *trimmed*
rm *unpaired*
rm ${sampleid}_R*
cd ..
mkdir rsem_hg19
mkdir rsem_hg19lincRNA
ln -s /home/dingxm/genome_reference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.* .
human_genes_gtf="/home/dingxm/genome_reference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
if [ $reads == 2 ]; then
  rsem-calculate-expression --bowtie2 -p 4 --bowtie2 --no-bam-output --paired-end $sampleid/*trimmed_R1* $sampleid/*trimmed_R2* /home/dingxm/genome_reference/hg19_rsem_index/hg19_ref_rsem rsem_hg19/${sampleid}
  #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --paired-end $sampleid/*R1* $sampleid/*R2* /u/project/lixuser/data/mm10lincRNA_rsem_index/mm10lincRNA_ref_rsem rsem_mm10lincRNA/$sampleid
  tophat -p 4 -G  $human_genes_gtf -g 1 -o ${sampleid}_thout genome  ${sampleid}/*trimmed_R1* ${sampleid}/*trimmed_R2*
else
    rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --fragment-length-mean 170.0 --fragment-length-sd 60.0 $sampleid/*trimmed_R1*  /home/dingxm/genome_reference/hg19_rsem_index/hg19_ref_rsem rsem_hg19/$sampleid
    rsem-calculate-expression --bowtie2 -p 4 --no-bam-output $sampleid/*trimmed_R1*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /home/dingxm/genome_reference/hg19lincRNA_rsem_index/hg19lincRNA_ref_rsem rsem_hg19lincRNA/$sampleid
    tophat -p 4 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome $sampleid/*trimmed_R1*
fi
cd ${sampleid}_thout
java -Xmx8g -jar /home/dingxm/NGS_software/picard-tools-1.119/CollectRnaSeqMetrics.jar RIBOSOMAL_INTERVALS=/home/dingxm/genome_reference/hg19_ribsome_interval.txt INPUT=accepted_hits.bam REF_FLAT=/home/dingxm/genome_reference/hg19_refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
ngs.plot.r -G hg19 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
#rm accepted_hits.bam
rm unmapped.bam
cd ..
number=`ls -dl *_thout/*heatmap.pdf | wc -l`
if [ $number == $sample_number ]
then
  mkdir report_of_genebody_coverage
  mkdir report_of_alignment_statistics
  ~/scripts/summary_rna_seq_alignment.sh $length
  cp *thout/*pdf report_of_genebody_coverage
  cp summary_of_alignment.xls report_of_alignment_statistics
  cd rsem_*
  ~/scripts/rsem_combine.sh
fi
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
mkdir rsem_hg19
mkdir rsem_hg19lincRNA
ln -s /u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.* .
human_genes_gtf="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
#cd $sampleid
#seqtk trimfq -b 3 *R1* | gzip > ${sampleid}_R1_trim.fastq.gz
#rm *R1.fastq.gz
#cd ..
if [ $reads == 2 ]; then
  rsem-calculate-expression --bowtie2 -p 8 --bowtie2 --no-bam-output --paired-end $sampleid/*R1* $sampleid/*R2* /u/nobackup/lixuser/data/mm10_rsem_index/mm10_ref_rsem rsem_mm10/${sampleid}
  #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --paired-end $sampleid/*R1* $sampleid/*R2* /u/nobackup/lixuser/data/mm10lincRNA_rsem_index/mm10lincRNA_ref_rsem rsem_mm10lincRNA/$sampleid
  tophat -p 8 -G  $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome  ${sampleid}/*R1* ${sampleid}/*R2*
else
    #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --fragment-length-mean 170.0 --fragment-length-sd 60.0 $sampleid/*R1*  /u/nobackup/lixuser/data/hg19_rsem_index/hg19_ref_rsem rsem_hg19/$sampleid
    #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output ${sampleid}/*R1*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/nobackup/lixuser/data/mm10lincRNA_rsem_index/mm10lincRNA_ref_rsem rsem_mm10lincRNA/$sampleid
    tophat -p 4 -G $human_genes_gtf -g 10 -o ${sampleid}_thout genome ${sampleid}/*R2*
fi
cd ${sampleid}_thout
java -Xmx4g -jar /u/home/d/dingxm/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=/u/home/d/dingxm/hg19_refFlat.txt RIBOSOMAL_INTERVALS=/u/home/d/dingxm/hg19_ribsome_interval.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
ngs.plot.r -G hg19 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
#mkdir ../bamfiles
#samtools index accepted_hits.bam
#mv accepted_hits.bam ../bamfiles/${sampleid}.bam
#mv accepted_hits.bam.bai ../bamfiles/${sampleid}.bam.bai
rm accepted_hits.bam
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
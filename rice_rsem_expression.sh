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
#rsem-prepare-reference --gtf genes.gtf --bowtie2 --bowtie2-path /u/local/apps/bowtie2/2.1.0 genome.fa /u/home/lixuser/data/hg19_ref_rsem
#rsem-prepare-reference --gtf $genes_gtf --bowtie2 --bowtie2-path /u/local/apps/bowtie2/2.1.0 $genome_fa /u/home/lixuser/data/cmcp6_rsem_index/cmcp6_ref_rsem
#mkdir rsem_cmcp6
#rsem-calculate-expression --bowtie2 -p 8 --no-bam-output --paired-end $sampleid/*R1.fastq $sampleid/*R2.fastq /u/home/lixuser/data/hg19_rsem_index/hg19_ref_rsem rsem_hg19/$sampleid

#rsem-calculate-expression --bowtie2 -p 8 --no-bam-output $sampleid/${sampleid}_*.fastq  /u/home/lixuser/data/cmcp6_rsem_index/cmcp6_ref_rsem rsem_cmcp6/$sampleid

mkdir rsem_rice

ln -s /u/project/lixuser/data/HiSeq2000_data/GenomeReference/Oryza_sativa_japonica/Ensembl/MSU6/Sequence/Bowtie2Index/genome.* .
rice_genes_gtf="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Oryza_sativa_japonica/Ensembl/MSU6/Annotation/Archives/archive-2015-07-17-14-34-10/Genes/genes2.gtf"
ref_flat=/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Oryza_sativa_japonica/Ensembl/MSU6/Annotation/Archives/archive-2015-07-17-14-34-10/Genes/refFlat.txt
#mkdir rsem_hg19lincRNA
#Next500 
if [ $reads == 2 ]; then
    
    rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --paired-end $sampleid/*R1* $sampleid/*R2* /u/project/lixuser/data/hg19_rsem_index/hg19_ref_rsem rsem_hg19/$sampleid
    #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --paired-end $sampleid/*R1.fastq* $sampleid/*R2.fastq* /u/project/lixuser/data/hg19lincRNA_rsem_index/hg19lincRNA_ref_rsem rsem_hg19lincRNA/$sampleid
    tophat -p 4 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/*R1* ${sampleid}/*R2*
else

    rsem-calculate-expression --bowtie2 -p 4 --no-bam-output ${sampleid}/*R1*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/rice_rsem_index/rice_ref_rsem rsem_rice/$sampleid
    #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output ${sampleid}/*R1*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/hg19lincRNA_rsem_index/hg19lincRNA_ref_rsem rsem_hg19lincRNA/$sampleid
    #tophat -p 4 -G $rice_genes_gtf -g 1 -o ${sampleid}_thout genome  ${sampleid}/${sampleid}*R1*
fi
# alignment QC
cd ${sampleid}_thout
#java -Xmx8g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=$ref_flat ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
#java -Xmx8g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=${sampleid}.dup.txt INPUT=accepted_hits.bam ASSUME_SORTED=True OUTPUT=${sampleid}.remove.bam
#ngs.plot.r -G hg19 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
rm accepted_hits.bam
rm unmapped.bam
#mkdir ../bamfiles
#samtools index accepted_hits.bam
#mv accepted_hits.bam ../bamfiles/${sampleid}.bam
#mv accepted_hits.bam.bai ../bamfiles/${sampleid}.bam.bai
#cd ..
#number=`ls -dl *_thout//*heatmap.pdf | wc -l`
#if [ $number == $sample_number ]
#then
#  mkdir report_of_genebody_coverage
#  mkdir report_of_alignment_statistics
#  ~/scripts/summary_rna_seq_alignment.sh $length
#  cp *thout/*pdf report_of_genebody_coverage
#  cp summary_of_alignment.xls report_of_alignment_statistics
#  cd rsem_*
#  ~/scripts/rsem_combine.sh
#fi





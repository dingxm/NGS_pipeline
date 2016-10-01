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

mkdir rsem_rn5
ln -s /u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Sequence/Bowtie2Index/genome.* .
rat_genes_gtf="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Annotation/Archives/archive-2013-03-06-19-22-38/Genes/genes.gtf"
if [ $reads == 2 ]; then
  #rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --paired-end $sampleid/*R1* $sampleid/*R2* /u/project/lixuser/data/rn5_rsem_index/rn5_ref_rsem rsem_rn5/${sampleid}
  tophat -p 4 -G  $rat_genes_gtf -g 1 -o ${sampleid}_thout2 genome  ${sampleid}/*R1* ${sampleid}/*R2*
else
    rsem-calculate-expression --bowtie2 -p 4 --no-bam-output --fragment-length-mean 170.0 --fragment-length-sd 60.0 $sampleid/*R1*  /u/nobackup/lixuser/data/rn5_rsem_index/rn5_ref_rsem rsem_rn5/$sampleid
    tophat -p 4 -G $rat_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/*R1*
fi
cd ${sampleid}_thout
java -Xmx4g -jar /u/home/d/dingxm/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam RIBOSOMAL_INTERVALS=~/rn5_ribsome_interval.txt REF_FLAT=/u/home/d/dingxm/rat_refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
ngs.plot.r -G rn5 -R genebody -C accepted_hits.bam -O $sampleid -T $sampleid -F rnaseq
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

#ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Sequence/Bowtie2Index/genome.* .
#rat_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Annotation/Archives/archive-2014-05-23-16-07-35/Genes/genes.gtf"
#tophat -p 8 -G $rat_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
#cd ${sampleid}_thout
#java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=~/rat_refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
#rm accepted_hits.bamhaet
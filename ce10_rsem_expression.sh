#!/bin/bash
sampleid=$1
mkdir rsem_ce10

#rsem-calculate-expression --bowtie2 -p 8 --no-bam-output ${sampleid}_*.fastq  /u/home/lixuser/data/mm10_rsem_index/mm10_ref_rsem rsem_mm10/$sampleid
#gene expression estimation
#single read
#rsem-calculate-expression --bowtie2 -p 6 --no-bam-output --fragment-length-mean 170.0 --fragment-length-sd 60.0 $sampleid/${sampleid}_R1.fastq  /u/home/lixuser/data/mm10_rsem_index/mm10_ref_rsem rsem_mm10/$sampleid
#paired read
/u/home/d/dingxm/rsem-1.2.15/rsem-calculate-expression --bowtie2 -p 6 --no-bam-output --forward-prob 0 --paired-end $sampleid/*_R1.fastq $sampleid/*_R2.fastq /u/home/lixuser/data/ce10_rsem_index/ce10_ref_rsem rsem_ce10/${sampleid}_2

#QC
ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index/genome.* .
celegans_genes_gtf=/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Caenorhabditis_elegans/UCSC/ce10/Annotation/Archives/archive-2014-05-23-16-02-13/Genes/genes.gtf
#single read
#tophat -p 6 -G $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
#paired read
/u/local/apps/tophat/2.0.9/tophat --no-coverage-search -p 6 -G  --fr-firststrand $celegans_genes_gtf -g 1 -o ${sampleid}_thout1 genome ${sampleid}/*_R1.fastq ${sampleid}/*_R2.fastq
cd ${sampleid}_thout1
java -Xmx4g -jar /u/home/d/dingxm/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Caenorhabditis_elegans/UCSC/ce10/Annotation/Archives/archive-2014-05-23-16-02-13/Genes/refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
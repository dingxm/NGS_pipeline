#!/bin/bash
/u/home/d/dingxm/bbmap/bbnorm.sh in=all.fastq.gz out=normalized.fastq.gz ecc
/u/home/d/dingxm/trinityrnaseq-2.0.6/Trinity --seqType fq --max_memory 120G --single normalized.fastq.gz --output trinity_output --CPU 4 
cd trinity_output
u/home/d/dingxm/trinityrnaseq-2.0.6/util/TrinityStats.pl Trinity.fasta > assembly_statistics.txt
/u/home/d/dingxm/trinityrnaseq-2.0.6/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > gene_to_trans_map
rsem-prepare-reference Trinity.fasta --transcript-to-gene-map gene_to_trans_map --bowtie2 --bowtie2-path /u/local/apps/bowtie2/2.1.0 /u/project/lixuser/data/Data_analysis/LQ_analysis/assembly/rsem_index/ref_rsem
sampleid=$1
mkdir rsem
sample_name=${sampleid/-treated-/-}
rsem-calculate-expression --bowtie2 -p 4 --no-bam-output ${sampleid}/${sampleid}*R1*  --fragment-length-mean 180.0 --fragment-length-sd 70.0 /u/project/lixuser/data/Data_analysis/LQ_analysis/assembly/rsem_index/ref_rsem rsem/$sample_name

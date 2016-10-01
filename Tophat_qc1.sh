#!/bin/bash
sampleid=$1
#dir_sic="/u/home/d/dingxm/workshop2_data/software/sickle-master"
ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.* .
#ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Sequence/Bowtie2Index/genome.* .
#ln -s /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.* .
human_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
#mouse_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2014-05-23-16-05-10/Genes/genes.gtf"
#rat_genes_gtf="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Rattus_norvegicus/UCSC/rn5/Annotation/Archives/archive-2014-05-23-16-07-35/Genes/genes.gtf"
#$dir_sic/sickle se -f ${sampleid}_R1.fastq -t sanger -o ${sampleid}_trimmed_R1.fastq -q 20 -l 25 
tophat -p 8 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_R1.fastq
cd ${sampleid}_thout
java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar RIBOSOMAL_INTERVALS=~/hg19_ribsome_interval.txt INPUT=accepted_hits.bam REF_FLAT=~/hg19_refFlat.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
rm accepted_hits.bam
#tophat -p 8 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/${sampleid}_*.fastq
#java -Xmx4g -jar ~/picard-tools-1.115/CollectRnaSeqMetrics.jar INPUT=accepted_hits.bam REF_FLAT=~/mm10_refFlat.txt RIBOSOMAL_INTERVALS=~/mm10_ribsome_interval.txt ASSUME_SORTED=True OUTPUT=${sampleid}_RNAseq_metrics STRAND_SPECIFICITY=NONE
#tophat -p 8 -G $human_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/*R1.fastq ${sampleid}/*R2.fastq
#tophat -p 8 -G $mouse_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}/*R1.fastq ${sampleid}/*R2.fastq
#tophat -p 8 -G $rat_genes_gtf -g 1 -o ${sampleid}_thout genome ${sampleid}_*.fastq
#samtools view -f 2 -bh ${sampleid}_thout/accepted_hits.bam > bamfiles/${sampleid}.bam
#java -Xmx4g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=bamfiles/${sampleid}_dup.txt INPUT=${sampleid}_thout/accepted_hits.bam ASSUME_SORTED=True OUTPUT=bamfiles/${sampleid}_remove.bam
#mv ${sampleid}_thout/accepted_hits.bam bamfiles/${sampleid}_accepted.bam
#cufflinks -p 8 -o ${sampleid}_clout2 -G genes.gtf -M mm_rRNA_tRNA.gtf ${sampleid}_thout/accepted_hits.bam
#bedtools -bams bamfiles/*.bam -bed mouse_gene.bed > ow_gene_counts
#cuffdiff -o AKT_Nor -b genome.fa -p 8 -u genes.gtf -L AKT,prostate bamfiles/Sample_062813NICD-AKT_accepted.bam,bamfiles/Sample_082413NICD-AKT_accepted.bam  bamfiles/Sample_MouseProstate_accepted.bam,bamfiles/Sample_MouseProstate2_accepted.bam
#cuffdiff -o KRas_Nor -b genome.fa -p 8 -u genes.gtf -L KRas,prostate bamfiles/Sample_082213NICD-kRas_accepted.bam,bamfiles/Sample_082413NICD-kRas_accepted.bam  bamfiles/Sample_MouseProstate_accepted.bam,bamfiles/Sample_MouseProstate2_accepted.bam
#cuffdiff -o MYC_Nor -b genome.fa -p 8 -u genes.gtf -L MYC,prostate bamfiles/Sample_082313NICD-Myc_accepted.bam,bamfiles/Sample_110113NICD-Myc_accepted.bam  bamfiles/Sample_MouseProstate_accepted.bam,bamfiles/Sample_MouseProstate2_accepted.bam
#python ~/local/bin/chimerascan_run.py -v --quals solexa /u/home/lixuser/data/hg19_chimera_index ~/chimerascan-0.4.5/tests/vcap_pe_53bp/TMPRSS2-ERG_1.fq ~/chimerascan-0.4.5/tests/vcap_pe_53bp/TMPRSS2-ERG_2.fq /u/home/lixuser/data/test1
#rsem-calculate-expression --bowtie2 -p 8 --no-bam-output --paired-end $sampleid/*R1.fastq $sampleid/*R2.fastq /u/home/lixuser/data/mm10_rsem_index/mm10_ref_rsem RSEM_result/$sampleid
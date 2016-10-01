#!/bin/bash
sampleid=$1
sample_number=$2
dir=`pwd`
human_genes_gtf="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
genome_path="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/STARindex/hg19"
human_repeats="/u/nobackup/lixuser/data/repeats/hg19_repeats.gtf"
human_anno="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
human_genome_fa="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#mkdir report_of_quality_of_raw_data
#cd $sampleid
#java -jar /u/home/d/dingxm/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 *R1* *R2* ${sampleid}_trimmed_R1.fastq.gz ${sampleid}_trimmed_R1_unpaired.fastq.gz ${sampleid}_trimmed_R2.fastq.gz ${sampleid}_trimmed_R2_unpaired.fastq.gz ILLUMINACLIP:/u/home/d/dingxm/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
##java -jar ~/NGS_software/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 -threads 8 ${sampleid}*R1*.gz ${sampleid}_trimmed_R1.fastq.gz ILLUMINACLIP:/home/dingxm/NGS_software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#rm *unpaired*
#rm ${sampleid}_R*
#~/FastQC/fastqc -o ../report_of_quality_of_raw_data *trimmed*
#cd ..
#mkdir ${sampleid}_stout
#STAR --runThreadN 8   --genomeDir $genome_path  --outSAMtype BAM SortedByCoordinate --readFilesIn $sampleid/*trimmed_R1* $sampleid/*trimmed_R2*   --readFilesCommand zcat  --outFileNamePrefix ${sampleid}_stout/$sampleid --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15
#
#STAR --runThreadN 8   --genomeDir $genome_path  --outSAMtype None --readFilesIn $sampleid/*trimmed_R1*  --readFilesCommand zcat   --outFileNamePrefix ${sampleid}_stout/${sampleid}_1 --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15
#
#STAR --runThreadN 8   --genomeDir $genome_path  --outSAMtype None --readFilesIn $sampleid/*trimmed_R2*  --readFilesCommand zcat   --outFileNamePrefix ${sampleid}_stout/${sampleid}_2 --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15
cd ${sampleid}_stout
samtools index *.bam
mkdir circ_result
cd circ_result
echo -e "${dir}/${sampleid}_stout/${sampleid}Chimeric.out.junction" >> samplesheet
echo -e "${dir}/${sampleid}_stout/${sampleid}_1Chimeric.out.junction" >> mate1
echo -e "${dir}/${sampleid}_stout/${sampleid}_2Chimeric.out.junction" >> mate2
number=`wc -l mate2 | cut -f1 -d " "`
if [ $number == $sample_number ]
then
DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -D -R $human_repeats -an $human_anno -Pi -F -M -Nr 2 1 -fg -A $human_genome_fa
fi

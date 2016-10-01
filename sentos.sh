#!/bin/bash
sample=$1
group=$2
fastq1=$3
fastq2=$4

#Set up environment
export SENTIEON_LICENSE=lm:8980
PATH=$PATH:/u/local/apps/sentieon-genomics/201603/bin/
export PATH

ref="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
dbsnpMills="/u/nobackup/lixuser/data/Mills_and_1000G_gold_standard.indels.hg19.vcf"


echo "Start Time:"; date


#Align reads
#bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:ILLUMINA" -t 8 $ref $fastq1 $fastq2 | sentieon util sort -o $sample.sorted.bam -t 8 --sam2bam -i -
#echo "BWA End Time:"; date

#Data metrics
sentieon driver -t 4 -r $ref -i $sample.sorted.bam --algo GCBias --summary $sample.GCsummary.txt $sample.GCmetrics.txt --algo MeanQualityByCycle $sample.MeanQualMetrics.txt --algo QualDistribution $sample.QualDistribMetrics.txt --algo InsertSizeMetricAlgo $sample.InsertSizeMetrics.txt --algo AlignmentStat $sample.AlignMetrics.txt
#sentieon plot metrics -o $sample.Metrics.pdf gc=$sample.GCmetrics.txt mq=$sample.MeanQualMetrics.txt qd=$sample.QualDistribMetrics.txt isize=$sample.InsertSizeMetrics.txt
echo "Data Metrics End Time:"; date

#Remove duplicates
sentieon driver -t 8 -i $sample.sorted.bam --algo LocusCollector --fun score_info $sample.score.txt
sentieon driver -t 8 -i $sample.sorted.bam --algo Dedup --rmdup --score_info $sample.score.txt --metrics $sample.DedupMetrics.txt $sample.dedup.bam
echo "Dedup End Time:"; date

#Indel realignment
sentieon driver -t 8 -r $ref -i $sample.dedup.bam --algo Realigner -k $dbsnpMills $sample.realigned.bam
echo "Indel Realign End Time:"; date


#Base Quality Score Recalibration (BQSR)
sentieon driver -t 8 -r $ref -i $sample.realigned.bam --algo QualCal -k $dbsnpMills $sample.RecalData.table
sentieon driver -t 8 -r $ref -i $sample.realigned.bam -q $sample.RecalData.table --algo QualCal -k $dbsnpMills $sample.RecalData.table.post --algo ReadWriter $sample.recalibrated.bam
sentieon driver -t 8 --algo QualCal --plot --before $sample.RecalData.table --after $sample.RecalData.table.post $sample.RecalResults.csv
#sentieon plot bqsr -o $sample.Recal.pdf $sample.RecalResults.csv
echo "BQSR End Time:"; date

#Variant calling
#sentieon driver -t 4 -r $ref -i $sample.realigned.bam -q $sample.RecalData.table --algo Haplotyper $sample.vcf
#echo "Variant Calling End Time:"; date

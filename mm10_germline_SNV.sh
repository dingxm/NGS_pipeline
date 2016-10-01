#!/bin/bash
reference=/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.6.0/genome.fa
mouse_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
mouse_coding_bed="/u/home/d/dingxm/mm10_coding_exon_merge.bed"
mkdir BWA_Alignment
i=$1
sample_number=1
echo `date`
RG='@RG\tID:'$i'\tSM:'$i'\tLB:'$i'\tPL:illumina'
#bwa mem $reference -t 4 -M -R "$RG" $i/${i}*R1* $i/${i}*R2* | samtools view -bS - > BWA_Alignment/${i}.bam
#describe=`echo ${i} | sed 's/.bam//'`
samtools sort BWA_Alignment/${i}.bam BWA_Alignment/${i}.sorted
java -Xmx8g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=BWA_Alignment/${i}.dup.txt INPUT=BWA_Alignment/${i}.sorted.bam ASSUME_SORTED=True OUTPUT=BWA_Alignment/${i}.remove.bam
rm BWA_Alignment/${i}.bam
rm BWA_Alignment/${i}.sorted.bam
samtools index BWA_Alignment/${i}.remove.bam
echo `date`
cd BWA_Alignment
samtools mpileup -C50 -DSuf $mouse_genome_fa -l $mouse_coding_bed ${i}.remove.bam | bcftools view -bvcg - > ${i}_var_raw.bcf
bcftools view ${i}_var_raw.bcf | vcfutils.pl varFilter -D 1000 > ${i}_var_filter.vcf
#find *bam | parallel 'bedtools coverage -hist -abam {} -b ~/mm10_coding_exon_merge.bed | grep ^all > {}.hist.all.txt'
#Rscript ~/mouse_exome_coverage.R
#SNV annotation
convert2annovar.pl -format vcf4 ${i}_var_filter.vcf -outfile ${i}.avinput
table_annovar.pl ${i}.avinput /u/project/lixuser/data/annovar/mousedb -buildver mm10 -out ${i} -remove -protocol refGene,snp138 -operation g,f -nastring na
echo -e "Chr\tStart\tEnd\tRef\tEnd\tType\tQuality\tDepth" > aviheader
echo "ref-forward:ref-reverse:alt-forward:alt-reverse" > header
cat aviheader ${i}.avinput > ${i}_head
cut -f 7-8 ${i}_head | paste ${i}.mm10_multianno.txt - > ${i}.SNV_dep2.txt
cat ${i}_var_filter.vcf | grep -oP '(?<=DP4=)[^;]*' | sed 's/,/:/g' | cat header - | paste ${i}.SNV_dep2.txt - > ${i}.SNV_dep.txt
#rm ${i}.SNV_dep2.txt
##filter synonymous
#for i in *dep.txt; do cat $i | grep -v rs[0-9]* | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" )  && ($16=="Depth" || $16 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV.txt; done
cat ${i}.SNV_dep.txt |  awk 'BEGIN{FS=OFS="\t"}{ if  ($13=="Depth" || $13 >= 4)  print $0}'> ${i}.reported_SNV_all.txt
##filter synonymous

cat ${i}.SNV_dep.txt | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($13=="Depth" || $13 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i}.reported_SNV_nonsynonymous.txt
#cat *nonsynonymous.txt | sort -k 1,1d -k 2,2n -k 5,5d | cut -f 1-14 | uniq -c | sort -k1,1nr > common_nonsynonymous.txt

##filter synonymous and common snp
#cat ${i}.SNV_dep.txt | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" || $11 < 0.01)  && ($17=="Depth" || $17 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i}.reported_SNV_nonsynonymous_remove_common.txt
number=`ls -l *vcf | wc -l`
if [ $number == $sample_number ]
then
  find *bam | parallel 'bedtools coverage -hist -abam {} -b /u/home/d/dingxm/mm10_coding_exon_merge.bed | grep ^all > {}.hist.all.txt'
  Rscript ~/mouse_exome_coverage.R 
  mkdir report_of_sample_SNV
  mkdir report_of_sample_vcf
  mkdir report_of_exome_coverage
  cp *reported* report_of_sample_SNV
  cp *vcf report_of_sample_vcf
  cp exome-coverage* report_of_exome_coverage
fi
  
#!/bin/bash
reference=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa
dbsnp=/u/nobackup/lixuser/data/dbsnp_138.hg19.vcf
gold_indel=/u/nobackup/lixuser/data/Mills_and_1000G_gold_standard.indels.hg19.vcf
human_genome_fa="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human_coding_bed="/u/home/d/dingxm/hg19_coding_merge.bed"
cancer_genes="/u/home/d/dingxm/cancer_genes.txt"
mkdir BWA_Alignment
i=$1
sample_number=$2
echo `date`
RG='@RG\tID:'$i'\tSM:'$i'\tLB:'$i'\tPL:illumina'
bwa mem $reference -t 8 -M -R "$RG" $i/*R1* $i/*R2* | samtools view -bS - > BWA_Alignment/${i}.bam
#single read
#bwa mem $reference -t 4 -M -R "$RG" $i/${i}*R1* | samtools view -bS - > BWA_Alignment/${i}.bam

#  #describe=`echo ${i} | sed 's/.bam//'`
samtools sort BWA_Alignment/${i}.bam BWA_Alignment/${i}.sorted
#i=$1
#  java -Xmx8g -jar ~/picard-tools-1.115/AddOrReplaceReadGroups.jar INPUT=BWA_Alignment/${i}.sorted.bam OUTPUT=BWA_Alignment/${i}.sorted.RG.bam SORT_ORDER=coordinate RGLB=8 RGPL=Illumina RGPU=1 RGSM=$i # for tophat add group
#  java -Xmx8g -jar ~/picard-tools-1.115/ReorderSam.jar INPUT=BWA_Alignment/${i}.sorted.RG.bam OUTPUT=BWA_Alignment/${i}.sorted.RE.bam REFERENCE=$reference #RNA_seq mutaion

java -Xmx12g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=BWA_Alignment/${i}.dup.txt INPUT=BWA_Alignment/${i}.sorted.bam ASSUME_SORTED=True OUTPUT=BWA_Alignment/${i}.remove.bam
  #java -Xmx12g -jar ~/picard-tools-1.141/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true METRICS_FILE=BWA_Alignment/${i}.dup.txt INPUT=BWA_Alignment/${i}.sorted.bam ASSUME_SORTED=True OUTPUT=BWA_Alignment/${i}.remove.bam
# java -jar ~/GenomeAnalysisTK.jar -T SplitNCigarReads -R $reference -I BWA_Alignment/${i}.remove.bam -o BWA_Alignment/${i}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS #RNA_seq mutaion
 rm BWA_Alignment/${i}.bam
 
#  
  samtools index BWA_Alignment/${i}.remove.bam
##  samtools index BWA_Alignment/${i}.split.bam
  echo `date`
  java -Xmx12g -jar ~/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -known $gold_indel -L ~/hg19_coding_re.bed -ip 100 -o BWA_Alignment/${i}.remove.bam.list -I BWA_Alignment/${i}.remove.bam 
  echo `date`
  java -Xmx12g -jar ~/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -known $gold_indel -I BWA_Alignment/${i}.remove.bam -targetIntervals BWA_Alignment/${i}.remove.bam.list -o BWA_Alignment/${i}.realigned.bam 
  echo `date`
  java -Xmx12g -jar ~/GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference -L ~/hg19_coding_re.bed -ip 100 -knownSites $dbsnp -knownSites $gold_indel -I BWA_Alignment/${i}.realigned.bam  -o BWA_Alignment/${i}_recal_table -nct 4
  echo `date`
  java -Xmx12g -jar ~/GenomeAnalysisTK.jar -T PrintReads -R $reference -I BWA_Alignment/${i}.realigned.bam  -o BWA_Alignment/${i}.realigned.recal.bam -BQSR BWA_Alignment/${i}_recal_table -nct 4
  echo `date`
  
 
 rm BWA_Alignment/${i}.remove.bam
 rm BWA_Alignment/${i}.sorted.bam
 rm BWA_Alignment/${i}.realigned.bam
 rm BWA_Alignment/${i}.remove.bam.list
 rm BWA_Alignment/${i}_recal_table
 samtools index BWA_Alignment/${i}.realigned.recal.bam
##find *bam | parallel 'bamToBed -i BAM {} | awk '{print $1,$2,$3,$6}' > {}_reads.bed'
cd BWA_Alignment
samtools mpileup -C50 -DSuf $human_genome_fa -l $human_coding_bed ${i}.realigned.recal.bam | bcftools view -bvcg - > ${i}_var_raw.bcf
#samtools mpileup -C50 -DSuf $mouse_genome_fa -l $mouse_coding_bed *.bam | bcftools view -bvcg - > var_raw.bcf
bcftools view ${i}_var_raw.bcf | vcfutils.pl varFilter -d 4 -D 1000 > ${i}_var_filter.vcf
#SNV annotation
convert2annovar.pl -format vcf4 ${i}_var_filter.vcf -outfile ${i}.avinput
table_annovar.pl ${i}.avinput /u/nobackup/lixuser/data/annovar/humandb -buildver hg19 -out ${i} -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
echo -e "Chr\tStart\tEnd\tRef\tEnd\tType\tQuality\tDepth" > aviheader
echo "ref-forward:ref-reverse:alt-forward:alt-reverse" > header
echo "Genotype" > header2 
cat aviheader ${i}.avinput > ${i}_head
cut -f 7-8 ${i}_head | paste ${i}.hg19_multianno.txt - > ${i}.SNV_dep2.txt
cat ${i}_var_filter.vcf | grep -v '^#' | cut -f 10 | cut -f 1 -d ":" | cat header2 - | paste ${i}.SNV_dep2.txt - > ${i}.SNV_dep3.txt
cat ${i}_var_filter.vcf | grep -oP '(?<=DP4=)[^;]*' | sed 's/,/:/g' | cat header - | paste ${i}.SNV_dep3.txt - > ${i}.SNV_dep.txt
rm ${i}.SNV_dep2.txt
rm ${i}.SNV_dep3.txt
##filter synonymous
#for i in *dep.txt; do cat $i | grep -v rs[0-9]* | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" )  && ($16=="Depth" || $16 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV.txt; done
cat ${i}.SNV_dep.txt |  awk 'BEGIN{FS=OFS="\t"}{ if  ($17=="Depth" || $17 >= 4)  print $0}'> ${i}.reported_SNV_all.txt
##filter synonymous

cat ${i}.SNV_dep.txt | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($17=="Depth" || $17 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i}.reported_SNV_nonsynonymous.txt
#get cancer genes SNV
awk 'BEGIN{FS=OFS="\t"} NR==FNR {list1[$1]=$1;next}{if ($7 in list1) print $0}' $cancer_genes ${i}.reported_SNV_nonsynonymous.txt > ${i}.reported_SNV_nonsynonymous_cancer_genes.txt

#cat *nonsynonymous.txt | sort -k 1,1d -k 2,2n -k 5,5d | cut -f 1-14 | uniq -c | sort -k1,1nr > common_nonsynonymous.txt

##filter synonymous and common snp
#cat ${i}.SNV_dep.txt | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" || $11 < 0.01)  && ($17=="Depth" || $17 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i}.reported_SNV_nonsynonymous_remove_common.txt
number=`ls -l *vcf | wc -l`
if [ $number == $sample_number ]
then
  find *recal.bam | parallel '/u/nobackup/lixuser/data/CAP_miRNA/bin/bedtools coverage -hist -abam {} -b ~/hg19_coding_merge.bed | grep ^all > {}.hist.all.txt'
  Rscript ~/exome_coverge_plot.R 
  mkdir report_of_sample_SNV
  mkdir report_of_sample_vcf
  mkdir report_of_exome_coverage
  cp *reported* report_of_sample_SNV
  cp *vcf report_of_sample_vcf
  cp exome-coverage* report_of_exome_coverage
fi
#cat 4.reported_SNV_nonsynonymous.txt 7.reported_SNV_nonsynonymous.txt 11.reported_SNV_nonsynonymous.txt | sort -k 1,1d -k 2,2n -k 5,5d | cut -f 1-15,18 | uniq -c | sort -k1,1nr | grep "^[[:space:]]*3" | cut -f 8- -d " " > dis_com_snp.txt
#cat An.hg19_multianno.vcf | grep -oP '(?<=ANNOVAR_DATE=2014-07-22;).*' | grep -oP '.*(?=ALLELE_END)' > genotype_anno.txt











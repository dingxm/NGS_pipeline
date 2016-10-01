#!/bin/bash

i=$1
sampleid_control=$2
tumor_purity=$3
reference=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa
human_genome_fa="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human_coding_bed="/u/home/d/dingxm/hg19_coding_merge.bed"
#mouse_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
#mouse_coding_bed="/u/home/d/dingxm/mm10_coding_exon.bed"
dbsnp=/u/nobackup/lixuser/data/dbsnp_138.hg19.vcf
gold_indel=/u/nobackup/lixuser/data/Mills_and_1000G_gold_standard.indels.hg19.vcf
cosmic_vcf=/u/home/d/dingxm/b37_cosmic_v54_120711.vcf

#mkdir BWA_Alignment
#for i in $1 $2
#do 
# echo `date`
# RG='@RG\tID:'$i'\tSM:'$i'\tLB:'$i'\tPL:illumina'
# bwa mem $reference -t 4 -M -R "$RG" $i/${i}*R1* $i/${i}*R2* | samtools view -bS - > BWA_Alignment/${i}.bam
##single read
##bwa mem $reference -t 4 -M -R "$RG" $i/${i}*R1* | samtools view -bS - > BWA_Alignment/${i}.bam
##  #describe=`echo ${i} | sed 's/.bam//'`
# samtools sort BWA_Alignment/${i}.bam BWA_Alignment/${i}.sorted
##  java -Xmx8g -jar ~/picard-tools-1.115/AddOrReplaceReadGroups.jar INPUT=BWA_Alignment/${i}.sorted.bam OUTPUT=BWA_Alignment/${i}.sorted.RG.bam SORT_ORDER=coordinate RGLB=8 RGPL=Illumina RGPU=1 RGSM=$i # for tophat add group
##  java -Xmx8g -jar ~/picard-tools-1.115/ReorderSam.jar INPUT=BWA_Alignment/${i}.sorted.RG.bam OUTPUT=BWA_Alignment/${i}.sorted.RE.bam REFERENCE=$reference
#  java -Xmx8g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=BWA_Alignment/${i}.dup.txt INPUT=BWA_Alignment/${i}.sorted.bam ASSUME_SORTED=True OUTPUT=BWA_Alignment/${i}.remove.bam
## java -jar ~/GenomeAnalysisTK.jar -T SplitNCigarReads -R $reference -I BWA_Alignment/${i}.remove.bam -o BWA_Alignment/${i}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
# rm BWA_Alignment/${i}.bam
##  
#  samtools index BWA_Alignment/${i}.remove.bam
###  samtools index BWA_Alignment/${i}.split.bam
#  echo `date`
#  java -Xmx8g -jar ~/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -known $gold_indel -L ~/hg19_coding_re.bed -ip 100 -o BWA_Alignment/${i}.remove.bam.list -I BWA_Alignment/${i}.remove.bam 
#  echo `date`
#  java -Xmx8g -jar ~/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -known $gold_indel -I BWA_Alignment/${i}.remove.bam -targetIntervals BWA_Alignment/${i}.remove.bam.list -o BWA_Alignment/${i}.realigned.bam 
#  echo `date`
#  java -Xmx8g -jar ~/GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference -L ~/hg19_coding_re.bed -ip 100 -knownSites $dbsnp -knownSites $gold_indel -I BWA_Alignment/${i}.realigned.bam  -o BWA_Alignment/${i}_recal_table
#  echo `date`
#  java -Xmx8g -jar ~/GenomeAnalysisTK.jar -T PrintReads -R $reference -I BWA_Alignment/${i}.realigned.bam  -o BWA_Alignment/${i}.realigned.recal.bam -BQSR BWA_Alignment/${i}_recal_table
#  echo `date`
#  
## 
# rm BWA_Alignment/${i}.remove.bam
# rm BWA_Alignment/${i}.sorted.bam
# rm BWA_Alignment/${i}.realigned.bam
# rm BWA_Alignment/${i}.remove.bam.list
# rm BWA_Alignment/${i}_recal_table
# samtools index BWA_Alignment/${i}.realigned.recal.bam
#done 
cd BWA_Alignment
mkdir ${i}_${sampleid_control}
cd ${i}_${sampleid_control}
samtools mpileup -B -q 1 -f $human_genome_fa -l $human_coding_bed ../${sampleid_control}.realigned.recal.bam > ${sampleid_control}.realigned.recal.bam.mpileup
samtools mpileup -B -q 1 -f $human_genome_fa -l $human_coding_bed ../${i}.realigned.recal.bam > ${i}.realigned.recal.bam.mpileup
mkdir ${i}_${sampleid_control}_vcf   
   java -Xmx8g -jar ~/VarScan.v2.3.7.jar somatic ${sampleid_control}.realigned.recal.bam.mpileup ${i}.realigned.recal.bam.mpileup ${i}_${sampleid_control}_vcf/${i} --min-coverage 10 --min-var-freq 0.20 --somatic-p-value 0.05 --tumor-purity $tumor_purity --output-vcf 1
cd ${i}_${sampleid_control}_vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar processSomatic ${i}.snp.vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar processSomatic ${i}.indel.vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.snp.Somatic.hc.vcf --min-var-freq 0.06 --indel-file ${i}.indel.vcf --output-file ${i}.snp.filter.vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.snp.LOH.hc.vcf --indel-file ${i}.indel.vcf --output-file ${i}.LOH.filter.vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.indel.Somatic.hc.vcf --output-file ${i}.indel.filter.vcf
  java -Xmx8g -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.indel.LOH.hc.vcf --output-file ${i}.indel.LOH.filter.vcf
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.snp.filter.vcf > ${i}.snp.filter.vcf.count
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.LOH.filter.vcf > ${i}.LOH.filter.vcf.count
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.indel.filter.vcf > ${i}.indel.filter.vcf.count
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.indel.LOH.filter.vcf > ${i}.indel.LOH.filter.vcf.count
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.snp.filter.vcf.count ../../${i}.realigned.recal.bam > ${i}.snp.readcount
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.LOH.filter.vcf.count ../../${sampleid_control}.realigned.recal.bam > ${i}.LOH.readcount
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.indel.filter.vcf.count ../../${i}.realigned.recal.bam > ${i}.indel.readcount
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.indel.LOH.filter.vcf.count ../../${sampleid_control}.realigned.recal.bam > ${i}.indel.LOH.readcount
  perl ~/scripts/fpfilter.pl --var-file ${i}.snp.filter.vcf --readcount-file ${i}.snp.readcount --output-file ${i}.snp.fpfilter
  perl ~/scripts/fpfilter.pl --var-file ${i}.LOH.filter.vcf --readcount-file ${i}.LOH.readcount --output-file ${i}.LOH.fpfilter
  perl ~/scripts/fpfilter.pl --var-file ${i}.indel.filter.vcf --readcount-file ${i}.indel.readcount --output-file ${i}.indel.fpfilter
  perl ~/scripts/fpfilter.pl --var-file ${i}.indel.LOH.filter.vcf --readcount-file ${i}.indel.LOH.readcount --output-file ${i}.indel.LOH.fpfilter
  convert2annovar.pl -format vcf4 ${i}.snp.filter.vcf -allsample -outfile ${i}.snp
  convert2annovar.pl -format vcf4 ${i}.LOH.filter.vcf -allsample -outfile ${i}.LOH
  convert2annovar.pl -format vcf4 ${i}.indel.filter.vcf -allsample -outfile ${i}.indel
  convert2annovar.pl -format vcf4 ${i}.indel.LOH.filter.vcf -allsample -outfile ${i}.indel.LOH
  table_annovar.pl ${i}.LOH.TUMOR.avinput /u/nobackup/lixuser/data/annovar/humandb/ -buildver hg19 -out ${i}_LOH -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
  table_annovar.pl ${i}.snp.TUMOR.avinput /u/nobackup/lixuser/data/annovar/humandb/ -buildver hg19 -out ${i}_snp -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
  table_annovar.pl ${i}.indel.TUMOR.avinput /u/nobackup/lixuser/data/annovar/humandb/ -buildver hg19 -out ${i}_indel -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
  table_annovar.pl ${i}.indel.LOH.TUMOR.avinput /u/nobackup/lixuser/data/annovar/humandb/ -buildver hg19 -out ${i}_indel_LOH -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
#  annotate_variation.pl -geneanno -buildver hg19 ${i}.snp.TUMOR.avinput /home/dingxm/NGS_software/annovar/humandb/
#  annotate_variation.pl -geneanno -buildver hg19 ${i}.LOH.TUMOR.avinput /home/dingxm/NGS_software/annovar/humandb/
  echo -e "Depth\tRAF\tVAF" > header1
  grep -w "Start" ${i}_snp.hg19_multianno.txt | paste - header1 > header
  
  grep -vw "Start" ${i}_snp.hg19_multianno.txt >  ${i}.snp.result
  grep -vw "Start" ${i}_LOH.hg19_multianno.txt >  ${i}.LOH.result
  grep -vw "Start" ${i}_indel.hg19_multianno.txt >  ${i}.indel.result
  grep -vw "Start" ${i}_indel_LOH.hg19_multianno.txt >  ${i}.indel.LOH.result
#  cat ~/mutation_header ${i}.snp.result > ${i}_novel_mutation.txt
#  cat ~/mutation_header ${i}.LOH.result > ${i}_novel_LOH.txt
#  egrep -w 'PASS|FILTER' ${i}.fpfilter > ${i}.filter.pass
#  egrep -w 'PASS|FILTER' ${i}.LOH.fpfilter > ${i}.LOH.filter.pass
  grep -w "PASS" ${i}.snp.fpfilter | awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' > ${i}.pass
  grep -w "PASS" ${i}.LOH.fpfilter | awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' > ${i}.LOH.pass
  grep -w "PASS" ${i}.indel.fpfilter | awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' > ${i}.indel.pass
  grep -w "PASS" ${i}.indel.LOH.fpfilter | awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' > ${i}.indel.LOH.pass
  awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' ${i}.snp.result > ${i}.mutation
  awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' ${i}.LOH.result > ${i}.LOH.mutation
  awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' ${i}.indel.result > ${i}.indel.mutation
  awk 'BEGIN{FS=OFS="\t"}{print $0, $1 "-" $2}' ${i}.indel.LOH.result > ${i}.indel.LOH.mutation
  awk 'BEGIN{FS=OFS="\t"} NR==FNR {list1[$10]=$10;list2[$10]=$0;next}{if ($16 in list1) print $0,list2[$16]}' ${i}.pass ${i}.mutation > ${i}.snp.temp
  awk 'BEGIN{FS=OFS="\t"} NR==FNR {list1[$10]=$10;list2[$10]=$0;next}{if ($16 in list1) print $0,list2[$16]}' ${i}.LOH.pass ${i}.LOH.mutation > ${i}.LOH.temp
  awk 'BEGIN{FS=OFS="\t"} NR==FNR {list1[$10]=$10;list2[$10]=$0;next}{if ($16 in list1) print $0,list2[$16]}' ${i}.indel.pass ${i}.indel.mutation > ${i}.indel.temp
  awk 'BEGIN{FS=OFS="\t"} NR==FNR {list1[$10]=$10;list2[$10]=$0;next}{if ($16 in list1) print $0,list2[$16]}' ${i}.indel.LOH.pass ${i}.indel.LOH.mutation > ${i}.indel.LOH.temp
#  sort -k15 ${i}.mutation > ${i}.mutation.sort
#  sort -k15 ${i}.LOH.mutation > ${i}.LOH.mutation.sort
#  sort -k10 ${i}.pass > ${i}.pass.sort
#  sort -k10 ${i}.LOH.pass > ${i}.LOH.pass.sort
#  join -1 15 -2 10 ${i}.mutation.sort ${i}.pass.sort | sed "s/ /\t/g" > ${i}.snp.temp
#  join -1 15 -2 10 ${i}.LOH.mutation.sort ${i}.LOH.pass.sort | sed "s/ /\t/g" > ${i}.LOH.temp
  cut -f 1-15 ${i}.snp.temp > snptemp
  cut -f 21-23 ${i}.snp.temp | paste snptemp - > temp1
  cut -f 1-15 ${i}.LOH.temp > LOHtemp
  cut -f 21-23 ${i}.LOH.temp | paste LOHtemp - > temp2 
  cut -f 1-15 ${i}.indel.temp > indeltemp
  cut -f 21-23 ${i}.indel.temp | paste indeltemp - > temp3
  cut -f 1-15 ${i}.indel.LOH.temp > indelLOHtemp
  cut -f 21-23 ${i}.indel.LOH.temp | paste indelLOHtemp - > temp4 
  #mv temp1 ${i}.snp.temp
  #mv temp2 ${i}.LOH.temp
  cat header temp1 > ${i}_novel_mutation.txt
  cat header temp2 > ${i}_novel_LOH.txt
  cat header temp3 > ${i}_novel_indel.mutation.txt
  cat header temp4 > ${i}_novel_indle.LOH.txt
#  cp ${i}_novel_mutation.txt ../../../result
#  cp ${i}_novel_LOH.txt ../../../result
  #java -jar ~/VarScan.v2.3.7.jar processSomatic ${i%.mpileup}.indel.vcf
  cd ..
# MutTect call
 
  
#done
#for i in $sampleid
#do
#   mkdir ${i%.mpileup}_vcf
#   java -jar ~/VarScan.v2.3.7.jar somatic Sample_GM00954Fibroblasts.mpileup $i ${i%.mpileup}_vcf/${i%.mpileup} --min-coverage 10 --min-var-freq 0.08 --somatic-p-value 0.05 --output-vcf 1
#  cd ${i%.mpileup}_vcf
#  java -jar ~/VarScan.v2.3.7.jar processSomatic ${i%.mpileup}.snp.vcf
#  java -jar ~/VarScan.v2.3.7.jar somaticFilter ${i%.mpileup}.snp.Somatic.hc.vcf --indel-file ${i%.mpileup}.indel.vcf --output-file ${i%.mpileup}.snp.filter.vcf
#  convert2annovar.pl -format vcf4 ${i%.mpileup}.snp.filter.vcf -allsample -outfile ${i%.mpileup}.snp
#  annotate_variation.pl -geneanno -buildver hg19 ${i%.mpileup}.snp.TUMOR.avinput ~/annovar/humandb/
#  grep -v -w "synonymous" ${i%.mpileup}.snp.TUMOR.avinput.exonic_variant_function >  ${i%.mpileup}.snp.result
#  cat ~/mutation_header ${i%.mpileup}.snp.result > ${i%.mpileup}_novel_mutation.txt
#  cp ${i%.mpileup}_novel_mutation.txt ../../../result
#  #java -jar ~/VarScan.v2.3.7.jar processSomatic ${i%.mpileup}.indel.vcf
#  cd ..
#done
#for i in $sampleid
#do
#   cd ..
#   samtools index ${i%.mpileup}.remove.bam 
#  cd mpileup/${i%.mpileup}_vcf
# 
#  cp ${i%.mpileup}.fpfilter ../../../result
#  cd ..
#done

#for i in $sampleid
#do
#  cd ${i%.mpileup}_vcf
#  egrep -w 'PASS|FILTER' ${i%.mpileup}.fpfilter > ${i%.mpileup}.filter.pass
#  grep -w "PASS" ${i%.mpileup}.fpfilter | awk '{print $0, $1 "-" $2}' > ${i%.mpileup}.pass
#  awk '{print $0, $5 "-" $6}' ${i%.mpileup}.snp.result > ${i%.mpileup}.mutation
#  sort -k13 ${i%.mpileup}.mutation > ${i%.mpileup}.mutation.sort
#  sort -k10 ${i%.mpileup}.pass > ${i%.mpileup}.pass.sort
#  join -1 13 -2 10 ${i%.mpileup}.mutation.sort ${i%.mpileup}.pass.sort | sed "s/ /\t/g" > ${i%.mpileup}.temp
#  cut -f 3-11 ${i%.mpileup}.temp > temp
#  mv temp ${i%.mpileup}.temp
#  cat ~/mutation_header ${i%.mpileup}.temp > ${i%.mpileup}_novel_mutation.txt
#  cp ${i%.mpileup}.filter.pass ../../../result
#  cp ${i%.mpileup}_novel_mutation.txt ../../../result
#  cd .. 
#done

#for i in $sampleid
#do
##   mkdir ${i%.mpileup}_cnv
#    
#   java -jar ~/VarScan.v2.3.7.jar copynumber Sample_GM00954Fibroblasts.mpileup2 $i ${i%.mpileup}_cnv/${i%.mpileup}.new --min-coverage 10 --min-segment-size 20 --max-segment-size 100
#  cd ${i%.mpileup}_cnv
#  java -jar ~/VarScan.v2.3.7.jar copyCaller ${i%.mpileup}.new.copynumber --output-file ${i%.mpileup}.new.called --output-homdel-file ${i%.mpileup}.new.homdel
#  sed 's/tumor_depth/{i%.mpileup}_depth/' {i%.mpileup}.new.called > temp
#  mv temp {i%.mpileup}.new.called
#  cp ${i%.mpileup}.new.called ../../../result
#  cd ..
#done
#for i in $sampleid
#do
#  samtools mpileup -B -q 1 -f $genome_fa -l ~/hg19_coding.bed Sample_GM00954Fibroblasts.remove.bam $i > mpileup/${i%.remove.bam}_cnv/${i%.remove.bam}.mpileup2 
#  grep -wv 0 mpileup/${i%.remove.bam}_cnv/${i%.remove.bam}.mpileup2 > mpileup/${i%.remove.bam}_cnv/${i%.remove.bam}.re.mpileup
  #sed -i 's/	0	/	1	/g' mpileup/${i%.remove.bam}_cnv/${i%.remove.bam}.mpileup2
   #java -jar ~/VarScan.v2.3.7.jar copynumber  mpileup/${i%.remove.bam}_cnv/${i%.remove.bam}.re.mpileup ${i%.remove.bam}.re --mpileup 1
#   java -jar ~/VarScan.v2.3.7.jar copyCaller ${i%.remove.bam}.re.copynumber --output-file ${i%.remove.bam}.re.called --output-homdel-file ${i%.remove.bam}.re.homdel --recenter-up 0.05
#  sed -i 's/tumor_depth/{i%.remove.bam}_depth/' ${i%.remove.bam}.re.called 
#  cp ${i%.remove.bam}.re.called ../result
#done
#
#
#for i in $sampleid
#do
#sample_name=${i%.re.called.segments}
#grep -wv NA $i > ${sample_name}.re.called.segm
#perl ~/mergeSegments.pl ${sample_name}.re.called.segm --ref-arm-sizes ~/arm_sizes2.txt --output-basename ${sample_name}
#rm ${sample_name}.re.called.segm
#done

#samtools view -h INPUT.bam | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > OUTPUT.bam
#find *bam | parallel 'samtools mpileup -B -q 1 -f /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -l ~/hg19_coding_merge.bed {} > mpileup/{}.mpileup'



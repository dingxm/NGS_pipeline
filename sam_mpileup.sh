#!/bin/bash
#samtools mpileup -C50 -DSuf ../genome.fa -l ~/mm10_coding_exon.bed *remove.bam | bcftools view -bvcg - > ../result/var_raw.bcf
sampleid_control=$1
i=$2
tumor_purity=$3
human_genome_fa="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human_coding_bed="/u/home/d/dingxm/hg19_coding_merge.bed"
mouse_genome_fa="/u/home/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
mouse_coding_bed="/u/home/d/dingxm/mm10_coding_exon.bed"
#for i in *bam
#do
#java -Xmx16g -jar ~/picard-tools-1.115/MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=${i%_accepted.bam}.dup.txt INPUT=$i ASSUME_SORTED=True OUTPUT=${i%_accepted.bam}.remove.bam
#done 
#samtools mpileup -C50 -DSuf $human_genome_fa -l $human_coding_bed *.realigned.recal.bam | bcftools view -bvcg - > var_raw.bcf
#samtools mpileup -C50 -DSuf $mouse_genome_fa -l $mouse_coding_bed *.bam | bcftools view -bvcg - > var_raw.bcf
#bcftools view var_raw.bcf | vcfutils.pl varFilter -D 5000 > var_filter.vcf

#for i in $sampleid
#do
#  #samtools index ${i%.mpileup}.sorted.bam
#  samtools mpileup -B -q 1 -f $human_genome_fa -l $human_coding_bed $i > mpileup/${i%.sorted.bam}.mpileup2
#done
#for i in $sampleid
#do
mkdir ${i}_${sampleid_control}
cd ${i}_${sampleid_control}
#samtools mpileup -B -q 1 -f $human_genome_fa -l $human_coding_bed ../${sampleid_control}.realigned.recal.bam > ${sampleid_control}.realigned.recal.bam.mpileup
#samtools mpileup -B -q 1 -f $human_genome_fa -l $human_coding_bed ../${i}.realigned.recal.bam > ${i}.realigned.recal.bam.mpileup
   mkdir ${i}_${sampleid_control}_vcf   
   java -jar ~/VarScan.v2.3.7.jar somatic ${sampleid_control}.realigned.recal.bam.mpileup ${i}.realigned.recal.bam.mpileup ${i}_${sampleid_control}_vcf/${i} --min-coverage 10 --min-var-freq 0.04 --somatic-p-value 0.05 --tumor-purity $tumor_purity --output-vcf 1
   cd ${i}_${sampleid_control}_vcf
  java -jar ~/VarScan.v2.3.7.jar processSomatic ${i}.snp.vcf 
  java -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.snp.Somatic.hc.vcf --min-var-freq 0.05 --indel-file ${i}.indel.vcf --output-file ${i}.snp.filter.vcf
  java -jar ~/VarScan.v2.3.7.jar somaticFilter ${i}.snp.LOH.hc.vcf --indel-file ${i}.indel.vcf --output-file ${i}.LOH.filter.vcf
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.snp.filter.vcf > ${i}.snp.filter2.vcf.count
  perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${i}.LOH.filter.vcf > ${i}.LOH.filter.vcf.count
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.snp.filter2.vcf.count ../../${i}.realigned.recal.bam > ${i}.snp2.readcount
  bam-readcount -q 1 -b 20 -f $human_genome_fa -l ${i}.LOH.filter.vcf.count ../../${sampleid_control}.realigned.recal.bam > ${i}.LOH.readcount
  perl ~/scripts/fpfilter.pl --var-file ${i}.snp.filter.vcf --readcount-file ${i}.snp2.readcount --output-file ${i}.snp.fpfilter
  perl ~/scripts/fpfilter.pl --var-file ${i}.LOH.filter.vcf --readcount-file ${i}.LOH.readcount --output-file ${i}.LOH.fpfilter
  convert2annovar.pl -format vcf4 ${i}.snp.filter.vcf -allsample -outfile ${i}.snp
  convert2annovar.pl -format vcf4 ${i}.LOH.filter.vcf -allsample -outfile ${i}.LOH
  annotate_variation.pl -geneanno -buildver hg19 ${i}.snp.TUMOR.avinput /u/home/lixuser/data/annovar/humandb/
  annotate_variation.pl -geneanno -buildver hg19 ${i}.LOH.TUMOR.avinput /u/home/lixuser/data/annovar/humandb/
  grep -v -w "whatever" ${i}.snp.TUMOR.avinput.exonic_variant_function >  ${i}.snp.result
  grep -v -w "whatever" ${i}.LOH.TUMOR.avinput.exonic_variant_function >  ${i}.LOH.result
  cat ~/mutation_header ${i}.snp.result > ${i}_novel_mutation.txt
  cat ~/mutation_header ${i}.LOH.result > ${i}_novel_LOH.txt
#  egrep -w 'PASS|FILTER' ${i}.fpfilter > ${i}.filter.pass
#  egrep -w 'PASS|FILTER' ${i}.LOH.fpfilter > ${i}.LOH.filter.pass
  grep -w "PASS" ${i}.snp.fpfilter | awk '{print $0, $1 "-" $2}' > ${i}.pass
  grep -w "PASS" ${i}.LOH.fpfilter | awk '{print $0, $1 "-" $2}' > ${i}.LOH.pass
  awk '{print $0, $5 "-" $6}' ${i%}.snp.result > ${i}.mutation
  awk '{print $0, $5 "-" $6}' ${i%}.LOH.result > ${i}.LOH.mutation
  sort -k13 ${i}.mutation > ${i}.mutation.sort
  sort -k13 ${i}.LOH.mutation > ${i}.LOH.mutation.sort
  sort -k10 ${i}.pass > ${i}.pass.sort
  sort -k10 ${i}.LOH.pass > ${i}.LOH.pass.sort
  join -1 13 -2 10 ${i}.mutation.sort ${i}.pass.sort | sed "s/ /\t/g" > ${i}.snp.temp
  join -1 13 -2 10 ${i}.LOH.mutation.sort ${i}.LOH.pass.sort | sed "s/ /\t/g" > ${i}.LOH.temp
  cut -f 3-11 ${i}.snp.temp > temp1
  cut -f 3-11 ${i}.LOH.temp > temp2
  mv temp1 ${i}.snp.temp
  mv temp2 ${i}.LOH.temp
  cat ~/mutation_header ${i}.snp.temp > ${i}_novel_mutation.txt
  cat ~/mutation_header ${i}.LOH.temp > ${i}_novel_LOH.txt
#  cp ${i}_novel_mutation.txt ../../../result
#  cp ${i}_novel_LOH.txt ../../../result
  #java -jar ~/VarScan.v2.3.7.jar processSomatic ${i%.mpileup}.indel.vcf
  cd ..
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



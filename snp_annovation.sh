#!/bin/bash
#echo -n "Please type your first name: "
#read first_name
#echo -n "Please type your last name: "
#read last_name
#echo 
#echo "Your name is : $first_name $last_name "
#  convert2annovar.pl -format vcf4 JT_analysis.vcf -allsample -outfile JT
#  annotate_variation.pl -geneanno -buildver hg19 ${i%.mpileup}.snp.TUMOR.avinput ~/annovar/humandb/
#  grep -v -w "unknown|synonymous" ${i%.mpileup}.snp.TUMOR.avinput.exonic_variant_function >  ${i%.mpileup}.snp.result
#  cat ~/mutation_header ${i%.mpileup}.snp.result > ${i%.mpileup}_novel_mutation.txt
PI=$1   
#table_annovar.pl var_filter.vcf /u/project/lixuser/data/annovar/humandb/ -buildver hg19 -out An -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70,clinvar_20150330 -operation g,f,f,f,f,f -nastring na
#awk '{for(i=10;i<=NF;i++) {split($i,a,":");if (a[3] <3 | a[3]>500 ) break;}; if (i > NF) print  }' JT.hg19_multianno.vcf | wc -l
#
#grep -e "ExonicFunc.refGene=nonsynonymous_SNV;.*1000g2012apr_all=na.*snp138=na" JT_analysis_highdep.vcf > JT_analysis_filter.vcf
#
#awk  'BEGIN {ORS="\t"}{print $1"_"$2; for(i=10;i<=57;i++) print substr($i,1,3); printf "\n" }' test.vcf > genotype
##generage input file for fbat
#awk  'BEGIN {FS = "/";ORS="\t"};{for (i=1;i<=NF;i++) print $i; printf "\n"}' genotype
##transpose
#awk ' BEGIN {ORS="\t"};{for (i=1;i<=NF;i++) {array[NR,i]=$i; if(big <= NF) big=NF}};END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++)  {print array[j,i]}; printf "\n"} }'
#samtools view -h INPUT.bam | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > OUTPUT.bam

convert2annovar.pl -format vcf4 var_filter.vcf -allsample -outfile $PI
for i in *avinput; do table_annovar.pl $i /u/project/lixuser/data/annovar/humandb/ -buildver hg19 -out ${i%.avinput} -remove -protocol refGene,1000g2012apr_all,snp138NonFlagged,ljb23_sift,cosmic70 -operation g,f,f,f,f -nastring na; done
echo -e "Chr\tStart\tEnd\tRef\tEnd\tType\tQuality\tDepth" > aviheader
for i in *avinput; do cat aviheader $i > ${i}_head; done
for i in *head; do cut -f 7-8 $i | paste ${i%.avinput_head}.hg19_multianno.txt - > ${i%.avinput}.SNV_dep.txt; done
for i in *dep.txt; do cat $i | grep -v rs[0-9]* | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2012apr_all" )  && ($16=="Depth" || $16 >= 4) && ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV.txt; done
#cat lncipedia_3_1.gtf | cut -f 9 | awk 'BEGIN {ORS=" "}{for(i=1;i<37;i++) {if (i%3 ==1) print $i; if(i%3==2) print "\""$i"\"""\;";}; printf "\n"}' > linc_attr
#cat lncipedia_3_1.gtf | cut -f 1-8 | paste - linc_attr > hg19_lincRNA.gtf
#awk 'NR==FNR {list[$1]=$1;next}{if ($8 in list) print }' autism_gene.txt QCalhoun.samtools.snp.annovar.hg19_multianno.xls
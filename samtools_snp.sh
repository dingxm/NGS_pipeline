#!/bin/bash
human_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
human_coding_bed="/u/home/d/dingxm/hg19_coding_merge.bed"
mouse_genome_fa="/u/project/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
mouse_coding_bed="/u/home/d/dingxm/mm10_coding_exon.bed"
samtools mpileup -C50 -DSuf $human_genome_fa -l $human_coding_bed *.realigned.recal.bam | bcftools view -bvcg - > var_raw.bcf
#samtools mpileup -C50 -DSuf $mouse_genome_fa -l $mouse_coding_bed *.bam | bcftools view -bvcg - > var_raw.bcf
bcftools view var_raw.bcf | vcfutils.pl varFilter -D 5000 > var_filter.vcf
find *bam | parallel 'bedtools coverage -hist -abam {} -b ~/hg19_coding_merge.bed | grep ^all > {}.hist.all.txt'
Rscript ~/exome_coverge_plot.R 
#SNV annotation
convert2annovar.pl -format vcf4 var_filter.vcf -allsample -outfile Anno
for i in *avinput; do table_annovar.pl $i /u/project/lixuser/data/annovar/humandb -buildver hg19 -out ${i%.avinput} -remove -protocol refGene,1000g2014oct_all,snp138NonFlagged,ljb23_sift,cosmic70 -operation g,f,f,f,f -nastring na; done
echo -e "Chr\tStart\tEnd\tRef\tEnd\tType\tQuality\tDepth" > aviheader
for i in *avinput; do cat aviheader $i > ${i}_head; done
for i in *head; do cut -f 7-8 $i | paste ${i%.avinput_head}.hg19_multianno.txt - > ${i%.avinput}.SNV_dep.txt; done
##filter synonymous
#for i in *dep.txt; do cat $i | grep -v rs[0-9]* | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" )  && ($16=="Depth" || $16 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV.txt; done
for i in *dep.txt; do cat $i |  awk 'BEGIN{FS=OFS="\t"}{ if  ($16=="Depth" || $16 >= 4)  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV_all.txt; done
##filter synonymous

for i in *dep.txt; do cat $i | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($16=="Depth" || $16 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV_nonsynonymous.txt; done
cat *nonsynonymous.txt | sort -k 1,1d -k 2,2n -k 5,5d | cut -f 1-14 | uniq -c | sort -k1,1nr > common_nonsynonymous.txt

##filter synonymous and common snp
#for i in *dep.txt; do cat $i | grep -v rs[0-9]* | grep -vw synonymous | awk 'BEGIN{FS=OFS="\t"}{ if (($11 == "na" || $11 == "1000g2014oct_all" )  && ($16=="Depth" || $16 >= 4)&& ($6 == "Func.refGene" || $6 == "exonic" ))  print $0}'> ${i%.avinput_head.SNV_dep.txt}.reported_SNV_nonsynonymous.txt; done
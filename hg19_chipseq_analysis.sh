#!/bin/bash
hg19_reference=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
mm9_reference=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome
mm10_reference=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome

header=/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome.fa.fai

#for i in $1
#do 
#bowtie -S -m 1 -p 4 $mm9_reference -q $i/${i}* --sam-RG ID:$i --sam-RG SM:$i | samtools view -bS - > Bowtie_Alignment/${i}.bam
#bowtie -S -m 1 -p 4 $reference -q $i/${i}_R1* --sam-RG ID:$i --sam-RG SM:$i | samtools view -bS - > Bowtie_Alignment/${i}.bam
#done
target=$1
control=$2
mkdir ${target}_Bowtie_Alignment
#ratio=$3 
gunzip $target/*R1*
gunzip $control/*R1*
bowtie -S -m 1 -p 4 $hg19_reference -q $target/*R1* --sam-RG ID:$target --sam-RG SM:$target | samtools view -bS - > ${target}_Bowtie_Alignment/${target}.bam
bowtie -S -m 1 -p 4 $hg19_reference -q $control/*R1* --sam-RG ID:$control --sam-RG SM:$control | samtools view -bS - > ${target}_Bowtie_Alignment/${control}.bam
#bowtie -S -p 4 $mm10_reference -q *.fastq --sam-RG ID:$control --sam-RG SM:$control | samtools view -bS - > Bowtie_Alignment/${target}.bam
cd ${target}_Bowtie_Alignment
samtools sort ${target}.bam ${target}.sort
samtools sort ${control}.bam ${control}.sort
samtools index ${target}.sort.bam
samtools index ${control}.sort.bam
samstat ${target}.sort.bam
samstat ${control}.sort.bam
bedtools bamtobed -i ${target}.sort.bam | cut -f 1-3,6 > ${target}.read.bed
rm ${target}.bam
rm ${control}.bam
bedtools genomecov -ibam ${target}.sort.bam -g /u/home/d/dingxm/hg19_chrom_size.txt -bg > ${target}.bedgraph
~/scripts/bdf2bw.sh ${target}.bedgraph /u/home/d/dingxm/hg19_chrom_size.txt
#samtools sort ${control}.bam ${control}.bam.sort
#bedtools bamtobed -i ${target}.sort.bam > ${target}.bed
#awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' ${target}.bed >${target}_read.bed
#ngs.plot.r -G mm10 -R tss -C ${target}.bam -O $target -T $target -L 3000
#ngs.plot.r -G mm10 -R tss -C ${control}.bam -O ${control} -T $control -L 3000
macs2 callpeak -t ${target}.sort.bam -c ${control}.sort.bam -n ${target} -g hs -p 0.00001 --outdir MACS_result
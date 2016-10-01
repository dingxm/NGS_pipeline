#!/bin/bash
#pi=$1  # first augument
#dir_scy="/u/home/d/dingxm/workshop2_data/software/scythe-master"
#dir_sic="/u/home/d/dingxm/workshop2_data/software/sickle-master"
##mkdir ${pi}_fastqc
#for i in Sample_*/*trimmed*fastq
#do
#    ~/FastQC/fastqc -o ${pi}_fastqc $i
#done
## adaptor trimming
#cd ${pi}_analysis
#for i in Sample*/*.fastq
#do
#    pos1=`expr index "$i" "/"`
#    pos2=`expr index "$i" "."`
#    dir_sample=${i:0:pos1-1}
#    echo ${i:pos1:pos2-pos1-1}_scy.fastq
#    echo $dir_sample
##    $dir_scy/scythe -a ~/workshop2_data/adaptors.fa -q sanger -o $dir_sample/${i:pos1:pos2-pos1-1}_scy.fastq $i
#     java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 $i ${i:pos1:pos2-pos1-1}_trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:36
#done
sampleid=$1
mkdir LC3_alignment
mkdir ${sampleid}_meth
mkdir ${sampleid}_split
cd ${sampleid}
gunzip *.gz
mv *.fastq ${sampleid}_R1.fq
cd ..
/u/project/lixuser/data/methy-pipe2/cpp_prog/BSAligner -a $sampleid/*R1*.fq -D /u/project/lixuser/data/methy-pipe2/database/hg19 -r 0 -v 2 -o LC3_alignment/${sampleid}.bsalign -p 6 
/u/project/lixuser/data/methy-pipe2/cpp_prog/meth_call se /u/project/lixuser/data/methy-pipe2/database/hg19.W.ori.fa LC3_alignment/${sampleid}.bsalign ${sampleid}_meth/${sampleid}.W
perl -lane 'print if $F[-1]=~/C:G:/i' ${sampleid}_meth/${sampleid}.W.call > ${sampleid}_meth/${sampleid}.W.CpG.call 
/u/project/lixuser/data/methy-pipe2/cpp_prog/meth_call se /u/project/lixuser/data/methy-pipe2/database/hg19.C.ori.fa LC3_alignment/${sampleid}.bsalign ${sampleid}_meth/${sampleid}.C
/u/project/lixuser/data/methy-pipe2/cpp_prog/reverse_crick_metCall /u/project/lixuser/data/methy-pipe2/database/hg19.size ${sampleid}_meth/${sampleid}.C.call ${sampleid}_meth/${sampleid}.C.rev.call
perl -lane 'print if $F[-1]=~/C:G:/i'  ${sampleid}_meth/${sampleid}.C.rev.call > ${sampleid}_meth/${sampleid}.C.rev.CpG.call 
perl /u/project/lixuser/data/methy-pipe2/utils/split_meth_call/split_call.pl ${sampleid}_meth/${sampleid}.W.CpG.call ${sampleid}_meth/${sampleid}.C.rev.CpG.call ${sampleid}_split/${sampleid} /u/project/lixuser/data/methy-pipe2/bed_files
rm ${sampleid}_meth/${sampleid}.W.call
rm ${sampleid}_meth/${sampleid}.C.call
rm ${sampleid}_meth/${sampleid}.C.rev.call


#!/bin/bash
human_repeats="/home/dingxm/NGS_software/DCC-0.3.4/DCC/data/hg19_repeats.gtf"
human_anno="/home/dingxm/genome_reference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf"
human_genome_fa="/u/nobackup/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
mkdir circ_result
samplesheet=/home/dingxm/NGS_software/DCC-0.3.4/DCC/data/samplesheet
mate1=/home/dingxm/NGS_software/DCC-0.3.4/DCC/data/mate1
mate2=/home/dingxm/NGS_software/DCC-0.3.4/DCC/data/mate2
DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -D -R $human_repeats -an $human_anno -Pi -F -M -Nr 3 3 -fg -G -A $human_genome_fa
#!/bin/bash
python $HOME/local/bin/chimerascan_index.py --bowtie-dir=/u/local/apps/bowtie/0.12.7 /u/home/lixuser/data/HiSeq2000_data/GenomeReference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ~/hg19ref_gene_model.txt /u/home/lixuser/data/hg19_chimera_index

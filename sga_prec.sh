#!/bin/bash
read1=$1
read2=$2
sga preprocess --pe-mode 1 $read1 $read2 > mygenome.fastq
sga index -a ropebwt --no-reverse -t 8 mygenome.fastq
sga preqc -t 8 mygenome.fastq > mygenome.preqc
/u/home/d/dingxm/sga/src/bin/sga-preqc-report.py mygenome.preqc /u/home/d/dingxm/sga/src/examples/preqc/[hs]*.preqc
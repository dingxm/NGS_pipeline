#!/bin/bash
hmmscan --cpu 6 --domtblout TrinotatePFAM.out /u/project/lixuser/data/Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log
signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep
tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out
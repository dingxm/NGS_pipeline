#!/bin/bash
#blastx -query Trinity.fasta -db /u/project/lixuser/data/uniprot_sprot.trinotate.pep -num_threads 6 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db /u/project/lixuser/data/uniprot_sprot.trinotate.pep -num_threads 6 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
blastx -query Trinity.fasta -db /u/project/lixuser/data/uniprot_uniref90.trinotate.pep -num_threads 6 -max_target_seqs 1 -outfmt 6 > uniref90.blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db /u/project/lixuser/data/uniprot_uniref90.trinotate.pep -num_threads 6 -max_target_seqs 1 -outfmt 6 > uniref90.blastp.outfmt6
#hmmscan --cpu 6 --domtblout TrinotatePFAM.out /u/project/lixuser/data/Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log
#signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep
#tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out

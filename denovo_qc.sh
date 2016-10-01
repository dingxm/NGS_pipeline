#!/bin/bash
#~/SOAPec_v2.01/bin/KmerFreq_HA -k 23 -t 12 -p licochi -l readlist -L 150
#~/SOAPec_v2.01/bin/Corrector_HA -k 23 -l 2 -Q 33 -t 8 licochi.freq.gz readlist
~/SOAPec_v2.01/bin/KmerFreq_AR -k 17 -t 10 -p Lico -q 33 readlist 
~/SOAPec_v2.01/bin/KmerFreq_AR -k 17 -s 6 -t 10 -p Lico -q 33 readlist 
~/SOAPec_v2.01/bin/Corrector_AR -k 17 -l 2 -K 17 -s 6 -L 2 -r 50 -t 10 -Q 33 Lico.freq.cz Lico.freq.cz.len Lico.space.freq.cz Lico.space.freq.cz.len readlist 
	  

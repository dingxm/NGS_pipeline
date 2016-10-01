#!/bin/bash
config=$1
prefix=$2
/u/home/lixuser/data/SOAPdenovo2/SOAPdenovo-63mer all  -s $config -K 63  -o  assemble_output/$prefix -p 8
#!/bin/bash
config_file=$1
output_name=$2
/u/home/lixuser/data/SOAPdenovo2/SOAPdenovo-63mer all -s $config_file -K 63 -o assemble_output/$output_name -p 12
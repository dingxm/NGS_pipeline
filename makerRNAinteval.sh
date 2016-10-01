#!/usr/bin/env bash
# make_rRNA.sh
# Kamil Slowikowski
# Created: December 12, 2014
# Original: https://gist.github.com/slowkow/b11c28796508f03cdf4b
# Modified: January 22, 2014 
# António Domingues 
# Added arguments - easier to run and choose species
# rRNA coordinates can be in Bed format

##########################
## Usage
##########################
usage() {
  echo "   Usage: `basename $0` genomeVersion rRNA.bed"
  echo "   Creates a rRNA interval list suitable for RNASEQ-QC"
  echo ""
  echo "example: `basename $0` hg19 rRNA.[bed/GTF]"
  echo ""
  echo "   The expected inputs are:"
  echo "      arg1=species as listed by UCSC genome browser"
  echo "      arg2=genomic locations of rRNAs either in bed format, or a Gencode GTF"
  echo "      Gencode files can be obtained from (hg19) ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
  echo "      Look for the chromosome sizes file in the current folder, e.g. hg19.chrom.sizes. If not present it will download the information from UCSC."
  exit 2
}

## are 2 arguments supplied?
   if [ "$#" != "2" ]
   then
      usage
      exit 1
   fi


## set the input variables
export genome=$1
rRNA_file=$2
chrom_sizes=$genome.chrom.sizes


# 1. Download Chromosome sizes from the UCSC genome browser


if [[ ! -s $chrom_sizes ]]
then
echo "File with chromosome sizes was not found."
echo "Downloading from UCSC..."
	mysql -N --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
"SELECT chrom,size FROM chromInfo ORDER BY size DESC;" $genome \
> $chrom_sizes
	echo "Done."
fi


# 2. rRNA interval_list file

## is it a bed file or a GTF?
filename=$(basename "$rRNA_file")
extension="${filename##*.}"
echo $extension

if [[ $extension == "bed" ]]
then
	echo "rRNA is in $extension format"
	# Output file suitable for Picard tools
	rRNA=$genome.rRNA.interval_list

	# Sequence names and lengths. (Must be tab-delimited.)
	perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:$ENV{'genome'}"' $chrom_sizes | \
	grep -v _ \
	>> $rRNA

	# Intervals for rRNA transcripts.
	cat $rRNA_file | \
	awk -v OFS='\t' '{print $1, $2, $3, $6, $4}' | \
	sort -k1V -k2n -k3n \
	>> $rRNA

	echo "rRNA interval file created. Have fun."
elif [[ $extension == "gtf" ]]
then
	echo "rRNA is in $extension format"

	# Output file suitable for Picard tools
	rRNA=$genome.rRNA.interval_list

	# Sequence names and lengths. (Must be tab-delimited.)
	perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:$ENV{'genome'}"' $chrom_sizes | \
	grep -v _ \
	>> $rRNA

	# Intervals for rRNA transcripts.
	grep 'gene_type "rRNA"' $rRNA_file | \
	awk '$3 == "transcript"' | \
	cut -f1,4,5,7,9 | \
	perl -lane '
	/transcript_id "([^"]+)"/ or die "no transcript_id on $.";
	print join "\t", (@F[0,1,2,3], $1)
	' | \
	sort -k1V -k2n -k3n \
	>> $rRNA

	echo "rRNA interval file created. Have fun."
else
	echo "Please make sure that the file format with the rNA locations is either bed or gtf, and the extension match 'bed' or 'gtf'."
fi



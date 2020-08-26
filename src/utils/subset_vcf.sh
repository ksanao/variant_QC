#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1
onefile=$2

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path) vcf_file (with path)\n\
Extract target and freference regions from vcf files\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 2 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"

outfile=${onefile/.vcf/_target.vcf}
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKESAMPLE\n" > $outfile
bedtools intersect -wa -a "$onefile" -b "$intdir"/${targets_bed/.gz/} >> $outfile

outfile=${onefile/.vcf/_onref.vcf}
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKESAMPLE\n" > $outfile
bedtools intersect -wa -a "$onefile" -b "$intdir"/${ref_bed/.gz/} >> $outfile

exit 0

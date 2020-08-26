#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1
bamfile=$2

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config bamfile (specify config with path, bam file is taken from raw dir)\n\
Mark duplicated and fix read groups\n\
Output in intdir, same bamfile name\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 2 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source $config

bamfile=${bamfile##*/}
sample_ref=${bamfile/.bam/}
bqsr_file="$interim_gatk"/"$sample_ref"bqsr.map
bqsr_bam="$sample_ref"_bqsr.bam

mkdir -p $interim_gatk

# Base quality score recalibration using known sites from golden reference
gatk BaseRecalibrator -I "$intdir"/$bamfile -R "$intdir"/ref.fasta -O $bqsr_file --known-sites "$bqsr"/"$ref_highconf" --QUIET
gatk ApplyBQSR -I "$intdir"/$bamfile -O "$intdir"/"$bqsr_bam" -bqsr $bqsr_file

exit 0

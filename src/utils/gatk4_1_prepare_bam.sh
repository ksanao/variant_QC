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

mkdir -p $interim_gatk
mkdir -p $stats_gatk
mkdir -p $stats_picard

# mark duplicates and sort
gatk MarkDuplicates -I "$rawdir"/$bamfile -M "$stats_gatk"/${bamfile/.bam/_dedup.metrics} -O $intdir/$bamfile
samtools index $intdir/$bamfile

sample_ref=${bamfile/.bam/}
picard AddOrReplaceReadGroups -I $intdir/$bamfile --RGID $sample_ref --RGLB "$sample_ref" --RGPL "Unknown" --RGPU "$sample_ref" --RGSM "FAKESAMPLE" -O $intdir/temp_"$bamfile"
mv $intdir/temp_"$bamfile" $intdir/"$bamfile"

exit 0

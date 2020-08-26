#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

in_bam=$1
out_bam=$2
regions_bed=$3

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` in_bam out_bam regions (with paths)\n\
Reheader bam file to keep in the header chromosomes containing bed file regions\n\
in_bam: input bam\n\
out_bam: output bam\n\
regions_bed: bed file\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 3 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

reg=`cut -f1 "$regions_bed" | sort | uniq | tr '\n' '%'`
reg=${reg%?}

echo "Reheading $in_bam: keeping locations defined in $regions_bed"
samtools view -H $in_bam \
| awk -v R=$reg '{split(R,a,"%"); if ($1=="@SQ"){for (i in a){if ($2=="SN:"a[i]){print $0}}}else{print $0}}' \
| samtools reheader - "$in_bam" > "$out_bam"

exit 0

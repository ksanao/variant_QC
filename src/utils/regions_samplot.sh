#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path)\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 1 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"

in_bam=`for onefile in ${sample_bam[@]}; do echo "$rawdir"/"$onefile"; done`
in_pref=`for onefile in ${sample_bam[@]}; do echo "${onefile/.bam/}"; done`
sorted_gff="$intdir"/"${ref_gff/gff3.gz/sorted.gff3.gz}"
ref_fasta="$intdir"/"$ref.fasta"
regions="$intdir"/${targets_bed/.gz/}

mkdir -p "$prodir"/samplots
for i in $(seq 1 `cat $regions | wc -l`)
do 
    chr=`awk -v L=$i 'FNR==L' $regions | cut -f1`
    start=`awk -v L=$i 'FNR==L' $regions | cut -f2`
    end=`awk -v L=$i 'FNR==L' $regions | cut -f3`
    outfile="$prodir"/samplots/"$chr":"$start"-"$end".png
    samplot plot -n `echo ${in_pref[@]}` -b `echo ${in_bam[@]}` -o "$outfile" \
                 -c $chr -s $start -e $end -r $ref_fasta -T $sorted_gff --same_yaxis_scales \
                 --common_insert_size --coverage_only --separate_mqual 55 --coverage_tracktype stack
done

exit 0

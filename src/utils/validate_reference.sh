#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path) vcf_file\n\
Validate benchmark data against targets\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 1 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"


bedtools intersect -wa -a "$intdir"/${targets_bed_all/.gz/} -b "$intdir"/${ref_bed/.gz/} | sort | uniq > "$intdir"/$targets_bed
bedtools intersect -v -a "$intdir"/${targets_bed_all/.gz/} -b "$intdir"/${ref_bed/.gz/} | sort | uniq > "$intdir"/$targets_bed_missed

# collect stats on total targets, targets with benchmark and targets without benchmark data
tot_targets=`cat "$intdir"/${targets_bed_all/.gz/} | wc -l`
bm_targets=`cat "$intdir"/$targets_bed | wc -l`
missing_targets=`cat "$intdir"/$targets_bed_missed | wc -l`
stat_file="$prodir"/target_benchmark.stats
printf "Total_targets\tBenchmarked_targets\tMissing_targets\n" > $stat_file
printf "$tot_targets\t$bm_targets\t$missing_targets\n" >> $stat_file

exit 0

#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1
bamfile=$2

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config bamfile (config with path)\n"

# Print usage: scriptname -options, if script invoked with no command-line args
if [ $# -ne 2 ]
then
	echo -e $USAGE
	exit $ERROR_ARGUMENTS
fi

source "$config"

export "PATH=/home/src/utils:$PATH"

bamfile=${bamfile##*/}
sample_ref=${bamfile/.bam/}

# gatk hartd filtering workflow
gatk4_1_prepare_bam.sh "$config" "$bamfile"
gatk4_2_bqsr.sh "$config" "$bamfile"
gatk4_3_call_variants.sh "$config" "$bamfile"
gatk4_4_filter.sh "$config" "$bamfile"

# visualise in IGV
sortedgff="${ref_gff/gff3.gz/sorted.gff3.gz}"
create_report "$intdir"/${targets_bed/.gz/} "$intdir"/ref.fasta --tracks "$intdir"/"$sample_ref"_bqsr.bam "$rawdir"/"$ref_vcf" "$intdir"/"$sortedgff" --output "$prodir"/igv_viewer_gatk.html

exit 0

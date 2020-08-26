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
vcf_hc="$interim_gatk"/"$sample_ref"_hc.g.vcf.gz
gatk_vcf="$interim_gatk"/"$sample_ref"_gatk.vcf.gz

# Variant calling using GATK HaplotypeCaller
gatk HaplotypeCaller -I "$intdir"/"$bqsr_bam" -R "$intdir"/ref.fasta -O $gatk_vcf --output-mode EMIT_ALL_CONFIDENT_SITES

# Genotype GVCFs
# gatk GenotypeGVCFs -V $vcf_hc -R "$intdir"/ref.fasta -O $gatk_vcf

tabix -fp vcf $gatk_vcf

exit 0

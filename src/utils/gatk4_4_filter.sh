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
        Gexit $ERROR_ARGUMENTS
fi

source $config
export "PATH=/home/src/utils:$PATH"

bamfile=${bamfile##*/}
sample_ref=${bamfile/.bam/}
gatk_vcf="$interim_gatk"/"$sample_ref"_gatk.vcf.gz
gatk_snp_vcf="$interim_gatk"/"$sample_ref"_gatk_snp.vcf.gz
gatk_filt_vcf="$intdir"/"$sample_ref"_gatk_filt.vcf
gatk_final_vcf="$intdir"/"$sample_ref"_gatk.vcf

gatk SelectVariants \
        -R "$intdir"/ref.fasta \
        -V "$gatk_vcf" \
        --select-type-to-include SNP \
        -O $gatk_snp_vcf

gatk VariantFiltration \
        -R "$intdir"/ref.fasta \
        -V $gatk_snp_vcf \
        -O $gatk_filt_vcf \
        -filter-name "QD_filter" -filter "QD < 5.0" \
        -filter-name "FS_filter" -filter "FS > 6.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 3.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -0.25" \
        -filter-name "ReadPosRankSum_filter_low" -filter "ReadPosRankSum < -6.0"

gatk SelectVariants \
        -R "$intdir"/ref.fasta \
        -V "$gatk_filt_vcf" \
        --exclude-filtered \
        -O $gatk_final_vcf

subset_vcf.sh $config $gatk_final_vcf

exit 0

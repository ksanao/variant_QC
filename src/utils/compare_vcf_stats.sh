#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1
onefile=$2

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path)\n\
Overlap stats for reference and target vcf files\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 1 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"

# subset target regions from bed files
statfile="$prodir"/target_variant.stats
printf "File\tOn_target\tOff_target\n" > $statfile
 
for onefile in ${sample_vcf[@]} $ref_vcf
do
    on_target=`awk 'FNR>1' "$intdir"/${onefile/.vcf.gz/_target.vcf} | wc -l`
    off_target=`bedtools intersect -v -a "$intdir"/${onefile/.gz/} -b "$intdir"/${targets_bed/.gz/} | wc -l`
    printf "$onefile\t$on_target\t$off_target\n" >> $statfile
done

# subset reference regions from bed files
statfile=$prodir/reference_variant.stats
printf "File\tOn_ref\tOff_ref\n" > $statfile

i=0
for onefile in ${sample_vcf[@]} $ref_vcf
do
    ((i++))
    on_target=`awk 'FNR>1' "$intdir"/${onefile/.vcf.gz/_onref.vcf} | wc -l`
    off_target=`bedtools intersect -v -a "$intdir"/${onefile/.gz/} -b "$intdir"/${ref_bed/.gz/} | wc -l`
    printf "$onefile\t$on_target\t$off_target\n" >> $statfile
done

statfile=$prodir/common_variants.stats
ontarget=`cut -f1,2,4,5 $intdir/*_target.vcf | sort | uniq -c | awk -v N=$i '$1==N' | grep -v '#CHROM' | wc -l`
onref=`cut -f1,2,4,5 $intdir/*_onref.vcf | sort | uniq -c | awk -v N=$i '$1==N' | grep -v '#CHROM' | wc -l`
printf "Region\tCommon_variants\tFiles\n" > $statfile
printf "Target\t$ontarget\t$i\n" >> $statfile
printf "Reference\t$onref\t$i\n" >> $statfile

for onefile in ${sample_vcf[@]}
do
    onref=`cat $intdir/${onefile/.vcf.gz/_onref.vcf} $intdir/${ref_vcf/.vcf.gz/_onref.vcf} | cut -f1,2,4,5 \
           | sort | uniq -c | awk '$1==2' | grep -v '#CHROM' | wc -l`
    ontarget=`cat $intdir/${onefile/.vcf.gz/_target.vcf} $intdir/${ref_vcf/.vcf.gz/_target.vcf} | cut -f1,2,4,5 \
           | sort | uniq -c | awk '$1==2' | grep -v '#CHROM' | wc -l`

    printf "Target\t$ontarget\t$onefile\n" >> $statfile
    printf "Reference\t$onref\t$onefile\n" >> $statfile
done

exit 0

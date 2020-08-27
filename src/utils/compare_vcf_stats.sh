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
gatk_vcf=`for onefile in ${sample_bam[@]}; do echo ${onefile/.bam/_gatk.vcf}; done`

# subset target regions from bed files
statfile="$prodir"/target_variant.stats
printf "File\tOn_target\tOff_target\n" > $statfile
 
for onefile in ${sample_vcf[@]} $ref_vcf ${gatk_vcf[@]}
do  
    in_target=${onefile/.gz/}
    in_target="$intdir"/${in_target/.vcf/_target.vcf}
    echo $in_target
    if [ -f "$in_target" ]
    then
        on_target=`awk 'FNR>1' "$in_target" | wc -l`
        off_target=`bedtools intersect -v -a "$intdir"/${onefile/.gz/} -b "$intdir"/${targets_bed/.gz/} | wc -l`
        printf "$onefile\t$on_target\t$off_target\n" >> $statfile
    fi
done

# subset reference regions from bed files
statfile=$prodir/reference_variant.stats
printf "File\tOn_ref\tOff_ref\n" > $statfile

for onefile in ${sample_vcf[@]} $ref_vcf ${gatk_vcf[@]}
do
    in_target=${onefile/.gz/}
    in_target="$intdir"/${in_target/.vcf/_onref.vcf}
    if [ -f "$in_target" ]
    then
        on_target=`awk 'FNR>1' "$in_target" | wc -l`
        off_target=`bedtools intersect -v -a "$intdir"/${onefile/.gz/} -b "$intdir"/${ref_bed/.gz/} | wc -l`
        printf "$onefile\t$on_target\t$off_target\n" >> $statfile
    fi
done

statfile=$prodir/common_variants.stats
in_target=`for onefile in ${sample_vcf[@]}; do echo "$intdir"/${onefile/.vcf.gz/_target.vcf}; done | tr '\n' ' '`
in_target_ref="$intdir"/${ref_vcf/.vcf.gz/_target.vcf}
in_onref=`for onefile in ${sample_vcf[@]}; do echo "$intdir"/${onefile/.vcf.gz/_onref.vcf}; done | tr '\n' ' '`
in_onref_ref="$intdir"/${ref_vcf/.vcf.gz/_onref.vcf}
i=`echo ${in_onref[@]} $in_onref_ref | tr ' ' '\n' | wc -l`

cat `echo ${in_target[@]} $in_target_ref` > $intdir/temp_cat
ontarget=`cut -f1,2,4,5 $intdir/temp_cat | sort | uniq -c | awk -v N=$i '$1==N' | grep -v '#CHROM' | wc -l`
cat `echo ${in_onref[@]} $in_onref_ref` > $intdir/temp_cat
onref=`cut -f1,2,4,5 $intdir/temp_cat | sort | uniq -c | awk -v N=$i '$1==N' | grep -v '#CHROM' | wc -l`
rm $intdir/temp_cat
printf "Region\tCommon_variants\tFiles\n" > $statfile
printf "Target\t$ontarget\t$i\n" >> $statfile
printf "Reference\t$onref\t$i\n" >> $statfile


for onefile in `echo ${sample_vcf[@]} ${gatk_vcf[@]}`
do
        
    infile=$intdir/${onefile/.vcf.gz/_onref.vcf}
    if [ -f "$infile" ]
    then
       onref=`cat $infile $in_onref_ref | cut -f1,2,4,5 \
             | sort | uniq -c | awk '$1==2' | grep -v '#CHROM' | wc -l`
       printf "Target\t$ontarget\t$onefile, reference\n" >> $statfile
    fi

    infile=$intdir/${onefile/.vcf.gz/_target.vcf}
    if [ -f "$infile" ]
    then
        ontarget=`cat $infile $in_target_ref | cut -f1,2,4,5 \
           | sort | uniq -c | awk '$1==2' | grep -v '#CHROM' | wc -l`

        printf "Reference\t$onref\t$onefile, reference\n" >> $statfile
    fi
done



exit 0

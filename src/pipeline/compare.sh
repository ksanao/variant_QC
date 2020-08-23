#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path)\n"

# Print usage: scriptname -options, if script invoked with no command-line args
if [ $# -eq 0 ]
then
	echo -e $USAGE
	exit $ERROR_ARGUMENTS
fi

source "$config"

# create interim unzipped file extracts
for onefile in ${sample_vcf[@]} $ref_vcf $ref_bed $regions_bed $targets_bed
do
    outfile=${onefile/.gz/}
    gunzip -c "$rawdir"/"$onefile" > "$intdir"/"$outfile"
done

# annotate regions and targets bed files
for onefile in $regions_bed $targets_bed
do
   bed_annotation -g "$ref_version" -o "$prodir"/${onefile/.bed.gz/_annot.bed} "$intdir"/${onefile/.gz/}
done

# index vcf files
pushd $rawdir
for onefile in ${sample_vcf[@]} $ref_vcf
do
    tabix -fp vcf $onefile 
done
popd

# subset target regions from bed files
statfile="$prodir"/target_variant.stats
printf "File\tOn_target\tOff_target\n" > $statfile
 
for onefile in ${sample_vcf[@]} $ref_vcf
do
    outfile="$intdir"/${onefile/.vcf.gz/_target.vcf}
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKESAMPLE\n" > $outfile
    bedtools intersect -wa -a "$intdir"/${onefile/.gz/} -b "$intdir"/${targets_bed/.gz/} >> $outfile
    on_target=`awk 'FNR>1' $outfile | wc -l`
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
    outfile="$intdir"/${onefile/.vcf.gz/_onref.vcf}
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKESAMPLE\n" > $outfile
    bedtools intersect -wa -a "$intdir"/${onefile/.gz/} -b "$intdir"/${ref_bed/.gz/} >> $outfile
    on_target=`awk 'FNR>1' $outfile | wc -l`
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

# index bam files
for onefile in ${sample_bam[@]}
do
    samtools index $rawdir/$onefile
done

# fetch genome sequence
pushd $intdir
printf > ref.fasta
for onefile in ${ref_fasta[@]}
do
    wget "$ref_dna_link"/"$onefile"
    gunzip "$onefile"
    cat ${onefile/.gz/} >> ref.fasta
done
samtools faidx ref.fasta
bedtools getfasta -fi ref.fasta -bed regions.bed > regions.fasta
samtools faidx regions.fasta 

wget "$ref_gff_link"/"$ref_gff"
sortedgff="${ref_gff/gff3.gz/sorted.gff3.gz}"
bedtools sort -i "$ref_gff" | bgzip -c > "$sortedgff"
tabix "$sortedgff"

popd

# visualise coverage from bam files
bash /home/src/utils/regions_samplot.sh $config

# visualise in IGV
in_bam=`for onefile in ${sample_bam[@]}; do echo "$rawdir"/"$onefile"; done`
create_report "$intdir"/${targets_bed/.gz/} "$intdir"/ref.fasta --tracks `echo ${in_bam[@]}` "$rawdir"/"$ref_vcf" "$intdir"/"$sortedgff" --output "$prodir"/igv_viewer.html

# bam stats
bash /home/src/utils/bamsats.sh $config

exit 0

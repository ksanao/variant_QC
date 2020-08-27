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

export "PATH=/home/src/utils:$PATH"

# fetch genome sequence
fetch_ref.sh "$config"

# create interim unzipped file extracts
for onefile in ${sample_vcf[@]} $ref_vcf $ref_bed $regions_bed $targets_bed_all
do
    outfile=${onefile/.gz/}
    gunzip -c "$rawdir"/"$onefile" > "$intdir"/"$outfile"
done

# validate golden reference
validate_reference.sh "$config"

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

# subset vcf files
for onefile in ${sample_vcf[@]} $ref_vcf
do
    subset_vcf.sh $config "$intdir"/${onefile/.gz/}
done


# subset and compare vcfs
compare_vcf_stats.sh $config

# index bam files
for onefile in ${sample_bam[@]}
do
    samtools index $rawdir/$onefile
done

# visualise coverage from bam files
bash /home/src/utils/regions_samplot.sh $config

# visualise in IGV
in_bam=`for onefile in ${sample_bam[@]}; do echo "$rawdir"/"$onefile"; done`
sortedgff="${ref_gff/gff3.gz/sorted.gff3.gz}"
create_report "$intdir"/${targets_bed_all/.gz/} "$intdir"/ref.fasta --tracks `echo ${in_bam[@]}` "$rawdir"/"$ref_vcf" "$intdir"/"$sortedgff" --output "$prodir"/igv_viewer.html

# bam stats
bash /home/src/utils/bamsats.sh $config

exit 0

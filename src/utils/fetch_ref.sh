#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path)\n\
If refefence and benchmark files do not exist, then fetch respective files\n"

# Print usage: scriptname -options, if script invoked with wrong command-line args
if [ $# -ne 1 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"

# fetch genome sequence
pushd $intdir

if [ ! -f ref.fasta ]
then
    printf > ref.fasta
    for onefile in ${ref_fasta[@]}
    do
        wget "$ref_dna_link"/"$onefile"
        gunzip "$onefile"
        cat ${onefile/.gz/} >> ref.fasta
        rm ${onefile/.gz/}
    done
    samtools faidx ref.fasta
fi

if [ ! -f ref.fasta ]
then
    bedtools getfasta -fi ref.fasta -bed regions.bed > regions.fasta
    samtools faidx regions.fasta
fi

# fetch annotations
sortedgff="${ref_gff/gff3.gz/sorted.gff3.gz}"
if [ ! -f $sortedgff ]
then
    wget "$ref_gff_link"/"$ref_gff"
    bedtools sort -i "$ref_gff" | bgzip -c > "$sortedgff"
    tabix "$sortedgff"
    rm "$ref_gff"
fi

# fetch clinical vairants vcf
if [ ! -f $ref_clinvar ]
then
    wget "$ref_clinvar_path"/"$ref_clinvar"
    tabix -fp vcf $ref_clinvar
fi

# fetch high confidence variants
if [ ! -f $bqsr/$ref_highconf ]
then
    mkdir -p $bqsr
    cd $bqsr
    wget "$ref_highconf_path"/"$ref_highconf"
    tabix $ref_highconf
    cd ..
fi

popd

exit 0

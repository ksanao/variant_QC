#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Date: 2020.08.22

config=$1

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` config (with path)\n"

# Print usage: scriptname -options, if script invoked with incorrect command-line args
if [ $# -ne 1 ]
then
        echo -e $USAGE
        exit $ERROR_ARGUMENTS
fi

source "$config"

statfile="$prodir"/bam_alignment.stats
printf "File\tRegion\tAligned\tNotaligned\n" > $statfile
metrfile="$prodir"/bam_metric.stats
samstats="$prodir"/samtools.stats
printf "Filename\tAttribute\tValue\n" > $samstats

for onefile in ${sample_bam[@]}
do
   samtools stats "$rawdir"/"$onefile" | grep '^SN' | tr ' ' '_' | awk -v F=$onefile '{print F"\t"$2"\t"$3}' >> $samstats

   tempfile="$intdir"/temp.stats
   picard BamIndexStats -I "$rawdir"/"$onefile" > $tempfile
   notaligned=`grep NoCoordinateCount $tempfile | awk '{print $NF}'`
   printf "$onefile\tNo_region\t\t$notaligned\n" >> $statfile
   
   for oneregion in `cut -f1 "$intdir"/${regions_bed/.gz/} | sort | uniq`
   do
      aligned_r=`awk -v R="$oneregion" '$1==R' $tempfile | awk '{print $5}'`
      unaligned_r=`awk -v R="$oneregion" '$1==R' $tempfile | awk '{print $NF}'`
      printf "$onefile\t$oneregion\t$aligned_r\t$unaligned_r\n" >> $statfile
   done
   rm $tempfile

   tempfile="$intdir"/temp.metrics
   picard CollectDuplicateMetrics -I "$rawdir"/$onefile -M $tempfile
   cat $tempfile | awk 'NF>3' | tail -n2 >> $metrfile
   rm $tempfile
done

sort $metrfile | uniq > temp.metric
mv temp.metric $metrfile

exit 0


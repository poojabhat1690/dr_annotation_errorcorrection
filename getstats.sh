#!/bin/bash


#####

## Date : 3rd July 2018
##Author : Pooja Bhat


## purpose : this script gets the  stats for the number of priming sites of pre-processing of the 3' end annotation pipeline. 

####


mkdir -p  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/stats_stages/

file="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/stats_stages/primingSite_stats.txt"
if [ -f "$file" ]
then
	echo "$file found deleting old file"

	rm "$file"
else
	echo "$file not found."

	touch /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/stats_stages/primingSite_stats.txt
fi



file_merged="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/stats_stages/mergedSite_stats.txt"
if [ -f "$file_merged" ]
then
        echo "$file_merged found deleting old file"

        rm "$file_merged"
else
        echo "$file_merged not found."

        touch "$file_merged"
fi





timepoints=( TP1 TP2 TP3 TP4 TP5 TP6 TP7 TP8 TP9 )





for i in "${timepoints[@]}"

	do

		INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/"$i"//polyAmapping_allTimepoints/n_100_global_a0/

		minus_countsUnique=$(zcat "$INDIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed.gz | wc -l )
		plus_countsUnique=$(zcat "$INDIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed.gz | wc -l )



		sumuniqueSites=$((minus_countsUnique+plus_countsUnique))
		paste <(echo "$i") <(echo "$sumuniqueSites") -d '\t' >> /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/stats_stages/primingSite_stats.txt

	

		minus_merged=$(zcat "$INDIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed.gz | wc -l )
		plus_merged=$(zcat "$INDIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed.gz | wc -l )


		sumMerged=$((plus_merged+minus_merged))

		paste <(echo "$i") <(echo "$sumMerged") -d '\t' >> "$file_merged"


	done








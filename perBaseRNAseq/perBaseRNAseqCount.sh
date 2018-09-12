#!/bin/bash



### gets per base counts and formats
ml bedtools

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW//


### minus strand


#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49979_wt-1-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49979_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49980_wt-1-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49980_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49981_wt-1-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49981_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49982_wt-1-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49982_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49983_wt-2-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49983_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49984_wt-2-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49984_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49985_wt-2-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49985_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49986_wt-2-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49986_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49987_wt-3-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49987_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49988_wt-3-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49988_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49989_wt-3-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49989_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49990_wt-3-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus49990_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50007_wt-1-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50007_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50008_wt-1-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50008_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50009_wt-1-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50009_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50010_wt-1-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50010_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50011_wt-2-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50011_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50012_wt-2-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50012_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50013_wt-2-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50013_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50014_wt-2-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50014_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50015_wt-3-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50015_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50016_wt-3-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50016_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50017_wt-3-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50017_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_minus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50018_wt-3-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_minus50018_single.bed


## plus strand

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49979_wt-1-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49979_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49980_wt-1-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49980_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49981_wt-1-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49981_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49982_wt-1-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49982_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49983_wt-2-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49983_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49984_wt-2-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49984_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49985_wt-2-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49985_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49986_wt-2-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49986_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49987_wt-3-T1-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49987_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49988_wt-3-T2-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49988_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49989_wt-3-T3-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49989_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/49990_wt-3-T4-polyA_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus49990_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50007_wt-1-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50007_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50008_wt-1-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50008_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50009_wt-1-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50009_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50010_wt-1-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50010_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50011_wt-2-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50011_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50012_wt-2-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50012_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50013_wt-2-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50013_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50014_wt-2-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50014_single.bed

#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50015_wt-3-T1-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50015_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50016_wt-3-T2-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50016_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50017_wt-3-T3-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50017_single.bed
#bedtools coverage -d -split -S -a countingWindows_greaterThan5CPM_plus_200nt.bed -b /clustertmp/pooja/polyAribo0_fastq/50018_wt-3-T4-RiboCop_mappedAligned.sortedByCoord.out.bam > perBasee200nts_5cpmCountingWindows_plus50018_single.bed

#################### pasting the counts from each file





paste <(awk '{print $4}' perBasee200nts_5cpmCountingWindows_minus50018_single.bed ) <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49979_single.bed ) <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49980_single.bed )  <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49981_single.bed )   <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49982_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49983_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49984_single.bed )   <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49985_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49986_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49987_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49988_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49989_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus49990_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50007_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50008_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50009_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50010_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50011_single.bed )      <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50012_single.bed )       <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50013_single.bed )       <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50014_single.bed )        <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50015_single.bed )        <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50016_single.bed )         <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50017_single.bed )         <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_minus50018_single.bed ) > counts_perBase_3replicates_minus.bed


paste <(awk '{print $4}' perBasee200nts_5cpmCountingWindows_plus50018_single.bed ) <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49979_single.bed ) <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49980_single.bed ) <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49981_single.bed )   <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49982_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49983_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49984_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49985_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49986_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49987_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49988_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49989_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus49990_single.bed )    <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50007_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50008_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50009_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50010_single.bed )     <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50011_single.bed )      <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50012_single.bed )       <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50013_single.bed )       <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50014_single.bed )        <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50015_single.bed )        <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50016_single.bed )         <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50017_single.bed )         <(awk '{print $8}' perBasee200nts_5cpmCountingWindows_plus50018_single.bed ) >counts_perBase_3replicates_plus.bed







########  get mean of the counts 

awk 'NR==0{next} { printf("%.2f\n", ($2 + $3 + $4 + $5 + $6 + $7 + $8 + $9 + $10 $10 + $11 + $12 + $13 + $14 + $15 + $16 + $17 + $18 + $19 $20 + $21 + $22 + $23 + $24 + $25)) }' counts_perBase_3replicates_minus.bed >meanCounts_minus.bed
awk 'NR==0{next} { printf("%.2f\n", ($2 + $3 + $4 + $5 + $6 + $7 + $8 + $9 + $10 $10 + $11 + $12 + $13 + $14 + $15 + $16 + $17 + $18 + $19 $20 + $21 + $22 + $23 + $24 + $25)) }' counts_perBase_3replicates_plus.bed >meanCounts_plus.bed












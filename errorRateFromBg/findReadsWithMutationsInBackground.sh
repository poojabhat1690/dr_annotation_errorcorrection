#!/bin/bash
ml samtools
mkdir /clustertmp/pooja/backgroundMutations

cd /clustertmp/pooja/backgroundMutations/

#samtools view -b -L /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed  /clustertmp/pooja/allUntreated.bam > /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam


#### getting reads with different number of mutations...removing all multimapping reads
samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:0 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_0_mut.txt 

#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:1 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_1_mut.txt 
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:2 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_2_mut.txt 
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:3 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_3_mut.txt #
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:4 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_4_mut.txt 
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:5 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_5_mut.txt 
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:6 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_6_mut.txt 
#samtools view /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam | grep TC:i:7 | awk ' $5==60 {print $0}' > /clustertmp/pooja/backgroundMutations/TC_7_mut.txt 

#samtools view -H /clustertmp/pooja/backgroundMutations/mappingToUTRs.bam > /clustertmp/pooja/backgroundMutations/header.txt 

#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_1_mut.txt  > /clustertmp/pooja/backgroundMutations/TC_1_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_2_mut.txt  >  /clustertmp/pooja/backgroundMutations/TC_2_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_3_mut.txt > /clustertmp/pooja/backgroundMutations/TC_3_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_4_mut.txt > /clustertmp/pooja/backgroundMutations/TC_4_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_5_mut.txt > /clustertmp/pooja/backgroundMutations/TC_5_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_6_mut.txt > /clustertmp/pooja/backgroundMutations/TC_6_mut.bam
#cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_7_mut.txt > /clustertmp/pooja/backgroundMutations/TC_7_mut.bam
cat /clustertmp/pooja/backgroundMutations/header.txt /clustertmp/pooja/backgroundMutations/TC_0_mut.txt > /clustertmp/pooja/backgroundMutations/TC_0_mut.bam

#samtools view -Sb TC_1_mut.bam > TC_1_mut_1.bam
#samtools view -Sb TC_2_mut.bam > TC_2_mut_1.bam
#samtools view -Sb TC_3_mut.bam > TC_3_mut_1.bam
#samtools view -Sb TC_4_mut.bam > TC_4_mut_1.bam
#samtools view -Sb TC_5_mut.bam > TC_5_mut_1.bam
#samtools view -Sb TC_6_mut.bam > TC_6_mut_1.bam
#samtools view -Sb TC_7_mut.bam > TC_7_mut_1.bam 
samtools view -Sb TC_0_mut.bam > TC_0_mut_1.bam


#samtools index TC_1_mut_1.bam
#samtools index TC_2_mut_1.bam
#samtools index TC_3_mut_1.bam
#samtools index TC_4_mut_1.bam
#samtools index TC_5_mut_1.bam
#samtools index TC_6_mut_1.bam
#samtools index TC_7_mut_1.bam
samtools index TC_0_mut_1.bam


### already getting the RA tags from the 0 mutation file as i cannot read this huge file into R
cut -f22 TC_0_mut.txt | cut -c 6- - | sed 's/,/\t/g' - > TC_0_RAtags.txt

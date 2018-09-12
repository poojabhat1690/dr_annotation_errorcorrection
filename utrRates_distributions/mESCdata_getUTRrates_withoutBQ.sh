#!/bin/bash


module purge
module load python
pip install --user IntervalTree
module load joblib/0.9.4-python2.7.3  
module load pysam
module load R/3.2.2
module load samtools/1.3.1
#module load python/2.7.13-foss-2017a

OUTDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/utrRates_noBq/
mkdir -p "$OUTDIR"

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/filter/

for i in *.bam

do

/groups/ameres/Veronika/bin/slamdunk/bin/alleyoop utrrates -o "$OUTDIR" -r /groups/ameres/bioinformatics/references/mmu/mm10/mmu_mm10_whole_genome.fa -mq 0 -b /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/mESC_subset/output/new/final90percent/allAnnotations.bed -l 38 -t 15 "$i"
########################################/groups/ameres/Veronika/bin/slamdunk/bin/alleyoop utrrates -o "$OUTDIR" -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -mq 0 -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations.bed -l 38 "$i"

done





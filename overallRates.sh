#!/bin/bash


#$ -pe smp 2-24



#/groups/ameres/Veronika/bin/slamdunk/bin/alleyoop utrrates -o "$OUTDIR"/filter/ -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -mq 0 -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed -l 88 "$index"_adapterTrimmed_slamdunk_mapped_filtered.bam



INDIR=/clustertmp/pooja/mapping_stageSpecificConuntingWindows_merged_includingPrimiRNA/filter/

module purge
module load python
pip install --user IntervalTree
module load joblib
module load pysam
module load R/3.2.2
module load samtools/1.3.1


cd "$INDIR"

for i in *.bam


do

/groups/ameres/Veronika/bin/slamdunk/bin/alleyoop rates -o /clustertmp/pooja/mapping_stageSpecificConuntingWindows_merged_includingPrimiRNA/Rates/ -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -mq 27  "$i"


done



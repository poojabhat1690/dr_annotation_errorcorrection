#!/bin/bash
#$ -pe smp 2-48



INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/combinedFiles/
ADAPTERTRIMMED=/clustertmp/pooja/mapping_dr10_06022018/

ml cutadapt

### chckeing if the adapter trimmed files exist by checking one of the files... 

file="/clustertmp/pooja/mapping_dr10_06022018/combinedFile_AACGCC.fastq.gz_adapterTrimmed.fastq"

if [ -f "$file" ]
then
	echo "$file found."
else

	echo "$file not found."

		### if the adapter trimmed files do not exist, then generate them...
		cd "$INDIR"
		for i in *.gz
			
			do
							
				cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -o "$ADAPTERTRIMMED"/"$i"_adapterTrimmed.fastq -m 18 --trim-n "$INDIR"/"$i"

			
			done
	echo "the adapter trimmed files have been generated using cutadapt"
fi





#module purge
#module load python
#pip install --user IntervalTree
#module load joblib
#module load pysam
#module load R/3.2.2
#module load samtools/1.3.1

#OUTDIR=/clustertmp/pooja/mapping_stageSpecificAnnotation_includingMT/
#mkdir -p "$OUTDIR"


#while read line; do 




#	timepoint=`echo "$line" | cut -f5`
#	sample=`echo "$line" | cut -f4`
#	adapterTrimmed=`echo "_adapterTrimmed.fastq"`
#	adapterTrimmedFile=$sample$adapterTrimmed
#	echo "$timepoint"
#	echo "$sample"
#	echo "$adapterTrimmedFile"

#TIMEDIR="$OUTDIR"/"$timepoint"
#mkdir -p "$TIMEDIR"

#/groups/ameres/Veronika/bin/slamdunk/bin/slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$TIMEDIR" -b /clustertmp/pooja/SLAMannotation/dr/"$timepoint"//output/final90percent/allAnnotations.bed  -fb /clustertmp/pooja/SLAMannotation/dr/"$timepoint"//output/final90percent//countingWindows_transcriptionalOutput.bed -a 5 -5 12 -n 100 -mv 0.2 -t 15 -mq 0 -mi 0.95 -m -rl 88 "$ADAPTERTRIMMED"/"$adapterTrimmedFile"

	

#done < "$INDIR"/sampleList_timepoints.txt 






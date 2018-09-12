#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-2:30:00     # 2 minutes
#SBATCH --output=/scratch/pooja/runSLAMdunk
#SBATCH --job-name=numberOfConversions
########SBATCH --array=1-63
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at



module load samtools/1.4-foss-2017a

OUTDIR=/scratch/pooja/differnetNumberOfTC/
mkdir -p "$OUTDIR"

for NTC in 1 2 3 4 5; do 

	var=TC:i:"$NTC"
	
	mkdir -p "$OUTDIR"/"$NTC"/
cd  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/
	for i in *.bam

		do
			samtools view -h  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/"$i" | grep $var - > "$OUTDIR"/"$NTC"/"$i"_TC"$NTC"reads.txt

			samtools view -H  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/"$i" > "$OUTDIR"/"$NTC"/header_"$i".txt


			cat "$OUTDIR"/"$NTC"/"$i"_TC"$NTC"reads.txt "$OUTDIR"/"$NTC"/header_"$i".txt > "$OUTDIR"/"$NTC"/"$i"_"$NTC".sam



		done





done

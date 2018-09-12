#!/usr/bin/env bash
qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
###########SBATCH --cpus-per-task=7
#SBATCH --mem=250G
#SBATCH --time=2-5:30:00     # 2 minutes
#SBATCH --output=/scratch/pooja/runSLAMdunk
#SBATCH --job-name=cleanUpInjectedSamples
########SBATCH --array=1-63
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at


module load r/3.4.1-foss-2017a-x11-20170314



Rscript /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/scripts/highConfidenceTCcalling/simplified_vectorized.R







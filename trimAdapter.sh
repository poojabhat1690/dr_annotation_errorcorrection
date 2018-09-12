#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-8:00:00     # 2 minutes
#SBATCH --output=/scratch/pooja/runSLAMdunk
#SBATCH --job-name=runSLAMdunk
############SBATCH --array=1-87
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at

module load cutadapt/1.9.1-foss-2017a-python-2.7.13

INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/combinedFiles/
ADAPTERTRIMMED=/scratch/pooja/mapping_dr10_06022018/

mkdir -p "$ADAPTERTRIMMED"



file="/scratch/pooja/mapping_dr10_06022018/combinedFile_AACGCC.fastq.gz_adapterTrimmed.fastq"

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





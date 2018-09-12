#!/bin/bash


ml bedtools

bedtools getfasta -s -fi /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -bed /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/refseq/peaks_10_120bps.bed  > /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/refseq//peaks_10_120bps.fa




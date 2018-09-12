
#### get the refSeq and refFlat annotations 

refSeq_unmerged = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
refSeqTranscr = refSeq_unmerged
refSeqTranscr <- unique(refSeqTranscr[,c("V1", "V2", "V3", "V4", "V6","V7")])
colnames(refSeqTranscr) <- c("chr", "start", "end", "gid", "strand", "tid")

peaks = refSeqTranscr
library(ggplot2)
library(reshape)
library(GenomicFeatures)
options(scipen=999)

peaks = peaks[,c(1:6)]
colnames(peaks) = paste0("V",c(1:ncol(peaks)))
peaks$V7 = peaks$V2
peaks$V8 = peaks$V3


peaks_plus = peaks[which(peaks$V5 == "+"),]
peaks_minus = peaks[which(peaks$V5 == "-"),]


peaks_plus[,2] = peaks_plus$V8 - 60
peaks_plus[,3] = peaks_plus$V8 + 60

peaks_minus[,2] = peaks_minus$V7 - 60
peaks_minus[,3] = peaks_minus$V7 + 60


peaks_total = rbind(peaks_minus,peaks_plus)

peaks_total_rearrange =cbind.data.frame(peaks_total[,"V1"],peaks_total[,"V2"],peaks_total[,"V3"],peaks_total[,"V4"],0,peaks_total[,"V5"],peaks_total[,"V6"],peaks_total[,"V7"],peaks_total[,"V8"])
colnames(peaks_total_rearrange) = c("chr","start","end","peak_name","count","strand","transcriptName","start_original","end_original")

write.table(peaks_total_rearrange,"/Volumes/groups/ameres/Pooja/Projects/finalUTRannotationPipeline/testsEachStep/internalPrimingEvents/data/refSeqAnnotation//peaks_10_120bps.bed",sep="\t",quote = F,row.names = F,col.names = F)


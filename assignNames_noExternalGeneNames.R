### some of the genes that do not have an external gene name, so not have any assigned gene... I will assign these by overlapping them with ensembl gene annotations.. 
library(GenomicRanges)


countingWindows_merged = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations.bed",sep="\t",stringsAsFactors = F)
countingWindows_merged_noName  = countingWindows_merged[which(countingWindows_merged$V4 == ""),]
ensemblGenes = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/genesEnsembl.bed",stringsAsFactors = F,sep="\t")
ensemblGenes = ensemblGenes[complete.cases(ensemblGenes),]
colnames(ensemblGenes) = paste0("V",c(1:ncol(ensemblGenes)))
ensemblGenes_GR = makeGRangesFromDataFrame(ensemblGenes,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
countingWindows_merged_GR = makeGRangesFromDataFrame(df = countingWindows_merged_noName,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",starts.in.df.are.0based = T,strand.field = "V6")

overlappingWindows = findOverlaps(query = countingWindows_merged_GR,ensemblGenes_GR)
countingWindows_merged_overlapping = countingWindows_merged_noName[queryHits(overlappingWindows),]
ensemblGenes_GR_overlapping = ensemblGenes[subjectHits(overlappingWindows),]

cw_ensemblGenes = cbind.data.frame(countingWindows_merged_overlapping,ensemblGenes_GR_overlapping)


cw_ensemblGenes = cw_ensemblGenes[!duplicated(cw_ensemblGenes[,c(1:6)]),]
countingWindows_merged[which(countingWindows_merged$V4 == ""),]$V4 <- cw_ensemblGenes[,10]

write.table(countingWindows_merged,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allCountingWindows_assigned.bed",sep="\t",quote = F)



######## also doing this for the high confidence mRNA 3' ends
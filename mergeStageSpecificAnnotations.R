library(dplyr)

##### combining stage specific annotation based on presence in two consecutive stages... 
timepoints = paste0("TP",c(1:9))
directories = paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/",timepoints,"/")
finalAnnotations =paste0(directories,"/final90percent/ends_greater90percent_intergenic_n100.txt")
finalAnnotations_data  = lapply(finalAnnotations,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(finalAnnotations_data) = timepoints

##### creating +/- 5 nt widths for the ends...


combinations_samples = vector("list",8)
names(combinations_samples) = c("TP1_TP2","TP2_TP3","TP3_TP4","TP4_TP5","TP5_TP6","TP6_TP7","TP7_TP8","TP8_TP9")



for(i in 1:length(combinations_samples)){
  a = finalAnnotations_data[[i]]
  b = finalAnnotations_data[[i+1]]
  
  a = a %>% mutate (start_shifted = startPeak-5) %>% mutate(end_shifted = endPeak + 5)
  b = b %>% mutate (start_shifted = startPeak-5) %>% mutate(end_shifted = endPeak + 5)
  
  library(GenomicRanges)
  a_ranges = makeGRangesFromDataFrame(df = a,keep.extra.columns = T,ignore.strand = F,start.field = "start_shifted",end.field = "end_shifted",strand.field = "strand",starts.in.df.are.0based = T,seqnames.field = "chr")
  b_ranges = makeGRangesFromDataFrame(df = b,keep.extra.columns = T,ignore.strand = F,start.field = "start_shifted",end.field = "end_shifted",strand.field = "strand",starts.in.df.are.0based = T,seqnames.field = "chr")
  
  a_bOverlap = findOverlaps(query = a_ranges,subject = b_ranges)
  a_overlapping = a[queryHits(a_bOverlap),]
  b_overlapping = b[subjectHits(a_bOverlap),]
  a_overlapping = a_overlapping[!duplicated(a_overlapping),]
  b_overlapping = b_overlapping[!duplicated(b_overlapping),]
  a_bfinal = rbind.data.frame(a_overlapping,b_overlapping)
  a_bfinal = a_bfinal[!duplicated(a_bfinal),]
  combinations_samples[[i]] = a_bfinal
}

allCombinations = do.call(rbind,combinations_samples)
allCombinations = allCombinations[!duplicated(allCombinations),]
allCombinations$origin = row.names(allCombinations)

dir.create(path = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/")
write.table(allCombinations,"///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/highConfidenceEnds_merged.bed",sep="\t",quote = F,row.names = F)


##### now i want to create counting windows for counting and transcriptional output for this dataset... 


### counting windows for half life calculations
#BOut="//clustertmp/bioinfo/pooja/SLAMannotation/dr/output/"
annotation_custom=allCombinations
annotation_custom = annotation_custom[,c(1:8)]
colnames(annotation_custom) = paste("V",c(1:8),sep="")
ucscDir="///groups/ameres/Pooja/Projects/zebrafishAnnotation//dr10/refSeq_dr10_GRCh38_20170504/"
refFlat <- read.table(paste0(ucscDir, "/refFlat.txt"), stringsAsFactors = F)
refSeqTranscr <- refFlat[grep("NM_", refFlat$V2),]

#add ncRNA
ensemblDir="///groups/ameres/Pooja/Projects/zebrafishAnnotation//dr10/ensembl_dr10_Ensembl_Genes_88/"


allEnsembl = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation//dr10/ensembl_dr10_Ensembl_Genes_88/transcriptStartsAndEnds_all.txt",sep="\t",header=T)

classesToInclude = c("antisense", "bidrectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", "processed_transcript", "sense_intronic", "sense_overlapping","protein_coding")
allEnsembl = allEnsembl[allEnsembl$transcript_biotype %in% classesToInclude,]

### not including these anymore but just including long non-coding RNAs from ensembl

ncRNA <- read.table(paste0(ensemblDir, "/ncRNA_refSeq.txt"), stringsAsFactors=F, header=T)
ncRNA <- refFlat[is.element(refFlat$V2, ncRNA[,"refseq_ncrna"]),]

allEnsembl <- unique(allEnsembl[,c("chromosome_name", "transcript_start", "transcript_end", "external_gene_name", "strand", "ensembl_transcript_id")])
colnames(allEnsembl) <- c("chr", "start", "end", "gid", "strand", "tid")

allEnsembl <- cbind(allEnsembl[,1:4], 0,  allEnsembl[,5:6])

allEnsembl$V8 = "EnsemblOriginal"
colnames(allEnsembl) <- paste0("V", 1:8)



annotation_custom_positive = annotation_custom[which(annotation_custom$V6 == "+"),]
annotation_custom_negative = annotation_custom[which(annotation_custom$V6 == "-"),]
noEnsemblEntries = allEnsembl$V7[which(allEnsembl$V4 == "")]
allEnsembl$V4[which(allEnsembl$V4 == "")]<-noEnsemblEntries
allEnsembl_positive = allEnsembl[which(allEnsembl$V6 == "+"),]
allEnsembl_negative = allEnsembl[which(allEnsembl$V6 == "-"),]

mtTranscripts = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/chrMT.txt",sep="\t",stringsAsFactors = F)
mtTranscripts = mtTranscripts[,c("V1","V2","V3","V7","V4","V6")]
colnames(mtTranscripts) <- c("chr", "start", "end", "gid", "strand", "tid")

mtTranscripts_positive = mtTranscripts[which(mtTranscripts$strand == "+"),]
mtTranscripts_negative = mtTranscripts[which(mtTranscripts$strand == "-"),]


######################## positive strand ###########################
library(GenomicRanges)
annotation_custom_positive = rbind(annotation_custom_positive,allEnsembl_positive)


annotation_custom_positive$V2 = annotation_custom_positive$V3 -250
annotation_custom_positive$V2 = annotation_custom_positive$V2 + 1

annotation_custom_positive_split = split(annotation_custom_positive,f = annotation_custom_positive$V4,drop = T )

positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))

#positive_ranges <- with(annotation_custom_positive,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))


allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )

allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 


reducedToDf = function(reduced){
  reduced <- data.frame(seqnames=seqnames(reduced),
                        starts=start(reduced),
                        ends=end(reduced),
                        names=c(names(reduced)),
                        scores=0,strand = strand(reduced))
  return(reduced)
}

allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))

################## minus strand ############################

annotation_custom_negative = rbind(annotation_custom_negative,allEnsembl_negative)

annotation_custom_negative$V3 = annotation_custom_negative$V2 + 250

## changing to 1 based

annotation_custom_negative$V2 = annotation_custom_negative$V2 + 1

annotation_custom_negative_split = split(annotation_custom_negative,f = annotation_custom_negative$V4,drop = T )

negative_ranges = lapply(annotation_custom_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))


allAnnotations_minus_ranges_reduced = lapply(negative_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 



allAnnotations_minus_ranges_reduced_df = lapply(allAnnotations_minus_ranges_reduced,function(x) reducedToDf(x))


allAnnotations_plus_ranges_reduced_df = do.call(rbind,allAnnotations_plus_ranges_reduced_df)
allAnnotations_minus_ranges_reduced_df = do.call(rbind,allAnnotations_minus_ranges_reduced_df)

allAnnotations = rbind(allAnnotations_plus_ranges_reduced_df,allAnnotations_minus_ranges_reduced_df)

### converting back to 0 based annotations : 

allAnnotations$starts = allAnnotations$starts -1

#valid chromosomes
chromosomes_mm = paste0("chr",c(1:25,"M"))
allAnnotations = allAnnotations[is.element(allAnnotations$seqnames, chromosomes_mm),]
allAnnotations = allAnnotations %>% filter(seqnames != "chrM")

mtTranscripts_rearranged = cbind.data.frame(mtTranscripts$chr,mtTranscripts$start,mtTranscripts$end,mtTranscripts$gid,0,mtTranscripts$strand)
colnames(mtTranscripts_rearranged) = colnames(allAnnotations)
allAnnotations = rbind.data.frame(allAnnotations,mtTranscripts_rearranged)
###### i want to remove the chrM from this (from refSeq annotation) and I want to add the complete list of ensembl MT annotations

##### i also wany to add here info about mir430 clusters that are we manually annotated... 
#mir430sites = read.table("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/pri_mir430.txt",sep="\t",stringsAsFactors = F)
#allAnnotations = rbind.data.frame(allAnnotations,mir430sites)
write.table(allAnnotations,"///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations.bed",sep="\t",quote = F,row.names = F,col.names = F)





#########################################################################
## counting windows for transcriptional output and multimapping - for this we need the ensembl 3' utr annotations, refSeq 3' utr annotations
## and intergenic peaks that pass the 90% threshold
#########################################################################

refSeq = read.table(paste0(ucscDir, "/processed/refSeq_mrna_utrsPresent.bed"),stringsAsFactors = F)
ensembl = read.delim(paste0(ensemblDir, "/proteinCoding_annotatedUTRs.bed"),stringsAsFactors = F,header = F)

allAnnotations <- cbind(allAnnotations, "250CountWindow")
colnames(allAnnotations) <- paste0("V", 1:7)
refSeq_ensembl = rbind(refSeq,ensembl,allAnnotations)
refSeq_ensembl_positive = refSeq_ensembl %>% filter(V6=="+")
refSeq_ensembl_negative = refSeq_ensembl %>% filter(V6=="-")


#annotation_custom = read.table("/clustertmp/pooja/mESCinput/final90percent//ends_greater90percent_intergenic_n100.txt",stringsAsFactors = F,sep="\t",header = T)
#annotation_custom_intergenic = annotation_custom[is.element(annotation_custom$peakKind,"intergenic"),]
#annotation_custom_intergenic = annotation_custom_intergenic[,c(1:7)]
#colnames(annotation_custom_intergenic) = paste0("V",c(1:ncol(annotation_custom_intergenic)))
#annotation_custom_intergenic_positive = annotation_custom_intergenic %>% filter(V6 == "+") %>% mutate(V2 = V3-250)
#annotation_custom_intergenic_negative = annotation_custom_intergenic %>% filter(V6 == "-") %>% mutate(V3 = V2+250)

total_positive = refSeq_ensembl_positive
#total_positive = total_positive[!duplicated(total_positive[,c(1:5)]),]
total_positive$V2 = total_positive$V2 +1
total_negative = refSeq_ensembl_negative
#total_negative = total_negative[!duplicated(total_negative[,c(1:5)]),]
total_negative$V2 = total_negative$V2 + 1
###### now reducing Ovelrapping annotations per gene



### positive strand 

total_positive_split = split(total_positive,f = total_positive$V4,drop = T )

total_positive_split_ranges = lapply(total_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_positive_reduced = lapply(total_positive_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 



total_positive_reduced_df = lapply(total_positive_reduced,function(x) reducedToDf(x))

total_positive_reduced_df = do.call(rbind,total_positive_reduced_df)


####### negative strand 



total_negative_split = split(total_negative,f = total_negative$V4,drop = T )

total_negative_split_ranges = lapply(total_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_negative_reduced = lapply(total_negative_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 



total_negative_reduced_df = lapply(total_negative_reduced,function(x) reducedToDf(x))

total_negative_reduced_df = do.call(rbind,total_negative_reduced_df)

countingWindowsTranscriptionalOutput = rbind(total_positive_reduced_df,total_negative_reduced_df)
#colnames(mir430sites) = colnames(countingWindowsTranscriptionalOutput)
#countingWindowsTranscriptionalOutput = rbind.data.frame(countingWindowsTranscriptionalOutput,mir430sites)
countingWindowsTranscriptionalOutput$starts = countingWindowsTranscriptionalOutput$starts -1


write.table(countingWindowsTranscriptionalOutput,"/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined//countingWindows_transcriptionalOutput.bed",sep="\t",quote = F,row.names = F,col.names = F)




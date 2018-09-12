###### for this analysis i want to only take counting windows that are supported by high confidence mRNA 3' ends



highConfidenceEnds = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/highConfidenceEnds_merged.bed",sep="\t",stringsAsFactors = F,header =T)
countingWindowsMerged = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed",sep="\t",stringsAsFactors = F)

##### i want to first identify the counting windows that are supported by any high confidence ends... So i want to overlap the high confidence ends with counting windows... 


countingWindowsMerged_GR = makeGRangesFromDataFrame(df = countingWindowsMerged,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
highConfidenceEnds = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/highConfidenceEnds_merged.bed",sep="\t",stringsAsFactors = F,header =T)
highConfidenceEnds_GR = makeGRangesFromDataFrame(df = highConfidenceEnds,keep.extra.columns = T,ignore.strand = F,start.field = "startPeak",end.field = "endPeak",strand.field = "strand",starts.in.df.are.0based = T,seqnames.field = "chr")

highConfidenceEndOverlap = findOverlaps(query = highConfidenceEnds_GR,subject = countingWindowsMerged_GR)
overlappingHighConfidence  = highConfidenceEnds[queryHits(highConfidenceEndOverlap),] ### all should overlap 
overlappingCountingWindows = countingWindowsMerged[subjectHits(highConfidenceEndOverlap),]
overlappingCountingWindows = overlappingCountingWindows[!duplicated(overlappingCountingWindows),]
overlappingCountingWindows$id = paste0(overlappingCountingWindows$V1,overlappingCountingWindows$V2,overlappingCountingWindows$V3,overlappingCountingWindows$V4,overlappingCountingWindows$V6)

########### reading in the CPMs of reads... 

allCPMS = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F,header = T)
allCPMS$id = paste0(allCPMS$V1,allCPMS$V2,allCPMS$V3,allCPMS$V4,allCPMS$V6)

overlappingCountingWindows_CPMS = allCPMS[allCPMS$id %in% overlappingCountingWindows$id ,]

#### for diffeerent time points I want to get the samples which are >5 cpm
forEachTimePoint = vector("list",9)
timepoints = paste0("TP",c(1:9))
names(forEachTimePoint) = timepoints

for(i in 1:length(forEachTimePoint)){

TP1_samples  = overlappingCountingWindows_CPMS %>% dplyr:: select(matches(timepoints[i])) %>% dplyr::select(matches('Inj|Unt'))
### adding back the metadata 

metadata_CPMwindows = overlappingCountingWindows[,c(1:6)]
TP1_samples = cbind.data.frame(metadata_CPMwindows, TP1_samples)
nGreaterThat5 = apply(TP1_samples[,c(7:10)],1,function(x) length(which(x>=5)))
TP1_samples$nGreaterThan5 = nGreaterThat5
TP1_samples_aboveThreshold = TP1_samples[which(TP1_samples$nGreaterThan5 >= 2),]
TP1_samples_aboveThreshold_CountingWindows = TP1_samples_aboveThreshold[,c(1:6)]
forEachTimePoint[[i]] = TP1_samples_aboveThreshold_CountingWindows

}

forEachTimePoint_samples = do.call(rbind,forEachTimePoint)

forEachTimePoint_samples = forEachTimePoint_samples[!duplicated(forEachTimePoint_samples),]
write.table(forEachTimePoint_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW/combinedTimpoints_suportedByhighConfidenceEnds.bed",sep="\t",quote = F)

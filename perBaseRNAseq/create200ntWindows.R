### counting  windows for per base RNAseq singal 


options(scipen=999)
countingWindows_greaterThan5 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW/combinedTimpoints_suportedByhighConfidenceEnds.bed",stringsAsFactors = F,sep="\t")

#### cpms



library(dplyr)
countingWindows_greaterThan5_plus = countingWindows_greaterThan5 %>% filter(V6=="+")
countingWindows_greaterThan5_minus = countingWindows_greaterThan5 %>% filter(V6=="-")


countingWindows_greaterThan5_plus$V2 = countingWindows_greaterThan5_plus$V3 - 100
countingWindows_greaterThan5_plus$V3 = countingWindows_greaterThan5_plus$V3 + 100

countingWindows_greaterThan5_minus$V3 = countingWindows_greaterThan5_minus$V2 + 100
countingWindows_greaterThan5_minus$V2 = countingWindows_greaterThan5_minus$V2 - 100

### converting to one based coordinates
# countingWindows_greaterThan5_plus$V2 = countingWindows_greaterThan5_plus$V2 +1
# countingWindows_greaterThan5_minus$V2 = countingWindows_greaterThan5_minus$V2 +1
dir.create("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW")
write.table(countingWindows_greaterThan5_plus,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW/countingWindows_greaterThan5CPM_plus_200nt.bed",sep="\t",quote=F,row.names = F,col.names = F)
write.table(countingWindows_greaterThan5_minus,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW//countingWindows_greaterThan5CPM_minus_200nt.bed",sep="\t",quote=F,row.names = F,col.names = F)


#################################### 



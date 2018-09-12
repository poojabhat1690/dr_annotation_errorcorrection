#### adding primirnas 

priMirna = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/pri_mir430.txt",stringsAsFactors = F,header = F)
allAnnotations = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations.bed",sep="\t",stringsAsFactors = F)
withPrimiRNa = rbind.data.frame(priMirna,allAnnotations)
write.table(withPrimiRNa,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed",sep="\t",quote = F,row.names = F,col.names = F)



toCw = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/countingWindows_transcriptionalOutput.bed",sep="\t",stringsAsFactors = F)
toCw_pri = rbind.data.frame(priMirna,toCw)
##### Also i want to add the mitochrondrial genes only present in ENSEMBL 

write.table(toCw_pri,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/countingWindows_transcriptionalOutput_priMirna.bed",sep="\t",quote = F,row.names = F,col.names = F)

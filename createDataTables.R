library(plyr)


###### making TPMs and count data from SLAMseq samples ... 
countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/count/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis//annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/count//",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))

countData_timeCourse_datasets = countDataSets_quantSeq_data
splitNames = unlist(lapply(strsplit(countDataSets_quantSeq,"_",T),function(x) x[2]))
barcodes = unlist(lapply(strsplit(splitNames,".",T),function(x) x[1]))
names(countData_timeCourse_datasets) = barcodes

### reading in the sample information file

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t")
countData_timeCourse_datasets_names = as.data.frame(names(countData_timeCourse_datasets))
colnames(countData_timeCourse_datasets_names) = "V2"

countData_timeCourse_datasets_names = join(countData_timeCourse_datasets_names,sampleInfo)
names(countData_timeCourse_datasets) = countData_timeCourse_datasets_names$V3

countDataSets_quantSeq_data = countData_timeCourse_datasets


timepoints = paste0("TP",c(1:9))

reads_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadCount))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
reads_samples = reads_samples[,orderAll]
metadata_samples =read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed",sep="\t",stringsAsFactors = F)
reads_samples = cbind.data.frame(metadata_samples,reads_samples)



dir.create(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/")
write.table(reads_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/reads_allCountingWindows.txt",sep="\t",quote = F,row.names = F,col.names = T)



cpm_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_samples = cpm_samples[,orderAll]
cpm_samples = cbind.data.frame(metadata_samples,cpm_samples)
write.table(cpm_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables//RPMData_allCountingWindows.txt",sep="\t",quote = F,row.names = F)

#### cpms of T>C reads 

cpm_TCsamples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ConversionRate))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_TCsamples = cpm_TCsamples[,orderAll]
cpm_TCsamples = cbind.data.frame(metadata_samples,cpm_TCsamples)
write.table(cpm_TCsamples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables//allCountingWindows_TCConversionRate.txt",sep="\t",quote = F,row.names = F)


###### taking the counting windws thrt have 5 CPM in at least 2 stages...


Injection_R1 = cpm_samples[grep("Inj_R1",names(cpm_samples))]
Injection_R2 =  cpm_samples[grep("Inj_R2",names(cpm_samples))]
Injection_R3 =  cpm_samples[grep("Inj_R3",names(cpm_samples))]

Incubation_R1 = cpm_samples[grep("Inc_R1",names(cpm_samples))]
Incubation_R2 =  cpm_samples[grep("Inc_R2",names(cpm_samples))]
Incubation_R3 =  cpm_samples[grep("Inc_R3",names(cpm_samples))]

Untreated = cpm_samples[grep("Untreated",names(cpm_samples))]

#### recording the number of time points that have >5RPM in the different time points ... Injection, Incubation, Untreated samples... 

cpm_samples$inInjected_R1 = NA
cpm_samples$inInjected_R2 = NA
cpm_samples$inInjected_R3 = NA
cpm_samples$inIncubation_R1 = NA
cpm_samples$inIncubation_R2 = NA
cpm_samples$inIncubation_R3 = NA

cpm_samples$inInjected_R1 = apply(Injection_R1,1,function(x) length(which(x>5)))
cpm_samples$inInjected_R2 = apply(Injection_R2,1,function(x) length(which(x>5)))
cpm_samples$inInjected_R3 = apply(Injection_R3,1,function(x) length(which(x>5)))
cpm_samples$inIncubation_R1 = apply(Incubation_R1,1,function(x) length(which(x>5)))
cpm_samples$inIncubation_R2 = apply(Incubation_R2,1,function(x) length(which(x>5)))
cpm_samples$inIncubation_R3 = apply(Incubation_R3,1,function(x) length(which(x>5)))
cpm_samples$untreatedRep = apply(Untreated,1,function(x) length(which(x>5)))

##### now taking all the counting windows that have more that 2 stages above 5 CPM in all replicates of Injected and Incubated 
cpmSamples_Above5CPM = cpm_samples[which(cpm_samples$inInjected_R1>=2 & cpm_samples$inInjected_R2>=2 & cpm_samples$inInjected_R3>=2 & cpm_samples$untreatedRep>=2 ),]
write.table(cpmSamples_Above5CPM,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/CPMs_greaterThan5.txt",sep="\t",col.names = T,quote = F)



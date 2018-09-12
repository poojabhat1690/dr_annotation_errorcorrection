################## the number of genes that have 1,2,3,4 TC mutations at different time points - represents the number of detectable genes at early time points.. 
##### i want to compare this to the mitoghondrial gens that are known to be expressed at early time points.. 
theme_ameres <- function (type) {  ### plotting function
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}

library(Rsamtools)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
#source("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/scripts/highConfidenceTCcalling/functions_highConfidenceTCcalling.R")

### from the RPMS of the data 

rpms_allCountingWindows = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows = rpms_allCountingWindows %>% select(1:6) %>% mutate(id = paste(V1,V2,sep="_"))

#### reading in all the reads for different files... 
zygoticReads_samples = list.files(path = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/",pattern = "freqAllMutations_allTranscripts*")
#zygoticReads_samples = zygoticReads_samples[28]
zygoticReads_samples_path = paste0("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/",zygoticReads_samples)

zygoticReads_samples_data = lapply(zygoticReads_samples_path,function(x) read.table(x,stringsAsFactors = F))
names(zygoticReads_samples_data) = zygoticReads_samples

for(i in 1:length(zygoticReads_samples_data)){
  zygoticReads_samples_data_current = zygoticReads_samples_data[[i]]
  zygoticReads_samples_data_current =zygoticReads_samples_data_current %>% mutate(id_samples = str_extract(string = zygoticReads_samples_data[[i]]$id,pattern = "[^-]+"))
  
  
  zygoticReads_samples_data_current = zygoticReads_samples_data_current %>% group_by(id_samples ) %>% mutate (n = n()) %>% ungroup()
  minusTab = zygoticReads_samples_data_current %>% dplyr::filter(strand == "-") %>%
    dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,n,id_samples)
  plusTab = zygoticReads_samples_data_current %>% dplyr::filter(strand == "+") %>%
    dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,n,id_samples)
  
  totalMutations = bind_rows(plusTab,minusTab)
  
  
  A_Cmutns = totalMutations %>% group_by(id_samples) %>% select(A_C) %>% table()
  A_Gmutns = totalMutations %>% group_by(id_samples) %>% select(A_G) %>% table()
  A_Tmutns = totalMutations %>% group_by(id_samples) %>% select(A_T) %>% table()
  C_Amutns = totalMutations %>% group_by(id_samples) %>% select(C_A) %>% table()
  C_Gmutns = totalMutations %>% group_by(id_samples) %>% select(C_G) %>% table()
  C_Tmutns = totalMutations %>% group_by(id_samples) %>% select(C_T) %>% table()
  G_Amutns = totalMutations %>% group_by(id_samples) %>% select(G_A) %>% table()
  G_Cmutns = totalMutations %>% group_by(id_samples) %>% select(G_C) %>% table()
  G_Tmutns = totalMutations %>% group_by(id_samples) %>% select(G_T) %>% table()
  T_Amutns = totalMutations %>% group_by(id_samples) %>% select(T_A) %>% table()
  T_Cmutns = totalMutations %>% group_by(id_samples) %>% select(T_C) %>% table()
  T_Gmutns = totalMutations %>% group_by(id_samples) %>% select(T_G) %>% table()
  
  
  allMutations = list(A_Cmutns,A_Gmutns,A_Tmutns,C_Amutns,C_Gmutns,C_Tmutns,G_Amutns,G_Cmutns,G_Tmutns,T_Amutns,T_Cmutns,T_Gmutns)
  names(allMutations) = c("A_C","A_G","A_T","C_A","C_G","C_T","G_A","G_C","G_T","T_A","T_C","T_G")
  dir.create(path = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/")
  
  save(allMutations,file =paste0("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/" ,zygoticReads_samples[i],".Rdata"))
  
}


allRdata = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",pattern = "*.Rdata")
allRdata_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",allRdata)
names_samples = str_extract(allRdata_path,"([^/]+)$" ) %>% str_sub(start = -20)
 
allSamples_mutations = vector("list",length(names_samples))
names(allSamples_mutations) = names_samples



for(i in 1:length(allRdata_path)){
  load(allRdata_path[i])
  allSamples_mutations[[i]] = as.data.frame(allMutations$T_C) %>% mutate(T_Cnum = as.numeric(as.character(T_C)))  %>% filter(T_Cnum >1) %>% filter(Freq > 0)

}
 
 
 allSamples_mutations_allFiles = do.call(rbind.data.frame,allSamples_mutations)
# 
 rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
 
# #### i also need to add the information of th genes to this... 
rpms_allCountingWindows$id = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1 )
### i also want to subset the genes that are expressed at t1 and t2
TP1 = allSamples_mutations_allFiles[grep("TP1",rownames(allSamples_mutations_allFiles)),]
TP2 = allSamples_mutations_allFiles[grep("TP2",rownames(allSamples_mutations_allFiles)),]
TP1_TP2 = rbind(TP2)
names_genes_expression  = rpms_allCountingWindows[rpms_allCountingWindows$id %in% TP1_TP2$id_samples,]
unique(names_genes_expression$V4)
write.table(allSamples_mutations_allFiles,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/moreThan1mutations_TP2allTranscripts.txt")


#####################################################################



##### i want to know a cumuklative sum of the number of genes that have >=2TCs at each point... 


TP1 = do.call(rbind,allSamples_mutations[grep("TP1",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP2 = do.call(rbind,allSamples_mutations[grep("TP2",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP3 = do.call(rbind,allSamples_mutations[grep("TP3",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP4 = do.call(rbind,allSamples_mutations[grep("TP4",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP5 = do.call(rbind,allSamples_mutations[grep("TP5",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP6 = do.call(rbind,allSamples_mutations[grep("TP6",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP7 = do.call(rbind,allSamples_mutations[grep("TP7",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP8 = do.call(rbind,allSamples_mutations[grep("TP8",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)
TP9 = do.call(rbind,allSamples_mutations[grep("TP9",names(allSamples_mutations))]) %>% filter(.,!grepl("chrM",id_samples)) %>% distinct(id_samples)

allTPs = list(TP1,TP2,TP3,TP4,TP5,TP6,TP7,TP8,TP9)
names(allTPs) = c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9")
allTPs_melt = melt(lapply(allTPs,nrow)) %>% mutate(cumsum_val = cumsum(value)/max(cumsum(value))) 




##### only mitochondrial transcripts 



TP1 = do.call(rbind,allSamples_mutations[grep("TP1",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP2 = do.call(rbind,allSamples_mutations[grep("TP2",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP3 = do.call(rbind,allSamples_mutations[grep("TP3",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP4 = do.call(rbind,allSamples_mutations[grep("TP4",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP5 = do.call(rbind,allSamples_mutations[grep("TP5",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP6 = do.call(rbind,allSamples_mutations[grep("TP6",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP7 = do.call(rbind,allSamples_mutations[grep("TP7",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP8 = do.call(rbind,allSamples_mutations[grep("TP8",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
TP9 = do.call(rbind,allSamples_mutations[grep("TP9",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)

allTPs = list(TP1,TP2,TP3,TP4,TP5,TP6,TP7,TP8,TP9)
names(allTPs) = c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9")
allTPs_melt_mt = melt(lapply(allTPs,nrow)) %>% mutate(cumsum_val = cumsum(value)/max(cumsum(value))) 


plotNumberOftranscripts = list(allTPs_melt,allTPs_melt_mt)
names(plotNumberOftranscripts) = c("zygoticTranscripts","mt transcripts")
plotNumberOftranscripts = do.call(rbind,plotNumberOftranscripts)
plotNumberOftranscripts$type = c(rep("allTranscripts",9),rep("mt",9))
plotNumberOftranscripts$TP = paste0("TP",c(1:9))


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/NumberOfGenes_morethan1TC_allCws.pdf",height=5)
ggpubr::ggline(data = plotNumberOftranscripts,x = 'TP',y='value',group = 'type',col='type',size = 1) + ylab("Number of CW with >=2 TCs")  + theme_ameres(type = "barplot") + scale_color_brewer(palette = "Dark2")
dev.off()



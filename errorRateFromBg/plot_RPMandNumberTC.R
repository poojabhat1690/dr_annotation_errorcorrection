########## steady state gene expression versus number of reads with TCs... 
library(Rsamtools)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)



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

  ####### i want to plot this for the background and for the treatment samples... 
      #### load all the untreated samples... 
### i want to check if in the untreated samples, there is a correlation between the number of reads with 3TCs and Steady state gene expression      


allFiles = c("freqAllMutations_allTranscriptsCAGCGTUntreated_TP1.txt.Rdata","freqAllMutations_allTranscriptsGATCACUntreated_TP2.txt.Rdata","freqAllMutations_allTranscriptsACCAGTUntreated_TP3.txt.Rdata","freqAllMutations_allTranscriptsTGCACGUntreated_TP4.txt.Rdata", "freqAllMutations_allTranscriptsACATTAUntreated_TP5.txt.Rdata", "freqAllMutations_allTranscriptsGTGTAGUntreated_TP6.txt.Rdata","freqAllMutations_allTranscriptsCTAGTCUntreated_TP7.txt.Rdata", "freqAllMutations_allTranscriptsTGTGCAUntreated_TP8.txt.Rdata","freqAllMutations_allTranscriptsTCAGGAUntreated_TP9.txt.Rdata" )
allFiles_paths = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",allFiles)

rpms_allCountingWindows = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/reads_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows$id = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/scatter_expression_readsWith3TCs.pdf")

    for(i in 1:length(allFiles)){
      load(allFiles_paths[i])
      refSample = as.data.frame.matrix(allMutations$T_C) %>% tibble::rownames_to_column(var = "id") %>% filter(`3` > 0)
      samples_combined = plyr::join(refSample,rpms_allCountingWindows)
      samples_combined_untreated = samples_combined %>% select(contains("Unt"))
      p = ggplot(data = samples_combined,aes(x = samples_combined$`3`,y=log10(samples_combined_untreated[,i]))) + geom_point(alpha=.2) + theme_ameres(type="barplot") 
      p = p + ylab(paste0("reads timepoint",i)) + xlab("Number of reads with 3 TCs") + ggtitle(paste0("Untreated samples ",i))
      print(p)
    }

dev.off()



###### i want to check the same for injected samples...

allFiles = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",pattern = "freqAllMutations_allTranscripts*")

allFiles=  allFiles[grep("Inj",allFiles)]
allFiles_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",allFiles)

load(allFiles_path[[5]])
df_considering = as.data.frame.matrix(allMutations$T_C) 
df_considering$greaterTHa2 = apply(df_considering[,4:ncol(df_considering)],1,function(x) sum(x))
df_considering = df_considering %>% rownames_to_column(var = "id") %>% filter(greaterTHa2>1)



# ord = order(apply(allMutations$T_C,1,sum))
# allMutations$T_C = allMutations$T_C[ord,]
# numberOfreads_TC = as.data.frame.matrix(allMutations$T_C) %>% mutate(numReads = apply(.,1,sum)) %>% mutate(log10_numReads = log10(numReads)) %>%mutate(fractionOf3TCreads = `3`/numReads)
# numberOfreads_TC %>% filter(numReads == 1) %>% filter(`3` ==1)
# plot(numberOfreads_TC$log10_numReads,numberOfreads_TC$`3`)
# #### so probability of finding reads >= 3 TC = 1-(prob(0) + prob(1) + prob(2))
# 
# prob0= sum(numberOfreads_TC$`0`)/sum(numberOfreads_TC$numReads)
# prob1 = sum(numberOfreads_TC$`1`)/sum(numberOfreads_TC$numReads)
# prob2 = sum(numberOfreads_TC$`2`)/sum(numberOfreads_TC$numReads)
# probOf_greaterTHan2 = 1 - (prob0 + prob1 + prob2)



comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}
p =  probOf_greaterTHan2

a= 1
b=100
(p^a) * (1-p)^(b-a)

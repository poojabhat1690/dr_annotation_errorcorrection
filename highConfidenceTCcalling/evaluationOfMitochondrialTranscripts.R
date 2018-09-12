###### evaluation of SLAMseq signal at mitochondrial transcripts... 

##### I want to evaluate the following: 

      ###### steady state gene expression at mitochondrial transcripts. 
      ###### incorporation rates at mutochondrial transcripts 
      ###### Number of transcripts that have 2 TC reads, 3 TC reads, 4 TC reads at different time points. 
      ###### per transcript how many TC conversions are there per read
###################################################################################################################


##### steady state expression of mitochondrial transcripts... 


library(reshape)
theme_ameres <- function (type) { 
  
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



##################




RPMs_all = read.table("///Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F,header=T)
conversionRates_all = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017//analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables//allCountingWindows_TCConversionRate.txt",sep="\t",stringsAsFactors = F,header=T)

#### i want to read in the metadata from the annotation used... 

annotation = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed",sep="\t")


###### only keeping the injection samples 

library(dplyr)


splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit %>% select_if(grepl(condition,names(.)))
  dataFrameToSplit_condition_R1 = dataFrameToSplit_condition %>% select_if(grepl("R1",names(.)))
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition %>% select_if(grepl("R2",names(.)))
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition %>% select_if(grepl("R3",names(.)))
  mean_repl = (dataFrameToSplit_condition_R1+dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/3
  dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R1,dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R1","R2","R3","mean")
  return(splitReplicates)
}

conversionRates_injection_replicates = splitReplicates(dataFrameToSplit = conversionRates_all,condition = "Inj",metadata_add = annotation )
RPM_injection_replicates = splitReplicates(dataFrameToSplit = RPMs_all,condition = "Inj",metadata_add = annotation )

#### get the injected mean and and background .. 

untreatedTCrates =  conversionRates_all %>% select_if(grepl("Unt",names(.)))
injectedTC = conversionRates_injection_replicates$mean[,c(1:9)]

backgroundSubtracted_zygotic = injectedTC - untreatedTCrates
backgroundSubtracted_zygotic = as.matrix(backgroundSubtracted_zygotic)
backgroundSubtracted_zygotic[which(backgroundSubtracted_zygotic<0)]<-0
backgroundSubtracted_zygotic = as.data.frame(backgroundSubtracted_zygotic)

backgroundSubtracted_zygotic = cbind.data.frame(conversionRates_all[,c(1:6)],backgroundSubtracted_zygotic)
mtGenes = backgroundSubtracted_zygotic[grep("chrM",backgroundSubtracted_zygotic$V1),]

mtGenes_melt = melt(mtGenes[,c(7:15)])
mtGenes_melt$gene = mtGenes$V4
palette_custom = c(RColorBrewer::brewer.pal(n = 8,"Dark2"),RColorBrewer::brewer.pal(n = 5,"Set1"))
ggplot(mtGenes_melt,aes(x=variable,y=value,group=gene,col=gene)) + geom_line() + scale_color_manual(values = palette_custom)
p = ggpubr::ggline(data = mtGenes_melt,x = 'variable',y = 'value',group = 'gene',col='gene')  + scale_color_manual(values = palette_custom) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate")
dir.create(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/") 
library(ggplot2)
ggsave(filename = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/mitochondrialRNArates.pdf",p,height = 5,width = 5)

### plotting the gene expression of mt genes 

meanExpression_ChrM = RPM_injection_replicates$mean[grep("chrM",RPM_injection_replicates$mean$V1),]
row.names(meanExpression_ChrM) = meanExpression_ChrM$V4

normalizedMeanExpresision =t(apply(meanExpression_ChrM[,c(1:9)],1,function(x) x/max(x) ))
normalizedMeanExpresision_melt = melt(normalizedMeanExpresision)
p = ggpubr::ggline(data = normalizedMeanExpresision_melt,x = 'X2',y = 'value',group = 'X1',col='X1')  + scale_color_manual(values = palette_custom) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("Norm. expression") + ylim(c(0,1))
ggsave(filename = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/mitochondrialGeneExpression.pdf",p,height = 5,width = 5)



############################## i would also like to calculate the conversion rates of these transcripts... 

utrRates_samples = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/utrRates/","*.csv")
utrRates_samples_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/utrRates/",utrRates_samples)
allUTRRates = lapply(utrRates_samples_path,function(x) read.table(x,stringsAsFactors = F,header = T,sep="\t"))

splitNames = unlist(lapply(strsplit(utrRates_samples,"_",T),function(x) x[2]))
barcodes = unlist(lapply(strsplit(splitNames,".",T),function(x) x[1]))
names(allUTRRates) = barcodes

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t")
allUTRRates_names = as.data.frame(names(allUTRRates))
colnames(allUTRRates_names) = "V2"

allUTRRates_names = join(allUTRRates_names,sampleInfo)
names(allUTRRates) = allUTRRates_names$V3
allUTRRates_mt = lapply(allUTRRates,function(x) x[which(x$Chr == "chrM"),])
allUTRRates_Inj_R1 = allUTRRates_mt[grep("Inj_R1",names(allUTRRates_mt))]
allUTRRates_Inj_R2 = allUTRRates_mt[grep("Inj_R2",names(allUTRRates_mt))]
allUTRRates_Inj_R3 = allUTRRates_mt[grep("Inj_R3",names(allUTRRates_mt))]

allUTRRates_Untreated = allUTRRates_mt[grep("Unt",names(allUTRRates_mt))]

getTCrates = function(dataFrame){
  
  
  
  T_C_all = vector("list",length(dataFrame))
  names(T_C_all) = names(dataFrame)
  
  
  
  
  for(i in 1:length(dataFrame))
  {
    
    
    curTab = dataFrame[[i]]
    plusTab = curTab %>% dplyr::filter(Strand == "+")
    minusTab = curTab %>% dplyr::filter(Strand == "-") %>%
      dplyr::select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)
    
    plusTab = plusTab %>% dplyr::select(-contains("N")) %>%
      dplyr::select(-one_of("Chr","Start")) %>%
      filter(rowSums(.) > 0) %>%
      mutate(Asum = A_A + A_C + A_G + A_T) %>%
      mutate(Csum = C_A + C_C + C_G + C_T) %>% 
      mutate(Gsum = G_A + G_C + G_G + G_T) %>%
      mutate(Tsum = T_A + T_C + T_G + T_T) %>%
      mutate_each(funs(. / Asum),matches("A_")) %>%
      mutate_each(funs(. / Csum),matches("C_")) %>%
      mutate_each(funs(. / Gsum),matches("G_")) %>%
      mutate_each(funs(. / Tsum),matches("T_")) %>%
      dplyr::select(-one_of(c("Asum","Csum","Gsum","Tsum"))) %>%
      mutate_each(funs(. * 100))
    
    
    minusTab = minusTab %>% filter(rowSums(.) > 0) %>%
      mutate(Asum = A_A + A_C + A_G + A_T) %>%
      mutate(Csum = C_A + C_C + C_G + C_T) %>% 
      mutate(Gsum = G_A + G_C + G_G + G_T) %>%
      mutate(Tsum = T_A + T_C + T_G + T_T) %>%
      mutate_each(funs(. / Asum),matches("A_")) %>%
      mutate_each(funs(. / Csum),matches("C_")) %>%
      mutate_each(funs(. / Gsum),matches("G_")) %>%
      mutate_each(funs(. / Tsum),matches("T_")) %>%
      dplyr::select(-one_of(c("Asum","Csum","Gsum","Tsum"))) %>%
      mutate_each(funs(. * 100))
    
    plotTab = rbind(plusTab, minusTab)
    plotTab = plotTab %>% dplyr::select(A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G)
    quantiles = lapply(plotTab, function(x) {
      return(quantile(x, na.rm=TRUE, p=0.75) + 1.5 * IQR(x, na.rm=TRUE))
    })
    
    ymax = ceiling(max(unlist(quantiles)))
    plotTab = gather(plotTab,"class","values",A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G)
    
    plotTab$highlight = "no"
    plotTab$highlight[plotTab$class == "T_C"] = "yes"
    T_C = plotTab[which(plotTab$class == "T_C"),]
    T_C_all[[i]] = T_C
    
  }
  
  return(T_C_all)
}
allUTRRates_Inj_R3 = list(allUTRRates_Inj_R3$Inj_R3_TP1,allUTRRates_Inj_R3$Inj_R3_TP2,allUTRRates_Inj_R3$Inj_R3_TP3,allUTRRates_Inj_R3$Inj_R3_TP4,allUTRRates_Inj_R3$Inj_R3_TP5,allUTRRates_Inj_R3$Inj_R3_TP6,allUTRRates_Inj_R3$Inj_R3_TP7,allUTRRates_Inj_R3$Inj_R3_TP8,allUTRRates_Inj_R3$Inj_R3_TP9)
allUTRRates_Inj_R2 = list(allUTRRates_Inj_R2$Inj_R2_TP1,allUTRRates_Inj_R2$Inj_R2_TP2,allUTRRates_Inj_R2$Inj_R2_TP3,allUTRRates_Inj_R2$Inj_R2_TP4,allUTRRates_Inj_R2$Inj_R2_TP5,allUTRRates_Inj_R2$Inj_R2_TP6,allUTRRates_Inj_R2$Inj_R2_TP7,allUTRRates_Inj_R2$Inj_R2_TP8,allUTRRates_Inj_R2$Inj_R2_TP9)
allUTRRates_Inj_R1 = list(allUTRRates_Inj_R1$Inj_R1_TP1,allUTRRates_Inj_R1$Inj_R1_TP2,allUTRRates_Inj_R1$Inj_R1_TP3,allUTRRates_Inj_R1$Inj_R1_TP4,allUTRRates_Inj_R1$Inj_R1_TP5,allUTRRates_Inj_R1$Inj_R1_TP6,allUTRRates_Inj_R1$Inj_R1_TP7,allUTRRates_Inj_R1$Inj_R1_TP8,allUTRRates_Inj_R1$Inj_R1_TP9)
allUTRRates_Untreated = list(allUTRRates_Untreated$Untreated_TP1,allUTRRates_Untreated$Untreated_TP2,allUTRRates_Untreated$Untreated_TP3,allUTRRates_Untreated$Untreated_TP4,allUTRRates_Untreated$Untreated_TP5,allUTRRates_Untreated$Untreated_TP6,allUTRRates_Untreated$Untreated_TP7,allUTRRates_Untreated$Untreated_TP8,allUTRRates_Untreated$Untreated_TP9)

rates_Inj3  = getTCrates(allUTRRates_Inj_R3)
rates_Inj2  = getTCrates(allUTRRates_Inj_R2)
rates_Inj1  = getTCrates(allUTRRates_Inj_R1)
rates_Untreated  = getTCrates(allUTRRates_Untreated)





library(reshape)
library(ggplot2)
plottingFunction = function(data,condition){
  T_C_all_melt  = melt(data)
  #destination = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/initialEvaluationOfData//plots/T_Crates_utr_",condition,".pdf")
 # pdf(destination,height=5,width=9)
  p = ggplot2::ggplot(T_C_all_melt,aes(x=L1,y=value,group=L1)) + geom_boxplot(fill='red') + ggtitle(condition) + theme_bw() + theme_ameres(type = "barplot")
  p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) 
  p = p + xlab("Timepoint") + ylab("Mutation rate per UTR base [%]")
  return(p)
  #print(p)
  #dev.off()
  
}

dir.create("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/onlyMitochondrial/")

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/onlyMitochondrial/TCratesPerMt.pdf",height=4,width=5)

plottingFunction(data = rates_Inj1,condition = "Injection_R1")
plottingFunction(data = rates_Inj2,condition = "Injection_R2")
plottingFunction(data = rates_Inj3,condition = "Injection_R3")
plottingFunction(data = rates_Untreated,condition = "Untreated")

dev.off()

######################################################################################

##### so this rate seems lower than those obtained for purely zygotic transcripts... 
##### this is because there are still unlabelled maternal transcripts.. 



#############  we can look at the number fo transcripts that have 1 TC, 2 TC,3TC... 


rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows$id = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1 )





getGenes = function(allRdata_path,nMut){
  allSamples_mutations = vector("list",length(names_samples))
  names(allSamples_mutations) = names_samples
  
  for(i in 1:length(allRdata_path)){
    load(allRdata_path[i])
    allSamples_mutations[[i]] = as.data.frame(allMutations$T_C) %>% mutate(T_Cnum = as.numeric(as.character(T_C)))  %>% filter(T_Cnum >nMut) %>% filter(Freq > 0)
    
  }
  allSamples_mutations_allFiles =  do.call(rbind.data.frame,allSamples_mutations)
  TP1 = allSamples_mutations_allFiles[grep("TP1",rownames(allSamples_mutations_allFiles)),]
  TP2 = allSamples_mutations_allFiles[grep("TP2",rownames(allSamples_mutations_allFiles)),]
  TP1_TP2 = rbind(TP2)
  names_genes_expression  = rpms_allCountingWindows[rpms_allCountingWindows$id %in% TP1_TP2$id_samples,]
  unique(names_genes_expression$V4)
  #write.table(names_genes_expression,paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized///numberOfGenes_allTranscripts//moreThan",nMut,"mutationsTP2_zygoticTranscripts.txt"),sep="\t",quote = F)
  
  
  
  ##### only mitochondrial transcripts 
  
  
  
  TP1 = do.call(rbind,allSamples_mutations[grep("TP1",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples))  %>% distinct(id_samples)
  TP2 = do.call(rbind,allSamples_mutations[grep("TP2",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples))  %>% distinct(id_samples)
  TP3 = do.call(rbind,allSamples_mutations[grep("TP3",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples))  %>% distinct(id_samples)
  TP4 = do.call(rbind,allSamples_mutations[grep("TP4",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples))  %>% distinct(id_samples)
  TP5 = do.call(rbind,allSamples_mutations[grep("TP5",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
  TP6 = do.call(rbind,allSamples_mutations[grep("TP6",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
  TP7 = do.call(rbind,allSamples_mutations[grep("TP7",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
  TP8 = do.call(rbind,allSamples_mutations[grep("TP8",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
  TP9 = do.call(rbind,allSamples_mutations[grep("TP9",names(allSamples_mutations))]) %>% filter(.,grepl("chrM",id_samples)) %>% distinct(id_samples)
  
  allTPs = list(TP1,TP2,TP3,TP4,TP5,TP6,TP7,TP8,TP9)
  names(allTPs) = c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9")
  allTPs_melt_mt = melt(lapply(allTPs,nrow)) %>% mutate(cumsum_val = cumsum(value)/max(cumsum(value))) 
  
  
  return(allTPs_melt_mt)
  
  
}



library(dplyr)
library(stringr)
allRdata = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized///numberOfGenes_allTranscripts//",pattern = "*.Rdata")
allRdata_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized///numberOfGenes_allTranscripts//",allRdata)
allRdata_path = allRdata_path[grep("Inj",allRdata_path)]
names_samples = str_extract(allRdata_path,"([^/]+)$" ) %>% str_sub(start = -20)

numberOfGenesWithNTCs = vector("list",7)
names(numberOfGenesWithNTCs) = paste0("nMut_",c(1:7))

for(i in 1:7){
  numberOfGenesWithNTCs[[i]]  = getGenes(allRdata_path = allRdata_path,nMut = i )
  
}


### i want to now plot the number of mitochondrial transcripts that have > 1:7 TC conversions.., 

p = melt(numberOfGenesWithNTCs) %>% filter(variable == "value") %>%  
    mutate(time = rep(c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),7)) %>% ggpubr::ggline(data = .,x = 'time',y='value',group='L1',col='L1',size=1) + theme_ameres(type = "barplot") + scale_color_brewer(palette = "Dark2") +
        ylab("Number of mt CWs") + xlab("Time(hpf)")
ggsave(filename = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/onlyMitochondrial/numberOfCWs_numberOfMutations.pdf",plot = p,height = 4,width=5)



###################################################### what are the fraction of reads that have a TC conversion, 2 TC conversions etc... 

        ###### i.e what is the probability of finding a read with 1,2,3,4,5 TC conversion.... 

freqOfTranscripts = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes/",pattern = "freqAllMutations_zygoticTranscripts*")
freqOfTranscripts_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes/",freqOfTranscripts)
numberOfTCreads = vector("list",length(freqOfTranscripts_path))
names(numberOfTCreads) = freqOfTranscripts
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/onlyMitochondrial/numberOfTCs_fractionOfreads.pdf")
for(i in 1:length(freqOfTranscripts_path)){
  load(freqOfTranscripts_path[i])
  chrM_TC   = as.data.frame(allMutations$T_C[grep("chrM",row.names(allMutations$T_C)),] ) %>% 
    group_by(id_samples) %>% mutate(totReads = sum(Freq)) %>% mutate(fracTC = Freq/totReads)   %>%
      ggpubr::ggboxplot(data = .,x = 'T_C',y = 'fracTC',group='id_samples',fill='red') + xlab("number of TC per read") + ylab("Fraction of total reads") + theme_ameres(type = "barplot")+  stat_summary(fun.data = function(x) data.frame(y=1, label = round(median(x),2)), geom="text") + 
          theme(legend.position="none") + ggtitle(label = freqOfTranscripts[i] )
  print(chrM_TC)
}
dev.off()


#### most of the transcripts are unlabelled as there is a maternal contribution..that actually hangs around... so i would expect that
##### the number of reads with multiple TC reads is less in the mitochondrial transcripts.. 






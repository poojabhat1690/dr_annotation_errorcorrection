library(reshape)
theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15,angle = 90, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}


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



####### now I also want to check for genes that do not show transcription in the first stage but show transcription in the second time point (64 cell stage)

meanRPM = RPM_injection_replicates$mean


expressedFromTP2  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 > 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP3  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP4  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP5  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP6  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 == 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP7  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 == 0 & meanRPM$Inj_R1_TP6 == 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]


TC_expressedFromTP2 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 > 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP3 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP4 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP5 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]

TC_expressedFromTP2_melt = melt(TC_expressedFromTP2[,c(7:15)])
TC_expressedFromTP2_melt$gene = TC_expressedFromTP2$V4
#ggplot(TC_expressedFromTP2_melt,aes(x=variable,y=value,group=gene,col=gene)) + geom_line() + scale_color_manual(values = palette_custom) 
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcrption_differentStages.pdf")

ggpubr::ggboxplot(data = TC_expressedFromTP2_melt,x = 'variable',y = 'value',size=1) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate") + geom_text(aes(x=7,y=0.2,label = nrow(TC_expressedFromTP2_melt)/9))  #+ scale_color_manual(values = palette_custom) 


TC_expressedFromTP3_melt = melt(TC_expressedFromTP3[,c(7:15)])
TC_expressedFromTP3_melt$gene = TC_expressedFromTP3$V4
ggpubr::ggboxplot(data = TC_expressedFromTP3_melt,x = 'variable',y = 'value',size=1) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate") + geom_text(aes(x=7,y=0.2,label = nrow(TC_expressedFromTP3_melt)/9))  #+ scale_color_manual(values = palette_custom) 


TC_expressedFromTP4_melt = melt(TC_expressedFromTP4[,c(7:15)])
TC_expressedFromTP4_melt$gene = TC_expressedFromTP4$V4
ggpubr::ggboxplot(data = TC_expressedFromTP4_melt,x = 'variable',y = 'value',size=1) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate") + geom_text(aes(x=7,y=0.2,label = nrow(TC_expressedFromTP4_melt)/9))  #+ scale_color_manual(values = palette_custom) 


############### only for highly expressed














meanRPM = RPM_injection_replicates$mean
backgroundSubtracted_zygotic = backgroundSubtracted_zygotic[which(apply(meanRPM[,c(1:9)],1,max)>5),]
meanRPM = meanRPM[which(apply(meanRPM[,c(1:9)],1,max)>5),]


expressedFromTP2  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 > 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP3  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP4  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP5  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP6  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 == 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
expressedFromTP7  = meanRPM[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 == 0 & meanRPM$Inj_R1_TP6 == 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]


TC_expressedFromTP2 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 > 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP3 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 > 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP4 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 > 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]
TC_expressedFromTP5 = backgroundSubtracted_zygotic[which(meanRPM$Inj_R1_TP1 == 0 & meanRPM$Inj_R1_TP2 == 0 & meanRPM$Inj_R1_TP3 == 0 & meanRPM$Inj_R1_TP4 == 0 & meanRPM$Inj_R1_TP5 > 0 & meanRPM$Inj_R1_TP6 > 0 & meanRPM$Inj_R1_TP7 > 0 & meanRPM$Inj_R1_TP8 > 0 & meanRPM$Inj_R1_TP9 > 0),]

TC_expressedFromTP2_melt = melt(TC_expressedFromTP2[,c(7:15)])
TC_expressedFromTP2_melt$gene = TC_expressedFromTP2$V4
ggpubr::ggline(data = TC_expressedFromTP2_melt,x = 'variable',y = 'value',group = 'gene',col='gene',size=1) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate")  + scale_color_manual(values = palette_custom)
#ggpubr::ggboxplot(data = TC_expressedFromTP2_melt,x = 'variable',y = 'value',size=1) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate")  #+ scale_color_manual(values = palette_custom)


TC_expressedFromTP3_melt = melt(TC_expressedFromTP3[,c(7:15)])
TC_expressedFromTP3_melt$gene = TC_expressedFromTP3$V4
ggpubr::ggline(data = TC_expressedFromTP3_melt,x = 'variable',y = 'value',group = 'gene',col='gene',size=1)  + scale_color_manual(values = palette_custom) + theme_ameres(type = "barplot") + xlab("Timepoint") + ylab("TC conversion rate")

# TC_expressedFromTP4_melt = melt(TC_expressedFromTP4[,c(7:15)])
# TC_expressedFromTP4_melt$gene = TC_expressedFromTP4$V4
# ggplot(TC_expressedFromTP4_melt,aes(x=variable,y=value,group=gene,col=gene)) + geom_line() + scale_color_manual(values = palette_custom) 
dev.off()

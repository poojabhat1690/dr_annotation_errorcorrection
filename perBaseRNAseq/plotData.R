

theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}

library(dplyr)
library(reshape)
## per base analysis of RNAseq around counting windows >5cpm
## the input for this has been generated by creating200ntBins.R
## the counts are obtained from getCountsPerBase.sh

### mean counts 
plusStrand = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW/meanCounts_plus.bed",stringsAsFactors = F)
minusStrand = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW//meanCounts_minus.bed",stringsAsFactors = F)
plusStrand$V1 = plusStrand$V1/24
minusStrand$V1 = minusStrand$V1/24

### each counting window has 200 entries. The order of the counting windows is the same as used for the input for coverage. 

inputReference_plus = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW//countingWindows_greaterThan5CPM_plus_200nt.bed",stringsAsFactors = F,sep="\t")
inputReference_minus = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/RNAseqSignal_endOfCW//countingWindows_greaterThan5CPM_minus_200nt.bed",stringsAsFactors = F,sep="\t")

### read in the duplicated counting windows 


allNames = c(inputReference_minus$V4, inputReference_plus$V4)

#### plus strand :
## get the frequncies of each of the genes 
colnames(plusStrand) = "counts"


freqsGenes = as.data.frame(table(paste0(inputReference_plus$V4)))

library(splitstackshape)
library(mefa)
plusStrand$position = c(1:200)
plusStrand$counts = plusStrand$counts + 0.1

metaData_plus = inputReference_plus %>% slice(rep(1:nrow(inputReference_plus), each=200))

plusStrand = cbind.data.frame(plusStrand,metaData_plus)

plusStrand = merge(plusStrand,freqsGenes,by.x="V4",by.y="Var1")


###### minus strand


colnames(minusStrand) = "counts"

freqsGenes = as.data.frame(table(paste0(inputReference_minus$V4)))
minusStrand$position = seq(200,1,-1)
minusStrand$counts = minusStrand$counts + 0.1
metaData_minus = inputReference_minus %>% slice(rep(1:nrow(inputReference_minus), each=200))

minusStrand = cbind.data.frame(minusStrand,metaData_minus)


minusStrand = merge(minusStrand,freqsGenes,by.x="V4",by.y="Var1")


plusMinus =  rbind(plusStrand,minusStrand)




######################

#### 1 counting window

#######################

plusMinus_1 = plusMinus[which(plusMinus$Freq == 1),]

plusMinus_1_split = split(plusMinus_1,plusMinus_1$V4,T)
plusMinus_1_split = lapply(plusMinus_1_split,function(x) x[order(x$position),])
plusMinus_1_splitCounts = lapply(plusMinus_1_split,function(x) x$counts/sum(x$counts))
plusMinus_1_split_melt= melt(plusMinus_1_splitCounts)
plusMinus_1_split_melt$position = c(1:200)
plusMinus_1_split_melt_splitPosition= split(plusMinus_1_split_melt,plusMinus_1_split_melt$position,T)

perPosition_median_1 = melt(lapply(plusMinus_1_split_melt_splitPosition,function(x) median(x$value)))
perPosition_median_quantile25  = melt(lapply(plusMinus_1_split_melt_splitPosition,function(x) quantile(x$value)[2]))
perPosition_median_quantile75  = melt(lapply(plusMinus_1_split_melt_splitPosition,function(x) quantile(x$value)[4]))

perPosition_median_1_25_27 = list(perPosition_median_1,perPosition_median_quantile25,perPosition_median_quantile75)
names(perPosition_median_1_25_27)  = c("median","quantile25","quantile75")
perPosition_median_1_25_27 = melt(perPosition_median_1_25_27)
perPosition_median_1_25_27$position = c(1:200)
perPosition_median_1_25_27$value = perPosition_median_1_25_27$value * 100

perPosition_median_1_25_27_1_us = perPosition_median_1_25_27 %>% filter(position < 101)
perPosition_median_1_25_27_1_ds = perPosition_median_1_25_27 %>% filter(position >= 101)



p = ggplot(perPosition_median_1_25_27_1_us,aes(x=position,y=value,group=L1,col=L1)) + geom_line(size=1) + scale_color_brewer(palette = "Dark2")+ theme_ameres(type = "barplot")
ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_1_signal_us.pdf",plot = p,height=4,width=5)

q = ggplot(perPosition_median_1_25_27_1_ds,aes(x=position,y=value,group=L1,col=L1)) + geom_line(size=1) + scale_color_brewer(palette = "Dark2")+ theme_ameres(type = "barplot")+ylim(c(0,2))
ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_1_signal_ds.pdf",plot = q,height=4,width=5)


###### i want to plot the median value of signal at each point... 


 p = ggpubr::ggline(data = perPosition_median_1,x = 'L1',y = 'value',size = 1,col='red') + theme_ameres(type = "barplot") + scale_x_discrete(breaks = c(0,20,40,60,80,100,120,140,160,180,200,220)) + ggtitle(label = nrow(plusMinus_1_split_melt_splitPosition$`1`))
 dir.create(path = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/")
 ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_1_signal.pdf",plot = p,height=4,width=5)


fractionPerPosition  = plusMinus_1_splitCounts
fractionPerPosition_df = do.call(cbind.data.frame,fractionPerPosition)
#write.table(fractionPerPosition_df,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/finalEvaluationOfEnds/perBaseRNAseq/dataFrames/1CountingWindowPergene.txt",sep="\t",quote=F)
#write.table(fractionPerPosition_df,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/figs7/fig_f/1CountingWindowPergene.txt",sep="\t",quote=F)





############## 2 priming sites

plusMinus_2 = plusMinus[which(plusMinus$Freq == 2),]
## here
plusMinus_2$id = paste0(plusMinus_2$V4,plusMinus_2$V1,plusMinus_2$V2,plusMinus_2$V3)

plusMinus_2_split = split(plusMinus_2,plusMinus_2$id,T)
plusMinus_2_split = lapply(plusMinus_2_split,function(x) x[order(x$position),])
plusMinus_2_splitCounts = lapply(plusMinus_2_split,function(x) x$counts/sum(x$counts))
plusMinus_2_split_melt= melt(plusMinus_2_splitCounts)
plusMinus_2_split_melt$position = c(1:200)
plusMinus_2_split_melt_splitPosition= split(plusMinus_2_split_melt,plusMinus_2_split_melt$position,T)

perPosition_median_2 = melt(lapply(plusMinus_2_split_melt_splitPosition,function(x) median(x$value)))


perPosition_median_2_quantile25  = melt(lapply(plusMinus_2_split_melt_splitPosition,function(x) quantile(x$value)[2]))
perPosition_median_2_quantile75  = melt(lapply(plusMinus_2_split_melt_splitPosition,function(x) quantile(x$value)[4]))

perPosition_median_2_25_75 = list(perPosition_median_2,perPosition_median_2_quantile25,perPosition_median_2_quantile75)
names(perPosition_median_2_25_75)  = c("median","quantile25","quantile75")
perPosition_median_2_25_75 = melt(perPosition_median_2_25_75)
perPosition_median_2_25_75$position = c(1:200)
perPosition_median_2_25_75$value = perPosition_median_2_25_75$value * 100

ggpubr::ggline(data = perPosition_median_2_25_75,x = 'position',y='value',group = 'L1',col='L1') + scale_color_brewer(palette = 'Dark2') + theme_ameres(type = 'barplot') 

ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_12_allTogether_signal.pdf",plot = p,height=4,width=5)



### split plus and minus 

plus_2 = plusMinus_2 %>% filter(V6 == "+")
minus_2 = plusMinus_2 %>% filter(V6 == "-")

plus_2 = plus_2[order(plus_2$position),]
minus_2 = minus_2[order(minus_2$position),]

plus_2_split = split(plus_2,plus_2$V4,T)

## order proximal and distal

plus_2_split = lapply(plus_2_split,function(x) x[order(x$V3,decreasing = F),])

plus_2_ordered_df = do.call(rbind.data.frame,plus_2_split)

plus_2_ordered_df = as.data.frame(plus_2_ordered_df  %>% dplyr::group_by(V4) %>% mutate(fraction=counts/sum(counts)))


plus_2_ordered_df$category = rep(c("proximal","distal"),each=200)
plus_2_ordered_df$positionGene = c(1:400)
#ggplot(plus_2_ordered_df,aes(x=positionGene,y=fraction,group=positionGene)) +geom_boxplot(outlier.shape=NA,coef = 0) + ylim(c(0,0.02))


### ssame fo minus 

minus_2_split = split(minus_2,minus_2$V4,T)

## order proximal and distal


minus_2_split = lapply(minus_2_split,function(x) x[order(x$V2,decreasing = T),])


###

minus_2_ordered_df = do.call(rbind.data.frame,minus_2_split)
minus_2_ordered_df = as.data.frame(minus_2_ordered_df  %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)))
minus_2_ordered_df$category = rep(c("proximal","distal"),each=200)

minus_2_ordered_df$positionGene = c(seq(1,400))
#ggplot(minus_2_ordered_df,aes(x=positionGene,y=fraction,group=positionGene)) +geom_boxplot(outlier.shape=NA,coef = 0) + ylim(c(0,0.02))



allOrdered = rbind(minus_2_ordered_df,plus_2_ordered_df)
allOrdered_df = as.data.frame(allOrdered  %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)))

all2_boxplot  = ggplot(allOrdered_df,aes(x=positionGene,y=fraction,group=positionGene)) +geom_boxplot(outlier.shape=NA,coef = 0) + ylim(c(0,0.0075))+ ggtitle(paste0("number of genes (2 counting windows)=",nrow(plusMinus_2)/400))+ xlab("position")  + ylab("Relative RNAseq signal") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

allOrdered_df_split = split(allOrdered_df,allOrdered_df$positionGene,T)
allOrdered_df_split_melt = melt(lapply(allOrdered_df_split,function(x) median(x$fraction)))


allOrdered_df_split_melt_quantile25  = melt(lapply(allOrdered_df_split,function(x) quantile(x$fraction)[2]))
allOrdered_df_split_melt_quantile75  = melt(lapply(allOrdered_df_split,function(x) quantile(x$fraction)[4]))

allOrdered_df_split_melt_2_25_75 = list(allOrdered_df_split_melt,allOrdered_df_split_melt_quantile25,allOrdered_df_split_melt_quantile75)
names(allOrdered_df_split_melt_2_25_75)  = c("median","quantile25","quantile75")
allOrdered_df_split_melt_2_25_75 = melt(allOrdered_df_split_melt_2_25_75)
#perPosition_median_2_25_75$position = c(1:200)
allOrdered_df_split_melt_2_25_75$value = allOrdered_df_split_melt_2_25_75$value * 100

allOrdered_df_split_melt_2_25_75$position = c(1:400)
plot(allOrdered_df_split_melt$L1,allOrdered_df_split_melt$value)
p = ggplot(allOrdered_df_split_melt_2_25_75,aes(x=position,y=value,group=L1,col=L1)) + geom_line(size=1) + scale_color_brewer(palette = "Dark2")+ theme_ameres(type = "barplot") 
p = p + theme(axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
#p = ggpubr::ggline(data = allOrdered_df_split_melt_2_25_75,x = 'position',y = 'value',group = 'L1',col= 'L1') + theme_ameres(type = 'barplot')
ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_2_spearate_signal.pdf",plot = p,height=4,width=9)


#### 3 counting windows


plusMinus_3 = plusMinus[which(plusMinus$Freq == 3),]
plusMinus_3$id = paste0(plusMinus_3$V4,plusMinus_3$V1,plusMinus_3$V2,plusMinus_3$V3)

plusMinus_3_split = split(plusMinus_3,plusMinus_3$id,T)
plusMinus_3_split = lapply(plusMinus_3_split,function(x) x[order(x$position),])



## here
plusMinus_3 = as.data.frame(plusMinus_3 %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)) )
plusMinus_3_split = split(plusMinus_3,plusMinus_3$position,T)

plusMinus_3_split_melt= melt(lapply(plusMinus_3_split,function(x) median(x$fraction)))


##### write table 

fractionPerPosition  = lapply(plusMinus_3_split,function(x) x$fraction)
fractionPerPosition_df = do.call(cbind.data.frame,fractionPerPosition)
#write.table(fractionPerPosition_df,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/finalEvaluationOfEnds/perBaseRNAseq/dataFrames/3CountingWindowPergene.txt",sep="\t",quote=F)
#write.table(fractionPerPosition_df,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/figs7/fig_l//3CountingWindowPergene.txt",sep="\t",quote=F)





all3 = ggplot(plusMinus_3_split_melt,aes(x=as.numeric(L1),y=value))+ geom_line()



p = ggpubr::ggline(data = plusMinus_3_split_melt,x = 'L1',y = 'value',size = 1,col='red') + theme_ameres(type = "barplot") + scale_x_discrete(breaks = c(0,20,40,60,80,100,120,140,160,180,200,220)) + ggtitle(label = nrow(plusMinus_1_split_melt_splitPosition$`1`))
ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_3_together_signal.pdf",plot = p,height=4,width=5)


### split plus and minus 

plus_3 = plusMinus_3 %>% filter(V6 == "+")
minus_3 = plusMinus_3 %>% filter(V6 == "-")


plus_3 = plus_3[order(plus_3$position),]
minus_2 = plus_3[order(plus_3$position),]

plus_3_split = split(plus_3,plus_3$V4,T)

## order proximal and distal

plus_3_split = lapply(plus_3_split,function(x) x[order(x$V3,decreasing = F),])

plus_3_ordered_df = do.call(rbind.data.frame,plus_3_split)

plus_3_ordered_df = as.data.frame(plus_3_ordered_df  %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)))


plus_3_ordered_df$category = rep(c("proximal","distal","distal1"),each=200)
plus_3_ordered_df$positionGene = c(1:600)
#ggplot(plus_3_ordered_df,aes(x=positionGene,y=fraction,group=positionGene)) +geom_boxplot(outlier.shape=NA,coef = 0) + ylim(c(0,0.02))


### ssame fo minus 

minus_3_split = split(minus_3,minus_3$V4,T)

## order proximal and distal


minus_3_split = lapply(minus_3_split,function(x) x[order(x$V2,decreasing = T),])


minus_3_ordered_df = do.call(rbind.data.frame,minus_3_split)
minus_3_ordered_df = as.data.frame(minus_3_ordered_df  %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)))
minus_3_ordered_df$category = rep(c("proximal","distal","distal1"),each=200)
#minus_3_ordered_df$positionGene = c(seq(200,1,-1),seq(400,201,-1),seq(600,401,-1))
minus_3_ordered_df$positionGene = c(1:600)

###ggplot(minus_3_ordered_df,aes(x=position,y=fraction,group=position)) +geom_boxplot(outlier.shape=NA,coef = 0) + ylim(c(0,0.02))



allOrdered = rbind(plus_3_ordered_df,minus_3_ordered_df)
allOrdered_df = as.data.frame(allOrdered  %>% group_by(V4) %>% mutate(fraction=counts/sum(counts)))
#allOrdered_df$category = rep(c("proximal","distal"),each=200)
all3_boxplot = ggplot(allOrdered_df,aes(x=positionGene,y=fraction,group=positionGene)) +geom_boxplot(outlier.shape=NA,coef = 0) + ggtitle(paste0("number of genes (3 counting windows)=",nrow(plusMinus_3)/600))+ xlab("position") + ylim(c(0,0.005)) + ylab("Relative RNAseq signal") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
allOrdered_df_split = split(allOrdered_df,allOrdered_df$positionGene,T)
allOrdered_df_split_median = melt(lapply(allOrdered_df_split,function(x) median(x$fraction)))
p = ggpubr::ggline(data = allOrdered_df_split_median,x = 'L1',y = 'value',size = 1,col='red') + theme_ameres(type = "barplot") + scale_x_discrete(breaks = c(0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600)) + ggtitle(label = length(unique(allOrdered_df$V4)))
ggsave(filename = "/Users/pooja.bhat/Dropbox/Paperoutline/analysis/UTRannotation/plots/perBaseRNAseq/cw_3_sepatate_signal.pdf",plot = p,height=4,width=7)


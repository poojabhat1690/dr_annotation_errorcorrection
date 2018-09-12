##### reading in the conversion rates per UTR for all the time points. 

library(plyr)

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



###### making TPMs and count data from SLAMseq samples ... 

countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/utrRates//",pattern = ".csv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis//annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/utrRates///",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T,stringsAsFactors = F))

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



library(dplyr)
library(tidyr)
getTCrates = function(dataFrame){
  
  
  
  T_C_all = vector("list",length(dataFrame))
  names(T_C_all) = names(dataFrame)
  
  allRates = vector("list",length(dataFrame))
  names(allRates) = names(dataFrame)
  
  
  
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
    allRates[[i]] =  plotTab
    
  }
  
  return(allRates)
}


TCrates_overlall = getTCrates(dataFrame = countData_timeCourse_datasets)


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/conversions_allData.pdf")
for(i in 1:length(TCrates_overlall)){
  
  p = ggplot2::ggplot(TCrates_overlall[[i]],aes(x=class,y=log10(values + 0.0001),group= class)) + geom_boxplot()+ theme_bw() + ggtitle(names(TCrates_overlall)[i])
  p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) 
  p = p + xlab("Timepoint") + ylab("Mutation rate per UTR base [%]")
  print(p)
  
}

dev.off()




########################### now i want to look at classified zygotic transcripts and their conversion rates... 

allTranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/countingWindows_classified.txt",sep="\t",stringsAsFactors = F)
allTranscripts_zygotic = allTranscripts[which(allTranscripts$V7  == "Z"),]

colnames(allTranscripts_zygotic) = c("Chr","Start","End","Name","V5","Strand","V7")
conversions_zygoticTranscripts = lapply(countData_timeCourse_datasets,function(x) left_join(allTranscripts_zygotic , x))
conversions_zygoticTranscripts = lapply(conversions_zygoticTranscripts,function(x) x %>% select(-c("V5","V7")))
TCrates_zygotic = getTCrates(dataFrame = conversions_zygoticTranscripts)



pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/conversions_zygoticTranscripts.pdf")
for(i in 1:length(TCrates_overlall)){
  
  p = ggplot2::ggplot(TCrates_zygotic[[i]],aes(x=class,y=values ,group= class)) + geom_boxplot(outlier.shape=NA)+ theme_bw() + ggtitle(names(TCrates_overlall)[i])
  p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) 
  p = p + xlab("Timepoint") + ylab("Mutation rate per UTR base [%]")
  print(p)
  
}

dev.off()




pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/conversions_zygoticTranscripts_densities.pdf")
for(i in 1:length(TCrates_zygotic)){
  
  p = ggplot2::ggplot(TCrates_zygotic[[i]],aes(log10(values) ,group= class,col=class)) + geom_density(size=1)+ theme_bw() + ggtitle(names(TCrates_overlall)[i])
  p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) 
  p = p + xlab("Timepoint") + ylab("Mutation rate per UTR base [%]")
  print(p)
  
}

dev.off()


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/conversions_zygoticTranscripts_densities_onlyTC.pdf")

for(i in 1:length(TCrates_zygotic)){
  currentFile = TCrates_zygotic[[i]]
  currentFile = currentFile[which(currentFile$class== "T_C"),]
  p = ggplot2::ggplot(currentFile,aes(log10(values) ,group= class,col=class)) + geom_density(size=1)+ theme_bw() + ggtitle(names(TCrates_overlall)[i])
  p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) 
  p = p + xlab("Timepoint") + ylab("Mutation rate per UTR base [%]")
  print(p)
}

dev.off()


########### i want to compare the background levels of TC conversions to conversions of other nucleotides in the treated sample.. 


timepoints = paste0("TP",c(1:9))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/untreatedTC_otherMutations_distribution.pdf")
for(i in 1:length(timepoints)){

TP1_untreated = TCrates_zygotic[grep(paste0("Untreated_",timepoints[i]),names(TCrates_zygotic))]
TP1_untreated_TC = TP1_untreated[[1]] %>% filter(class == "T_C")
TP1_untreated_notTC = TP1_untreated[[1]] %>% filter(class != "T_C")

TP1_injected = TCrates_zygotic[which(names(TCrates_zygotic) == paste0("Inj_R1_",timepoints[i])  | names(TCrates_zygotic) == paste0("Inj_R2_",timepoints[i]) | names(TCrates_zygotic) == paste0("Inj_R3_",timepoints[i]) )]

a = do.call(rbind,TP1_injected)
a_copy = a
a = a %>% filter(class != "T_C")
a_TC = a_copy %>% filter(class == "T_C")
a$category = "treated_other"
TP1_untreated_TC$category = "TC_untreated"
TP1_untreated_notTC$category= "other_untreated"
a_TC$category = "TC_treated"
total_samples = rbind(a,TP1_untreated_TC,a_TC,TP1_untreated_notTC)
total_samples$log10Vals = log10(total_samples$values)

p = ggpubr::ggdensity(data = total_samples, 'log10Vals',col='category',size=1) + scale_color_brewer(palette = "Dark2") + theme_ameres(type = "barplot") + ggtitle(i)
print(p)
my_comparisons <- list( c("treated_other", "TC_untreated"), c("treated_other", "other_untreated"), c("treated_other", "TC_treated"),c("TC_untreated","other_untreated"),c("TC_untreated","TC_treated"),c( "other_untreated","TC_treated") )
library(ggpubr)
q = ggplot2::ggplot(total_samples,aes(x=category,y=values,group=category)) + geom_boxplot(outlier.shape = NA)   + ggtitle(i)
q = q + stat_compare_means(comparisons = my_comparisons)
print(q)


}

dev.off()

####### i Also want to correlate the conversion rates in samples... in the GRANDSLAM paper, the other conversions are taken as proxy 
########### because three is a correlation between TC rates in the background and 


convertPlusMinus = function(dataFrame){
  
   T_C_all = vector("list",length(dataFrame))
  names(T_C_all) = names(dataFrame)
  
  allRates = vector("list",length(dataFrame))
  names(allRates) = names(dataFrame)
  
  
  
  for(i in 1:length(dataFrame))
  {
    
    curTab = dataFrame[[i]]
    plusTab = curTab %>% dplyr::filter(Strand == "+")
    minusTab = curTab %>% dplyr::filter(Strand == "-") %>%
      dplyr::select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)
    
    plusTab = plusTab %>% dplyr::select(-contains("N")) %>%
      dplyr::select(-one_of("Chr","Start")) %>%
      
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
    
    
    minusTab = minusTab %>% 
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
    plotTab_minusTC = plotTab[,-which(colnames(plotTab) == "T_C")]
    plotTab_minusTC = rowMeans(plotTab_minusTC)
    plotTab_new = cbind.data.frame(plotTab_minusTC,plotTab$T_C)
    colnames(plotTab_new) = c("otherRates","TC")
    allRates[[i]] = plotTab
    cat(i)
  }
  return(allRates)
}

countData_minusConverted = convertPlusMinus(dataFrame = countData_timeCourse_datasets)
countData_minusConverted_untreated = countData_minusConverted[grep("Unt",names(countData_minusConverted))]
corr_TC_others_untreatedSamples = vector("list",length(countData_minusConverted_untreated))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/correlation_untreatedsamples.pdf")
for(i in 1:length(countData_minusConverted_untreated)){
  sample_use = paste0("Untreated_",timepoints[i])
  corr_TC_others = cor(countData_minusConverted_untreated[[i]],method = "spearman",use = "complete.obs")
  corr_TC_others_untreatedSamples[[i]] = corr_TC_others
  corrplot::corrplot(corr_TC_others)
}
dev.off()

names(corr_TC_others_untreatedSamples) = names(countData_minusConverted_untreated)

#####now I also want to correlate the TC conversions in the background samples with other conversiosn in treated samples.. 
countData_minusConverted_treated = countData_minusConverted[grep("Inj",names(countData_minusConverted))]
cor(countData_minusConverted_treated$Inj_R1_TP9$C_T,countData_minusConverted_untreated$Untreated_TP9$T_C,use ="complete.obs",method = "spearman")

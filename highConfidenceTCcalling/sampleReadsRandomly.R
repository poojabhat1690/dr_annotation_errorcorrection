library(matrixStats)
#### random sampling of reads... 
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




getProbability = function(allReads_zygoticTranscripts,condition){
  
  
  ###### changing to the reverse complement on the minus strand
  minusTab = allReads_zygoticTranscripts %>% dplyr::filter(strand == "-") %>%
    dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)
  plusTab = allReads_zygoticTranscripts %>% dplyr::filter(strand == "+") %>%
    dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G )
  
  allReads = bind_rows(minusTab,plusTab)
  
  
  library(tidyverse)
  
  
  sampleSize = c(1,2,5,10,100,1000,5000,10000,100000)
  subsamplingVectors_2Ormore = vector("list",length(sampleSize))
  names(subsamplingVectors_2Ormore) = sampleSize
  subsamplingVectors_2Ormore =lapply(subsamplingVectors_2Ormore,function(x) vector("list",100))
  subsamplingVectors_1 = subsamplingVectors_2Ormore
  subsamplingVectors_0 = subsamplingVectors_2Ormore
  
  for(i in 1:length(sampleSize)){
    for(j in 1:100){
      nucleotideOccurance = sample_n(allReads, size = sampleSize[i], replace = FALSE, weight = NULL ) %>% select(1:12) 
      subsamplingVectors_2Ormore[[i]][[j]] = apply(nucleotideOccurance,2,function(x) length(which(x>=2))/sampleSize[i])
      subsamplingVectors_1[[i]][[j]] = apply(nucleotideOccurance,2,function(x) length(which(x==1))/sampleSize[i])
      subsamplingVectors_0[[i]][[j]] = apply(nucleotideOccurance,2,function(x) length(which(x==0))/sampleSize[i])
    }
    
  }
  
  cols_use = c(RColorBrewer::brewer.pal(name = "Dark2",8),RColorBrewer::brewer.pal(name = "Set1",4) )
  
  subsamplingVectors_2Ormore_perReadDepth = do.call(rbind,lapply(subsamplingVectors_2Ormore,function(x) colMeans(do.call(rbind,x) )))
  subsamplingVectors_1_perReadDepth = do.call(rbind,lapply(subsamplingVectors_1,function(x) colMeans(do.call(rbind,x) )))
  subsamplingVectors_1_perReadDepth_melt= melt(subsamplingVectors_1_perReadDepth)
  subsamplingVectors_2Ormore_perReadDepth_melt= melt(subsamplingVectors_2Ormore_perReadDepth)
  
  subsamplingVectors_0_perReadDepth = do.call(rbind,lapply(subsamplingVectors_0,function(x) colMeans(do.call(rbind,x) )))
  subsamplingVectors_0_perReadDepth_melt= melt(subsamplingVectors_0_perReadDepth)
  
  
  pdf(paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/probabilityOfTCreads_",condition,".pdf"),height = 5,width=7)
  p = ggpubr::ggline(data = subsamplingVectors_1_perReadDepth_melt,x = 'X1',y = 'value',group = 'X2',col='X2',size = 1)+ scale_color_manual(values = cols_use) + ggtitle("probability of 1 mutation") + theme_ameres(type = "barplot") + ylab("Probability of finding read") + xlab("Read depth")
  print(p)
  p = ggpubr::ggline(data = subsamplingVectors_2Ormore_perReadDepth_melt,x = 'X1',y = 'value',group = 'X2',col='X2',size = 1) + scale_color_manual(values = cols_use) + ggtitle("probability of >=2 mutation") + theme_ameres(type = "barplot")+ ylab("Probability of finding read")+ xlab("Read depth")
  print(p)
  p = ggpubr::ggline(data = subsamplingVectors_0_perReadDepth_melt,x = 'X1',y = 'value',group = 'X2',col='X2',size = 1) + scale_color_manual(values = cols_use) + ggtitle("probability of 0 mutation")+ theme_ameres(type = "barplot")+ ylab("Probability of finding read")+ xlab("Read depth")
  print(p)
  dev.off()
  
  
  
}


### zygotic + mitochondrial counting windows
allReads_zygoticTranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/randomProbabilities/onlyZygoticTranscriptsallcombinedBam.txt",stringsAsFactors = F)
allReads_zygoticTranscripts_chrM =allReads_zygoticTranscripts[grep("chrM",allReads_zygoticTranscripts$id),]
allReads_zygoticTranscripts = allReads_zygoticTranscripts[-grep("chrM",allReads_zygoticTranscripts$id),]
getProbability(allReads_zygoticTranscripts = allReads_zygoticTranscripts,condition = "onlyZygotic")
getProbability(allReads_zygoticTranscripts = allReads_zygoticTranscripts_chrM,condition = "mt")

#### all counting windows (no zygotic)

allReads_all = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/randomProbabilities/allcombinedBam.txt",stringsAsFactors = F)
allReads_all =  allReads_all[-grep("chrM",allReads_all$id),]
getProbability(allReads_zygoticTranscripts = allReads_all,condition = "allCw")


#### https://stackoverflow.com/questions/48612153/how-to-calculate-confidence-intervals-for-a-vector
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

confidence_interval(a,0.95)

### or just use the t distribution to find the condidence interval 
t.test(a,conf.level = 0.95)


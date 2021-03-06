##### to determine the number of TC reads to use.. I want calculate the fraction of reads with TC conversions in the background 
##### and compare this to the fraction of TC conversions in the treatment samples...

library(dplyr)


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


######### preparing the data tables...of the fraction of reads with 1,2,...n TCs... (run on cluster)

####### getting the fraction of reads separately for mitochondrial genes
# zygoticGenes = read.table("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)
# countingWindws = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined/allAnnotations_withPriMirna.bed",sep="\t",stringsAsFactors = F)
# countingWindws_zygotic = countingWindws[countingWindws$V4 %in% zygoticGenes$external_gene_name,]
# countingWindws_zygotic$V2 = countingWindws_zygotic$V2 + 1
# countingWindws_zygotic$id  = paste0(countingWindws_zygotic$V1,":",countingWindws_zygotic$V2)
# 
# 
# folder = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/"
# treatmentSamples = list.files(path = folder,pattern = "freqAllMutations_allTranscripts_totalSNP*")
# treatmentSamples_Inj = treatmentSamples[grep("Inj",treatmentSamples)] %>% paste0(folder,.)
# treatmentSamples_Unt = treatmentSamples[grep("Unt",treatmentSamples)] %>% paste0(folder,.)
# 
# getTabulatedTcs = function(sampleName){ ###3 get fraction of reads with 1TC, 2 TCs... etc...
#   sampleData = read.table(sampleName)
#   sampleData$id = unlist(lapply(strsplit(as.character(sampleData$id),"-",T),function(x) x[[1]]))
#   sampleData_mt = filter(sampleData, grepl("chrM",id))
#   sampleData = sampleData[sampleData$id %in% countingWindws_zygotic$id,]
#   minusTab = sampleData %>% dplyr::filter(strand == "-") %>%
#     dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,id,strand)
#   plusTab = sampleData %>% dplyr::filter(strand == "+") %>%
#     dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,id,strand)
#   sampleData  = rbind(plusTab,minusTab)
#   #sampleData = filter(sampleData, !grepl("chrM",id))
# 
#   ### all excluding mitochrondria
#   #fractionTC_samples = table(sampleData$T_C)/nrow(sampleData)
#   sampleData = sampleData %>% select(-matches('A_A|G_G|C_C|T_T'))
#   sampleData_mt = sampleData_mt %>% select(-matches('A_A|G_G|C_C|T_T'))
# 
#   fractionTC_samples = apply(sampleData[,1:12],2,function(x) table(x)/length(x))
#   totalReads =   apply(sampleData[,1:12],2,function(x) table(x))
#   fractionTC_samples_reads = vector("list",2)
#   names(fractionTC_samples_reads) = c("fractionTC","numberOfreads")
#   fractionTC_samples_reads[[1]] = fractionTC_samples
#   fractionTC_samples_reads[[2]] = totalReads
# 
#   ##### only mitochondria
#   #fractionTC_samplesMt = table(sampleData_mt$T_C)/nrow(sampleData_mt)
#   fractionTC_samplesMt = apply(sampleData_mt[,1:12],2,function(x) table(x)/length(x))
#   totalReads_mt =   apply(sampleData_mt[,1:12],2,function(x) table(x))
# 
#   #totalReads_mt = table(sampleData_mt$T_C)
#   fractionTC_samples_reads_mt = vector("list",2)
#   names(fractionTC_samples_reads_mt) = c("fractionTC","numberOfreads")
#   fractionTC_samples_reads_mt[[1]] = fractionTC_samplesMt
#   fractionTC_samples_reads_mt[[2]] = totalReads_mt
# 
#   totalFracs = list(fractionTC_samples_reads_mt,fractionTC_samples_reads)
#   names(totalFracs) = c("mt","allOther")
#   return(totalFracs)
# }
# 
# 
# sample_TCs = lapply(treatmentSamples_Unt,function(x) getTabulatedTcs(sampleName = x))
# names(sample_TCs) = treatmentSamples_Unt
# save(sample_TCs ,file = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples_zygoticGenes.Rdata")
# 
# sample_TCs_injection = lapply(treatmentSamples_Inj,function(x) getTabulatedTcs(sampleName = x))
# names(sample_TCs_injection) = treatmentSamples_Inj
# save(sample_TCs_injection ,file = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples_zygoticGenes.Rdata")


####### plotting these datasets...

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples_zygoticGenes.Rdata")
names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 211,stop = 224)

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples_zygoticGenes.Rdata")
names(sample_TCs) = substr(names(sample_TCs),start = 211,stop = 224)


#### for each timepoint... 
library(data.table)
library(ggplot2)
library(scales)


#########################
TPs = paste0("TP",c(1:9))
library(RColorBrewer)
colsUse = c(brewer.pal(n = 8,"Dark2"),brewer.pal(n = 4,"Set1"))

pdf("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/otherMutations_compareTC_zygoticGenes.pdf",height = 5,width = 5)
for(i in 1:length(TPs)){
  
  load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples_zygoticGenes.Rdata")
  names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 227,stop = 240)
  
  load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples_zygoticGenes.Rdata")
  names(sample_TCs) = substr(names(sample_TCs),start = 230,stop = 243)
  
  sample_TCs_injection_total = sample_TCs_injection
  sample_TCs_injection = lapply(sample_TCs_injection,function(x) x$allOther)
  
  sample_TCs_TP = sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))] 
  
  sample_TCs_TP = lapply(sample_TCs_TP,function(x) x$numberOfreads)
  
  p =  lapply(sample_TCs_TP,function(x) lapply(x,function(y) data.frame(y))) %>% lapply(. %>% plyr::ldply(., rbind)) %>% 
    lapply(. %>% mutate(id = paste0(.id,"_",x))) %>% purrr::reduce(., full_join, by = 'id') %>% replace(., is.na(.), 0) %>% 
    tidyr::separate(data = .,col = id,sep="_",into=c('firstNt','lastNt','numberNt')) %>% 
    mutate(totalReads = Freq.x + Freq.y + Freq)  %>% dplyr::select(-matches('x.x|x.y|.id.x|.id.y')) %>%
    mutate(id_final = paste(firstNt,lastNt,sep="_")) %>% dplyr::group_by(id_final) %>% mutate(fractionReads = totalReads/sum(totalReads)) %>%
    mutate(numberNt = as.numeric(numberNt)) %>% ggpubr::ggline(.,x = 'numberNt',y='fractionReads',group='id_final',col='id_final' ) +  
    scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = colsUse) + ggtitle(paste0("allGenes-",TPs[i])) + theme_ameres(type = "barplot") + xlab("Number of conversions") + ylab("Fraction of reads")
  print(p)
  
  sample_TCs_TP_mt = lapply(sample_TCs_injection_total,function(x) x$mt)
  
  
  sample_TCs_TP_mt = sample_TCs_TP_mt[grep(TPs[i],names(sample_TCs_TP_mt))] 
  
  sample_TCs_TP_mt = lapply(sample_TCs_TP_mt,function(x) x$numberOfreads)
  
  
  r = lapply(sample_TCs_TP_mt,function(x) lapply(x,function(y) data.frame(y))) %>% lapply(. %>% plyr::ldply(., rbind)) %>% 
    lapply(. %>% mutate(id = paste0(.id,"_",x))) %>% purrr::reduce(., full_join, by = 'id') %>% replace(., is.na(.), 0) %>% 
    tidyr::separate(data = .,col = id,sep="_",into=c('firstNt','lastNt','numberNt')) %>% 
    mutate(totalReads = Freq.x + Freq.y + Freq) %>% dplyr::select(-matches('x.x|x.y|.id.x|.id.y')) %>%
    mutate(id_final = paste(firstNt,lastNt,sep="_")) %>% dplyr::group_by(id_final) %>% mutate(fractionReads = totalReads/sum(totalReads)) %>%
    mutate(numberNt = as.numeric(numberNt)) %>% ggpubr::ggline(.,x = 'numberNt',y='fractionReads',group='id_final',col='id_final' ) +  
    scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))+
    scale_color_manual(values = colsUse)+ ggtitle(paste0("mt-",TPs[i])) + theme_ameres(type = "barplot")+ xlab("Number of conversions") + ylab("Fraction of reads")
  
  print(r)
  
  #### doing the same for the background data...
  
  sample_TCs_TP_bg = lapply(sample_TCs,function(x) x$allOther)
  sample_TCs_TP_bg = sample_TCs_TP_bg[grep(TPs[i],names(sample_TCs_TP_bg))] 
  sample_TCs_TP_bg = sample_TCs_TP_bg[[1]]$numberOfreads
  
  bg_all = lapply(sample_TCs_TP_bg,function(x) data.frame(x)) %>% plyr::ldply(., rbind)%>%
    dplyr::group_by(.data = .,.id) %>% mutate(fractionReads = Freq/sum(Freq)) %>%
    ggpubr::ggline( .,x = 'x',y='fractionReads',group = '.id',col='.id') + 
    scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))+
    scale_color_manual(values = colsUse)+ ggtitle(paste0("BG_allGenes-",TPs[i])) + theme_ameres(type = "barplot")  + xlab("Number of conversions") + ylab("Fraction of reads")
  print(bg_all)
  
  
  sample_TCs_TP_bg_mt = lapply(sample_TCs,function(x) x$mt)
  sample_TCs_TP_bg_mt = sample_TCs_TP_bg_mt[grep(TPs[i],names(sample_TCs_TP_bg_mt))] 
  sample_TCs_TP_bg_mt = sample_TCs_TP_bg_mt[[1]]$numberOfreads
  
  bg_mt = lapply(sample_TCs_TP_bg_mt,function(x) data.frame(x)) %>% plyr::ldply(., rbind)%>%
    dplyr::group_by(.data = .,.id) %>% mutate(fractionReads = Freq/sum(Freq)) %>%
    ggpubr::ggline( .,x = 'x',y='fractionReads',group = '.id',col='.id') + 
    scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))+
    scale_color_manual(values = colsUse)+ ggtitle(paste0("BG_mt-",TPs[i])) + theme_ameres(type = "barplot")  + xlab("Number of conversions") + ylab("Fraction of reads")
  print(bg_mt)
  
  
}

dev.off()





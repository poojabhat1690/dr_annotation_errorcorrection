#### depreciated....please see the rrscript numberOfTCConversionsPerStage.R


numberOfreadFiles = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/",pattern = "*_TC.txt")
numberOfreadFiles_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/",numberOfreadFiles)
numberOfreadFiles_data  = lapply(numberOfreadFiles_path,function(x) read.table(x,stringsAsFactors = F))
names(numberOfreadFiles_data) = numberOfreadFiles
totalNumberOfreads = do.call(rbind.data.frame,numberOfreadFiles_data)


minusTab = totalNumberOfreads %>% dplyr::filter(strand == "-") %>%
  dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,strand,id)
plusTab = totalNumberOfreads %>% dplyr::filter(strand == "+") %>%
  dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,strand,id)

totalReads = rbind.data.frame(plusTab,minusTab)
totalReads_tabulated  = as.data.frame(table(totalReads$T_C))

#### total number of reads in this bam file - 33612806 (from samtools flagstat mappingToUTRs.bam)


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


### reads with no TC conversions have not been considered in this analysis ...

options(scipen=999)
pdf("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/numberOfTcs_errorRate.pdf",height=4.5)
p  =totalReads_tabulated %>% mutate(fractionReads = (Freq * 100)/33612806) %>% filter(Var1 != 0)  %>% mutate(frac_round = round(fractionReads,5)) %>% ggpubr::ggline(x='Var1',y='fractionReads',size=1.5,color = 'red',label = 'frac_round') + scale_y_log10() + theme_ameres(type = "barplot") + ylab("error rate (percentage)") + xlab("Number of T>C") 
print(p)
dev.off()


##### i want to sample from this dataset and calculate a distribution of the probabilities..  ##### i need all the reads.. in this session i 
    ### have onlt imported the data that have atleast 1 TC 

#### readingin reads that have 0 TC ... I have to get the RA tags of these...
library(data.table)
TC_0 = fread("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/TC_0_RAtags.txt")
#TC_0 = read.table(file = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/TC_0_RAtags.txt")

sampled_1000Reads = vector("list",1000)
names(sampled_1000Reads) = paste0("iter_",c(1:1000))

for(i in 1:1000){
  sampled_1000Reads[[i]]= sample_n(tbl = totalReads,size = 100000,replace = F) %>% filter(T_C == 3) 
  cat(i)  
}


melt(lapply(sampled_1000Reads,nrow)) %>% mutate(fraction = value/100000)


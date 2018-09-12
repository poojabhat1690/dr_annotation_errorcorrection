##### reading in the conversion rates per UTR for all the time points. 
library(dplyr)
library(plyr)


theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15,angle = 90, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}



####### first reading in the count data to get the CPMs of the data... 

countDataSets_quantSeq = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/count/",pattern = "*.tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/count//",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T,stringsAsFactors = F))

CPM_counts = do.call(cbind,lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM))
CPM_counts  = cbind.data.frame(countDataSets_quantSeq_data[[1]][,c(1:6)],CPM_counts)
meanCPMs = apply(CPM_counts[,(7:ncol(CPM_counts))],1,mean)
CPM_counts = CPM_counts[which(meanCPMs>100),]
CPM_counts_cws = CPM_counts[,c(1:6)]

###### making TPMs and count data from SLAMseq samples ... 

countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/utrRates_noBq//",pattern = "mutationrates_utr.csv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/mESCdata_forCheckingUTRerrors/slamdunk_mESC_full_timecourse_1_revision/utrRates_noBq//",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T,stringsAsFactors = F))

countData_timeCourse_datasets = countDataSets_quantSeq_data
splitNames = unlist(lapply(strsplit(countDataSets_quantSeq,"_",T),function(x) x[1]))
#barcodes = unlist(lapply(strsplit(splitNames,".",T),function(x) x[1]))
names(countData_timeCourse_datasets) = splitNames
splitNames = as.data.frame(splitNames)
colnames(splitNames) = "number"

### reading in the sample information file

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/otherstuff/ucscTracks/SLAMseq_revision/quantSeq/timeCourse/summary_short.txt",stringsAsFactors = F)
sampleInfo = sampleInfo[order(sampleInfo$V1),]
sampleInfo$number = unlist(lapply(strsplit(sampleInfo$V1,"_",T),function(x) x[1]))


splitNames = plyr::join(splitNames,sampleInfo)
names(countData_timeCourse_datasets) = splitNames$V2


countDataSets_quantSeq_data = countData_timeCourse_datasets
colnames(CPM_counts_cws) = c("Chr","Start","End","Name","V5","Strand")

countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_data,function(x) left_join(CPM_counts_cws,x))


#countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_data,function(x) x[which(x$ReadCount>500),])

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
      dplyr::select(-one_of("Chr","Start","V5")) %>% 
      
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
    allRates[[i]] =  plotTab
  }
  
  return(allRates)
}



#### correlating TC in the abckground with other mutations... in the grand SLAM paper, the authors used T to A and T to G could explain the variation..

### picking the last time point considered for the pulse ... i.e the 12 hour time point... initially considering only the first replicate..
###### correlating the 

  ### picking all the 12 hr time points... 

countDataSets_quantSeq_data_12hr = countDataSets_quantSeq_data[grep("pulse_12h",names(countDataSets_quantSeq_data))]

#### i need to split the reads from the plus and minus strand... 

countData_12hrPulseCombined = countDataSets_quantSeq_data_12hr[[1]][,c(7:ncol(countDataSets_quantSeq_data_12hr[[1]]))] + countDataSets_quantSeq_data_12hr[[2]][,c(7:ncol(countDataSets_quantSeq_data_12hr[[2]]))] + countDataSets_quantSeq_data_12hr[[3]][,c(7:ncol(countDataSets_quantSeq_data_12hr[[3]]))]
countData_12hrPulseCombined$strand = countDataSets_quantSeq_data_12hr$pulse_12h_R1_34336$Strand

countData_12hrPulseCombined_plus = countData_12hrPulseCombined %>% filter(strand == "+") %>% select(A_A,G_G,C_C,T_T, A_C ,A_G, A_T, C_A, C_G, C_T, G_A, G_C, G_T, T_A,  T_C, T_G)
countData_12hrPulseCombined_minus = countData_12hrPulseCombined %>% filter(strand == "-") %>% select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)

countData_12hrPulseCombined = rbind(countData_12hrPulseCombined_plus,countData_12hrPulseCombined_minus)
samples_log10 = log10(countData_12hrPulseCombined[,c(8:ncol(countData_12hrPulseCombined))])

dir.create(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/mESCdistributions/")

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/mESCdistributions/TA_TG_TC_12hr_NumberOfReads.pdf",height=5,width=5)
ggpubr::ggscatterhist(data =samples_log10,x = 'T_A',y="T_G" ,alpha= 0.5,title = paste(" T_A v/s T_G 12h_R1 (cor = ",round(cor(samples_log10$T_A,samples_log10$T_G,method = "spearman"),2),")","n = ",nrow(samples_log10)),xlab = "log10( Number of T_A reads)", ylab= "log10( Number of T_G reads)",add="reg.line")  + geom_text(aes(label = "cor",x=1,y=1))
ggpubr::ggscatterhist(data =samples_log10,x = 'T_C',y="T_G" ,alpha = 0.5,title = paste(" T_C v/s T_G 12h_R1 (cor = ",round(cor(samples_log10$T_C,samples_log10$T_G,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Number of T_C reads)", ylab= "log10( Number of T_G reads)",add="reg.line")  + geom_text(aes(label = "cor",x=1,y=1))
ggpubr::ggscatterhist(data =samples_log10,x = 'T_C',y="T_G" ,alpha = 0.5,title = paste(" T_C v/s T_A 12h_R1 (cor = ",round(cor(samples_log10$T_C,samples_log10$T_A,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Number of T_C reads)", ylab= "log10( Number of T_A reads)",add="reg.line")  + geom_text(aes(label = "cor",x=1,y=1))
samples_log10 = samples_log10 %>% select(T_A,T_C,T_G)
samples_log10_melt = melt(samples_log10)
ggpubr::ggdensity(data =samples_log10_melt,'value',fill='variable',add="mean") + theme_ameres(type = "barplot")
dev.off()

######### doing the same for the background sample...


countDataSets_quantSeq_data_0hr = countDataSets_quantSeq_data[grep("pulse_0h",names(countDataSets_quantSeq_data))]
countData_0hrPulseCombined = countDataSets_quantSeq_data_0hr[[1]][,c(7:ncol(countDataSets_quantSeq_data_0hr[[1]]))] + countDataSets_quantSeq_data_0hr[[2]][,c(7:ncol(countDataSets_quantSeq_data_0hr[[2]]))] + countDataSets_quantSeq_data_0hr[[3]][,c(7:ncol(countDataSets_quantSeq_data_0hr[[3]]))]

countData_0hrPulseCombined$strand = countDataSets_quantSeq_data_0hr$pulse_0h_R1_34330$Strand

countData_0hrPulseCombined_plus = countData_0hrPulseCombined %>% filter(strand == "+") %>% select(A_A,G_G,C_C,T_T, A_C ,A_G, A_T, C_A, C_G, C_T, G_A, G_C, G_T, T_A,  T_C, T_G)
countData_0hrPulseCombined_minus = countData_0hrPulseCombined %>% filter(strand == "-") %>% select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)

countData_0hrPulseCombined = rbind(countData_0hrPulseCombined_plus,countData_0hrPulseCombined_minus)


samples_log10_0hr = log10(countData_0hrPulseCombined[,c(8:ncol(countData_0hrPulseCombined))])

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/mESCdistributions/TA_TG_TC_0hr_NumberOfReads.pdf",height = 5,width=5)
ggpubr::ggscatterhist(data =samples_log10_0hr,x = 'T_A',y="T_G" ,alpha = 0.5,title = paste(" T_A v/s T_G 0h (cor = ",round(cor(countData_0hrPulseCombined$T_A,countData_0hrPulseCombined$T_G,method = "spearman"),2),")","n = ",nrow(samples_log10)),xlab = "log10( Number of T_A reads)", ylab= "log10( Number of T_G reads)",add="reg.line")  
ggpubr::ggscatterhist(data =samples_log10_0hr,x = 'T_C',y="T_G",alpha = 0.5,title = paste(" T_C v/s T_G 0h (cor = ",round(cor(countData_0hrPulseCombined$T_C,countData_0hrPulseCombined$T_G,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Number of T_C reads)", ylab= "log10( Number of T_G reads)",add="reg.line") 
ggpubr::ggscatterhist(data =samples_log10_0hr,x = 'T_C',y="T_A",alpha = 0.5,title = paste(" T_C v/s T_A 0h (cor = ",round(cor(countData_0hrPulseCombined$T_C,countData_0hrPulseCombined$T_A,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Number of T_C reads)", ylab= "log10( Number of T_A reads)",add="reg.line")  
samples_log10_0hr = samples_log10_0hr %>% select(T_A,T_C,T_G)
samples_log10_0hr_melt = melt(samples_log10_0hr)
ggpubr::ggdensity(data =samples_log10_0hr_melt,'value',fill='variable',add="mean") + theme_ameres(type = "barplot")
dev.off()


### there is also a correlation of the 0 hour samples  between the total number of TA,TC, and TGs...



####################### similarly i would like to plot the fractions of TC conversions..

Tconversions_12hr = countData_12hrPulseCombined%>% select(c("T_T","T_A","T_G","T_C"))
sum_ts = apply(Tconversions_12hr,1,sum)
Tconversions_12hr_fraction= as.data.frame(Tconversions_12hr/sum_ts)

######################## now exactly for this dataset, i want to check the correlation rates of the datasets....

Tconversions_0hr = countData_0hrPulseCombined %>% select("T_A","T_C","T_G","T_T")

sum_ts_0h = apply(Tconversions_0hr,1,sum)

Tconversions_0hr= as.data.frame(Tconversions_0hr/sum_ts_0h)
Tconversions_0hr = Tconversions_0hr
Tconversions_0hr = Tconversions_0hr %>% select("T_A","T_C","T_G")
Tconversions_0hr_melt = melt(Tconversions_0hr)
Tconversions_0hr_melt$value = log10(Tconversions_0hr_melt$value)
Tconversions_0hr_log10 = log10(Tconversions_0hr)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/mESCdistributions/TA_TG_TC_0hr_FractionOfConversions.pdf",height = 5,width=5)
ggpubr::ggscatterhist(data =Tconversions_0hr_log10,x = 'T_A',y="T_G" ,alpha = 0.5,title = paste(" T_A v/s T_G 0h (cor = ",round(cor(Tconversions_0hr_log10$T_A,Tconversions_0hr_log10$T_G,method = "spearman"),2),")","n = ",nrow(samples_log10)),xlab = "log10( fraction of T_As)", ylab= "log10( Fraction of T_Gs)",add="reg.line")  
ggpubr::ggscatterhist(data =Tconversions_0hr_log10,x = 'T_C',y="T_G",alpha = 0.5,title = paste(" T_C v/s T_G 0h (cor = ",round(cor(Tconversions_0hr_log10$T_C,Tconversions_0hr_log10$T_G,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Fraction of T_Cs)", ylab= "log10( Fraction of  T_Gs)",add="reg.line") 
ggpubr::ggscatterhist(data =Tconversions_0hr_log10,x = 'T_C',y="T_A",alpha = 0.5,title = paste(" T_C v/s T_A 0h (cor = ",round(cor(Tconversions_0hr_log10$T_C,Tconversions_0hr_log10$T_A,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Fraction of T_Cs)", ylab= "log10( Fraction of T_As)",add="reg.line")  
ggpubr::ggdensity(data =Tconversions_0hr_melt,'value',fill='variable',add="mean") + theme_ameres(type = "barplot") 
dev.off()

Tconversions_12hr = countData_12hrPulseCombined %>% select("T_A","T_C","T_G","T_T")
sum_ts_12h = apply(Tconversions_12hr,1,sum)

Tconversions_12hr= as.data.frame(Tconversions_12hr)/sum_ts_12h
Tconversions_12hr = Tconversions_12hr
Tconversions_12hr = Tconversions_12hr %>% select("T_A","T_C","T_G")
Tconversions_12hr_melt = melt(Tconversions_12hr)
Tconversions_12hr_melt$value = log10(Tconversions_12hr_melt$value)
Tconversions_12hr_log10 = log10(Tconversions_12hr)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/distributionOfConversions/mESCdistributions/TA_TG_TC_12hr_FractionOfConversions.pdf",height = 5,width=5)
ggpubr::ggscatterhist(data =Tconversions_12hr_log10,x = 'T_A',y="T_G" ,alpha = 0.5,title = paste(" T_A v/s T_G 0h (cor = ",round(cor(Tconversions_12hr_log10$T_A,Tconversions_12hr_log10$T_G,method = "spearman"),2),")","n = ",nrow(samples_log10)),xlab = "log10( fraction of T_As)", ylab= "log10( Fraction of T_Gs)",add="reg.line")  
ggpubr::ggscatterhist(data =Tconversions_12hr_log10,x = 'T_C',y="T_G",alpha = 0.5,title = paste(" T_C v/s T_G 0h (cor = ",round(cor(Tconversions_12hr_log10$T_C,Tconversions_12hr_log10$T_G,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Fraction of T_Cs)", ylab= "log10( Fraction of  T_Gs)",add="reg.line") 
ggpubr::ggscatterhist(data =Tconversions_12hr_log10,x = 'T_C',y="T_A",alpha = 0.5,title = paste(" T_C v/s T_A 0h (cor = ",round(cor(Tconversions_12hr_log10$T_C,Tconversions_12hr_log10$T_A,method = "spearman"),2),")",nrow(samples_log10)),xlab = "log10( Fraction of T_Cs)", ylab= "log10( Fraction of T_As)",add="reg.line")  
ggpubr::ggdensity(data =Tconversions_12hr_melt,'value',fill='variable',add="mean") + theme_ameres(type = "barplot") 
dev.off()


ggpubr::ggdensity(data =Tconversions_12hr_melt,'value',col='variable',add="mean",alpha=0.5,size=1) + theme_ameres(type = "barplot") + scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2")





cdsQS_0h <-
  as.tibble(countDataSets_quantSeq_data_0hr[[1]]) %>%
  select(Chr, Start, End, Name, Strand, ReadCount, A_A, A_C, A_G, A_T, T_A, T_C, T_G, T_T) %>%
  mutate(trueT_A = ifelse(Strand == "+", T_A, A_T),
         trueT_C = ifelse(Strand == "+", T_C, A_G),
         trueT_G = ifelse(Strand == "+", T_G, A_C),
         trueT_T = ifelse(Strand == "+", T_T, A_A),
         Tsum = trueT_A + trueT_G + trueT_C + trueT_T) %>%
  mutate(fracT_A = trueT_A / Tsum,
         fracT_C = trueT_C / Tsum,
         fracT_G = trueT_G / Tsum,
         fracT_T = trueT_T / Tsum) %>%
  select(-matches("A_"), -(T_A:T_T))



# cdsQS_12h %>%
#   ggplot(aes(log10(trueT_A), log10(trueT_G))) +
#   geom_point(alpha = 0.5)

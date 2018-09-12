
library(RColorBrewer)

#### finds the number of reads with TCs in the untreated samples
 
####### it is better to run the first part on the cluster.

# untreatedFiles = list.files(path ="////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/",pattern ="freqAllMutations_allTranscripts*" )
# untreatedFiles = untreatedFiles[grep("Untreated",untreatedFiles)]
# untreatedFiles_path = paste0("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/",untreatedFiles)
# 
# untreatedFiles_data = lapply(untreatedFiles_path,function(x) read.delim(x,stringsAsFactors = F,header = T))
# names(untreatedFiles_data) = untreatedFiles
# 
# getTabulatedReads = function(untreatedSample_1){
#   minusTab = untreatedSample_1 %>% dplyr::filter(strand == "-") %>% 
#     dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,strand,id)
#   plusTab = untreatedSample_1 %>% dplyr::filter(strand == "+") %>%
#     dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,strand,id)
#   
#   totalReads = rbind.data.frame(plusTab,minusTab)
#   totalReads_tabulated  = as.data.frame(table(totalReads$T_C)) %>% mutate(frac_freq = Freq/sum(Freq))
#   return(totalReads_tabulated)
# }
# 
# 
# tabulatedSamples = lapply(untreatedFiles_data,function(x) getTabulatedReads(x))
# 
# #############################################################################################
############ 

#save(list = tabulatedSamples,file = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/Frequency_perTimePoint_numberOfTs.Rdata")
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

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/Frequency_perTimePoint_numberOfTs.Rdata")
tabulatedSamples_timepoints = unlist(lapply(strsplit(names(tabulatedSamples),"_",T) ,function(x) x[3]) )
names(tabulatedSamples) = substr(tabulatedSamples_timepoints,start = 1,3)
nTC = as.data.frame(factor(c(0:15)))
colnames(nTC) = "Var1"
tabulatedSamples = tabulatedSamples[paste0("TP",c(1:9))]
nTC$Freq = 0
library(plyr)
tabulatedSamples$df = nTC
library(purrr)
tabulatedSamples = tabulatedSamples  %>% purrr::reduce(plyr::join, by = "Var1",.init = tabulatedSamples[[10]])

colnames(tabulatedSamples) = paste(colnames(tabulatedSamples),c("ref","ref",rep(c(1:9),each = 2)),sep="_")
tabulatedSamples$Freq_ref = NULL
colsUse = c(brewer.pal(n = 8,name = "Dark2"),brewer.pal(3,name = "Set1"))

options( scipen = 20 )

pdf("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/numberOfTcs_errorRate.pdf",height=4.5)

p = tabulatedSamples %>% dplyr::select(contains("frac_freq")) %>% tidyr::gather() %>% plyr::mutate(nMut = c(0:15))%>%  mutate(time = paste0("TP_",rep(1:9,each=16)))  %>% dplyr::group_by(nMut) %>% dplyr::mutate(medianS = median(value,na.rm = T)) %>% ungroup() %>% 
    filter(nMut != 0)%>% filter(nMut != 1) %>%ggpubr::ggline(data = .,x = 'nMut',y = 'value',group='time',col='time',size=1,alpha=0.5) + scale_color_manual(values = colsUse) +
    theme_ameres(type="barplot") + scale_x_continuous(labels = c(1:15),breaks = c(1:15)) + xlab("number of T_C conversions") + ylab("log10(Fraction of TC reads)")+ scale_y_log10()
print(p)

p = tabulatedSamples %>% dplyr::select(contains("frac_freq")) %>% tidyr::gather() %>% plyr::mutate(nMut = c(0:15))%>%  mutate(time = paste0("TP_",rep(1:9,each=16)))  %>% dplyr::group_by(nMut) %>% dplyr::mutate(medianS = median(value,na.rm = T)) %>% ungroup() %>% 
  filter(nMut != 0)%>% filter(nMut != 1) %>% filter(row_number() <= 14) %>%mutate(frac_round = round(medianS,6)) %>% ggpubr::ggline(data = .,x='nMut',y='medianS',label = 'frac_round',col='red',size=1)  + scale_x_continuous(labels = c(1:15),breaks = c(1:15))+
    theme_ameres(type = "barplot") + scale_y_log10() + ylab("log10(median fraction TC)") + xlab("number of T_C conversions")

print(p)
dev.off()

######## so...I want to estimate the number of reads that will be wrong in the data based on the sequencing depth of the datasets... 


tabulatedSamples %>% dplyr::select(contains("frac_freq")) %>% tidyr::gather() %>% plyr::mutate(nMut = c(0:15)) %>%  mutate(time = paste0("TP_",rep(1:9,each=16)))  %>% dplyr::group_by(nMut) %>% dplyr::mutate(medianS = median(value,na.rm = T)) %>% ungroup() %>% 
  filter(nMut != 0)%>% filter(nMut != 1) %>% filter(row_number() <= 14) 

#### i also want to know the number of reads that have 3 TCs..... per library.... 

      #### i want to read all the data, of number of TCs per samples.....


Rdata_allSamples = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",pattern ="freqAllMutations_allTranscripts*")
Rdata_samplepaths  = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",Rdata_allSamples)

falsePositiveReads = c()
totalNumberOfreadsWithTC = c()
allReads = c()
for(i in 1:length(Rdata_samplepaths)){
  load(Rdata_samplepaths[i])
  a = sum(allMutations$T_C)  * 0.0000619
 falsePositiveReads = c(falsePositiveReads, a)
  totalNumberOfreadsWithTC = c(totalNumberOfreadsWithTC,apply(allMutations$T_C,2,sum)[4])
  allReads = c( allReads,sum(allMutations$T_C))
}

falsePositvies_samples = cbind.data.frame(falsePositiveReads,totalNumberOfreadsWithTC,Rdata_allSamples,allReads)
falsePositvies_samples = falsePositvies_samples  %>% mutate(falsePositiveRate = falsePositiveReads/totalNumberOfreadsWithTC)  %>% filter(grepl("Inj",Rdata_allSamples)) %>% 
  mutate(condition = substr(Rdata_allSamples,start = 42,stop = 47)) %>% mutate(timepoint = substr(Rdata_allSamples,start = 45,stop = 47)) %>% arrange(.,timepoint) %>% mutate(index = c(1:27)) 

  pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/expectedBackgroundVStrue.pdf")
  p  = ggpubr::ggscatter(data = falsePositvies_samples,x = 'index',y = 'falsePositiveRate',col='timepoint') + ggrepel::geom_text_repel(label=falsePositvies_samples$condition) + ylab("Expected bg reads/actual reads")
  p + theme_ameres(type = "barplot") + geom_hline(aes(yintercept = 0.5),col="red",linetype="dashed",size=0.5)
  dev.off()
  
  
  #### i also want to calculate an average of the replicates... 
  detach(package:plyr)
  library(dplyr)
  
  a = falsePositvies_samples %>% dplyr::group_by(timepoint) %>% mutate(sum_TotalReadsWithTC = sum(totalNumberOfreadsWithTC)) %>%
    mutate(sum_allreads = sum(allReads)) %>% mutate(totalFalsePositive = sum(falsePositiveReads) ) %>% mutate(falsePositive_combined = totalFalsePositive/sum_TotalReadsWithTC )%>% distinct(falsePositive_combined) 
  
  pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/expecterackgroundVsTrue_combinedReplicates.pdf")
  p  = ggpubr::ggscatter(data = a,x = 'timepoint',y = 'falsePositive_combined') + ylab("Expected bg reads/actual reads")
  p + theme_ameres(type = "barplot") + geom_hline(aes(yintercept = 0.5),col="red",linetype="dashed",size=0.5)
  dev.off()
  
  
  #####################
  ###### i can also get the specific mutation rates for each dataset... 
  #####################
  
  
  backgoroundErrorPerreplicateInBg = tabulatedSamples %>% filter(Var1_ref == 3) %>% dplyr::select(contains("frac_freq")) 
  falsePositiveReads = c()
  totalNumberOfreadsWithTC = c()
  allReads = c()
  for(i in 1:length(Rdata_samplepaths)){
    load(Rdata_samplepaths[i])
    sampleNum = as.numeric(substr(Rdata_allSamples[1],47,47))
    a = sum(allMutations$T_C)  * backgoroundErrorPerreplicateInBg[,sampleNum]
    falsePositiveReads = c(falsePositiveReads, a)
    totalNumberOfreadsWithTC = c(totalNumberOfreadsWithTC,apply(allMutations$T_C,2,sum)[4])
    allReads = c( allReads,sum(allMutations$T_C))
  }
  
  falsePositvies_samples = cbind.data.frame(falsePositiveReads,totalNumberOfreadsWithTC,Rdata_allSamples,allReads)
  falsePositvies_samples = falsePositvies_samples  %>% mutate(falsePositiveRate = falsePositiveReads/totalNumberOfreadsWithTC)  %>% filter(grepl("Inj",Rdata_allSamples)) %>% 
    mutate(condition = substr(Rdata_allSamples,start = 42,stop = 47)) %>% mutate(timepoint = substr(Rdata_allSamples,start = 45,stop = 47)) %>% arrange(.,timepoint) %>% mutate(index = c(1:27)) 
  a = falsePositvies_samples %>% dplyr::group_by(timepoint) %>% mutate(sum_TotalReadsWithTC = sum(totalNumberOfreadsWithTC)) %>%
    mutate(sum_allreads = sum(allReads)) %>% mutate(totalFalsePositive = sum(falsePositiveReads) ) %>% mutate(falsePositive_combined = totalFalsePositive/sum_TotalReadsWithTC )%>% distinct(falsePositive_combined) 
  
  
  
  pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/expecterackgroundVsTrue_correctedForBg_eachTimepoint.pdf")
  
  p  = ggpubr::ggscatter(data = falsePositvies_samples,x = 'index',y = 'falsePositiveRate',col='timepoint') + ggrepel::geom_text_repel(label=falsePositvies_samples$condition) + ylab("Expected bg reads/actual reads")
  p + theme_ameres(type = "barplot") + geom_hline(aes(yintercept = 0.5),col="red",linetype="dashed",size=0.5)
  
  p  = ggpubr::ggscatter(data = a,x = 'timepoint',y = 'falsePositive_combined') + ylab("Expected bg reads/actual reads")
  p + theme_ameres(type = "barplot") + geom_hline(aes(yintercept = 0.5),col="red",linetype="dashed",size=0.5)
  
  dev.off()  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


#### so R2_TP2 has some genes with Fp<0
r2_TP2 = Rdata_samplepaths[grep("freqAllMutations_allTranscriptsACCGTGInj_R2_TP2.txt.Rdata",Rdata_samplepaths)]
load(r2_TP2)

a = as.data.frame(allMutations$T_C[which(allMutations$T_C[,4] >=1),] )  %>% distinct(.,id_samples)

allHaving3 = allMutations$T_C[which(allMutations$T_C[,4] >=1),] 

### out of these 105 reads are probably wrong....wrong 
allHaving3 = as.data.frame.matrix(allHaving3)  
colnames(allHaving3) = paste0("nMut",colnames(allHaving3))
allHaving3$ids = row.names(allHaving3)

allHaving3 = allHaving3 %>% filter(nMut3 >1) %>% filter(nMut2>1)

table(allHaving3[,4])



### what are these genes 

rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows$id = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1)
gnees_exp = rpms_allCountingWindows[rpms_allCountingWindows$id %in% a$id_samples,]

genes_strongExp = rpms_allCountingWindows[rpms_allCountingWindows$id %in% allHaving3$ids,]


##### now i also want to check the stochasticity of these genes... 

#### i neeed to check if these genes have expression in each stage... i.e >=3TC reads in each sample....






#### readin in all the datasets of the genes that have T

library(dplyr)
library(stringr)
allRdata = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized///numberOfGenes_allTranscripts//",pattern = "*.Rdata")
allRdata_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized///numberOfGenes_allTranscripts//",allRdata)
allRdata_path = allRdata_path[grep("Inj",allRdata_path)]
names_samples = str_extract(allRdata_path,"([^/]+)$" ) %>% str_sub(start = -20)

rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows$id = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1 )

searchForGenes = function(TP,path,genes_exp){
  
  ids_greaterthan2 = vector("list",length(allRdata_path[grep(TP,allRdata_path)]))
      for(i in 1:length(ids_greaterthan2)){
          load(allRdata_path[grep(TP,allRdata_path)][i])
          tctemp = as.data.frame.matrix(allMutations$T_C)
          ids_greaterthan2[[i]] = row.names(tctemp[which(tctemp$`3` > 0),])
  
            }
        ids_greaterthan2_combine = unique(do.call(c,ids_greaterthan2))
        gnees_exp_inTP = ids_greaterthan2_combine[ids_greaterthan2_combine %in% gnees_exp$id] ### these are the genes that have atleas 3 TCss at TP1.
        gnees_exp_inTP = as.data.frame(as.character(gnees_exp_inTP))
        names(gnees_exp_inTP) = "id"
        gnees_exp_inTP$present = 1
        return(gnees_exp_inTP)

        }

inTP1 = searchForGenes(TP = "TP1",path = allRdata_path,genes_exp = gnees_exp)
inTP2 = searchForGenes(TP = "TP2",path = allRdata_path,genes_exp = gnees_exp)
inTP3 = searchForGenes(TP = "TP3",path = allRdata_path,genes_exp = gnees_exp)
inTP4 = searchForGenes(TP = "TP4",path = allRdata_path,genes_exp = gnees_exp)
inTP5 = searchForGenes(TP = "TP5",path = allRdata_path,genes_exp = gnees_exp)
inTP6 = searchForGenes(TP = "TP6",path = allRdata_path,genes_exp = gnees_exp)
inTP7 = searchForGenes(TP = "TP7",path = allRdata_path,genes_exp = gnees_exp)
inTP8 = searchForGenes(TP = "TP8",path = allRdata_path,genes_exp = gnees_exp)
inTP9 = searchForGenes(TP = "TP9",path = allRdata_path,genes_exp = gnees_exp)

reference = as.data.frame(as.character(gnees_exp$id))
colnames(reference) = "id"
m = list(reference,inTP1,inTP2,inTP3,inTP4,inTP5,inTP6,inTP7,inTP8,inTP9)
allStochastic = m %>% purrr::reduce(plyr::join, by = "id")
allStochastic[is.na(allStochastic) ]<-0
colnames(allStochastic) = c("id",paste0("TP",c(1:9)))
row.names(allStochastic) = allStochastic$id
allStochastic$id = NULL
ord <- hclust( dist((allStochastic), method = "euclidean"), method = "ward.D" )$order
allStochastic = allStochastic[ord,]
allStochastic_melt = melt(allStochastic)
allStochastic_melt$nGene = c(1:196)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/stochasticExpression_sampleAboveBackgroundTC.pdf")
p = ggplot(allStochastic_melt,aes(x=variable,y = nGene,fill=value)) + geom_tile() + scale_fill_gradientn(colours = c("white", "black"), values = c(0,0.1,1)) 
print(p)
dev.off()
allStochastic= allStochastic %>% rownames_to_column(var = "id")
#genes_consider = rpms_allCountingWindows[rpms_allCountingWindows$id %in% row.names(allStochastic),] 

allGenes_3TC = plyr::join(allStochastic,rpms_allCountingWindows,by="id")

########################## GO term analysis.... 

library(biomaRt)
library(org.Dr.eg.db)
library(topGO)
ensembl  = useMart(host='dec2017.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',dataset = "drerio_gene_ensembl")

all_genes = getBM(attributes = c("external_gene_name","ensembl_gene_id"),filters = "external_gene_name",values = rpms_allCountingWindows$V4,mart = ensembl)
colnames(all_genes) = c("V4","ensembl_gene_id")

all_genes_join  =plyr::join(all_genes,gnees_exp)
all_genes_join = all_genes_join[!duplicated(all_genes_join),]
head(all_genes)
all_genes_join$category =  0
all_genes_join[!complete.cases(all_genes_join),]$category <- 1



all_genes_category = all_genes_join$category
names(all_genes_category) = all_genes_join$ensembl_gene_id
GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes_category, geneSel = function(p) p <
                0.01, description = "Test", annot = annFUN.org, mapping = "org.Dr.eg.db",
              ID = "Ensembl")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
options(scipen=999)
table_goTerms=   GenTable(GOdata, classicFisher = resultFisher,topNodes = 10)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/GoTermAnalysis_allCws_moreThan2TCreads_allCws.pdf")
p = table_goTerms %>% mutate(log20P = -log10(as.numeric(classicFisher))  ) %>% ggplot() + geom_bar(aes(x=Term,y=log20P),stat="identity",col="black",fill="black") +  coord_flip()  + theme_ameres(type = "barplot") + ylab("-log2(pVal)")
print(p)
dev.off()

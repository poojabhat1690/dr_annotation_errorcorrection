#### cheking the distibution of fraction of reads with TC conversions..
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

### run below chunk on cluster

# TPs = paste0("TP",c(1:9))
# 
# allFiles = list.files("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",pattern = "freqAllMutations_allTranscripts*")
# files_injected = allFiles[grep("Inj",allFiles)]
# files_untreated = allFiles[grep("Unt",allFiles)]
# 
# files_injected_path = paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",files_injected)
# files_untreated_path  = paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",files_untreated)
# 
# 
# ####### 
# 
# getPvals = function(TP,files_untreated_path,files_injected_path){
#   load(files_untreated_path[grep(TP,files_untreated_path)])
#   fractionOfreadsWithGreaterThan2 = as.data.frame.matrix(allMutations$T_C) %>% select(.,c(4:ncol(.))) %>% mutate(rowSums = rowSums(.))
#   totalReads = as.data.frame.matrix(allMutations$T_C)  %>% mutate(sumReads = rowSums(.))
#   
#   ### injected samples 
#   injectedSamples = files_injected_path[grep(TP,files_injected_path)]
#   pVal_vecs = vector("list",length = length(injectedSamples))
#   names(pVal_vecs) = substr(injectedSamples,248,256)
#   
#   for(j in 1:length(injectedSamples)){
#     
#     load(injectedSamples[j])
#     fractionOfreads_test = as.data.frame.matrix(allMutations$T_C) %>% select(.,c(4:ncol(.))) %>% mutate(rowSums = rowSums(.))
#     totalReads_test = as.data.frame.matrix(allMutations$T_C)  %>% mutate(sumReads = rowSums(.))
#     fractionTOTest = fractionOfreads_test$rowSums/totalReads_test$sumReads
#     
#     p_vals_samples = c()
#     
#     for(i in 1:length(fractionTOTest)) {
#       p_vals_samples = c(p_vals_samples,t.test(x = fractionOfreadsWithGreaterThan2$rowSums /totalReads$sumReads,mu = fractionTOTest[i],alternative = "less",conf.level = 0.99)$p.value)
#     }
#     pVal_vecs[[j]] = p_vals_samples
#   }
#   
#   return(pVal_vecs)
# }
# 
# TP1 = getPvals(TP = "TP1",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP1,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP1.Rdata")
# 
# TP2 = getPvals(TP = "TP2",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP2,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP2.Rdata")
# 
# TP3 = getPvals(TP = "TP3",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP3,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP3.Rdata")
# 
# TP4 = getPvals(TP = "TP4",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP4,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP4.Rdata")
# 
# TP5 = getPvals(TP = "TP5",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP5,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP5.Rdata")
# 
# 
# TP6 = getPvals(TP = "TP6",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP6,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP6.Rdata")
# 
# TP7 = getPvals(TP = "TP7",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP7,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP7.Rdata")
# 
# TP8 = getPvals(TP = "TP8",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP8,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP8.Rdata")
# 
# TP9 = getPvals(TP = "TP9",files_untreated_path = files_untreated_path ,files_injected_path = files_injected_path)
# save(TP9,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP9.Rdata")
# 

###############

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP1.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP2.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP3.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP4.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP5.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP6.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP7.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP8.Rdata")
load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/pvals_TP9.Rdata")

######## combining all the data 


TP1_significant = TP1 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
 mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 
  
TP2_significant = TP2 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 

TP3_significant = TP3 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 


TP4_significant = TP4 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 


TP5_significant = TP5 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 


TP6_significant = TP6 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 

TP7_significant = TP7 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 

TP8_significant = TP8 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 

TP9_significant = TP9 %>% lapply(. %>% (function(x) data.frame(x))) %>%  lapply(. %>% (function(x) tibble::rownames_to_column(x) ))%>% plyr::join_all(.,by = "rowname",match = "all") %>%
  stats::setNames(c("rowname", "pvals1", "pvals2","pvals3")) %>%  replace(., is.na(.), 1) %>% mutate_at (.vars = vars(pvals1,pvals2,pvals3) ,.funs = funs(numSignificant = . < 0.01)) %>%
  mutate(sumSig = pvals1_numSignificant+pvals2_numSignificant+pvals3_numSignificant) %>% filter(sumSig >1) 


allSignificantTranscripts = list(TP1_significant,TP2_significant,TP3_significant,TP4_significant,TP5_significant,TP6_significant,TP7_significant,TP8_significant,TP9_significant)
names(allSignificantTranscripts) = paste0("TP",c(1:9))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs//numberOfSignificantGenes.pdf",height = 4.5,width=5)
melt(lapply(allSignificantTranscripts,nrow)) %>% ggpubr::ggline(.,x = 'L1',y = 'value',label = 'value') + ylab("number of transcripts significant in 2 replicates") + theme_ameres(type = "barplot") + xlab("timepoint")
melt(lapply(allSignificantTranscripts,function(x) length(which(x$sumSig == 3)))) %>% ggpubr::ggline(.,x = 'L1',y = 'value',label = 'value') + ylab("number of transcripts significant in 3 replicates") + theme_ameres(type = "barplot") + xlab("timepoint")
dev.off()



allSignificantTranscripts_3replicates = lapply(allSignificantTranscripts,function(x) x[(which(x$sumSig == 3)),])



#### checking these in the untreated samples


files_untreated = c("freqAllMutations_allTranscriptsCAGCGTUntreated_TP1.txt.Rdata","freqAllMutations_allTranscriptsGATCACUntreated_TP2.txt.Rdata","freqAllMutations_allTranscriptsACCAGTUntreated_TP3.txt.Rdata","freqAllMutations_allTranscriptsTGCACGUntreated_TP4.txt.Rdata", "freqAllMutations_allTranscriptsACATTAUntreated_TP5.txt.Rdata", "freqAllMutations_allTranscriptsGTGTAGUntreated_TP6.txt.Rdata", "freqAllMutations_allTranscriptsCTAGTCUntreated_TP7.txt.Rdata","freqAllMutations_allTranscriptsTGTGCAUntreated_TP8.txt.Rdata","freqAllMutations_allTranscriptsTCAGGAUntreated_TP9.txt.Rdata")
files_untreated_path  = paste0("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",files_untreated)
siginficantGenes_inBackground = c()
  for(i in 1:length(files_untreated_path)){
    load(files_untreated_path[i])
    fractionOfreadsWithGreaterThan2 = as.data.frame.matrix(allMutations$T_C) %>% dplyr::select(.,c(4:ncol(.))) %>% dplyr::mutate(rowSums = rowSums(.))
    totalReads = as.data.frame.matrix(allMutations$T_C)  %>% mutate(sumReads = rowSums(.))
    confidece_int = t.test(fractionOfreadsWithGreaterThan2$rowSums/totalReads$sumReads,conf.level = 0.99,alternative = "less")$conf.int
    singificantInBG=  which((fractionOfreadsWithGreaterThan2$rowSums/totalReads$sumReads) > confidece_int[2])
    singificantInBG = data.frame(rowname = names(singificantInBG))
   
    siginficantGenes_inBackground =  c(siginficantGenes_inBackground,length(intersect(allSignificantTranscripts_3replicates[[i]]$rowname,singificantInBG$rowname)))
     }


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs//numberOfSignificantGenes_inBg.pdf",height = 4.5,width=5)

significanInBg = data.frame(nGenes=  siginficantGenes_inBackground,TP = paste0("TP",c(1:9)))
ggpubr::ggline(data = significanInBg,x = 'TP',y = 'nGenes',label = 'nGenes')
dev.off()

############## checking the GO terms of the genes that are expressed at early time points... 


library(biomaRt)
library(org.Dr.eg.db)
library(topGO)
ensembl  = useMart(host='dec2017.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',dataset = "drerio_gene_ensembl")
rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
rpms_allCountingWindows$rowname = paste0(rpms_allCountingWindows$V1,":",rpms_allCountingWindows$V2+1)
transcriptsSignificantIn2replicates = lapply(allSignificantTranscripts,function(x) x[which(x$sumSig >= 2),])

 transcriptInfo = lapply(transcriptsSignificantIn2replicates,function(x) plyr::join(x,rpms_allCountingWindows) %>% distinct() )

 
 
all_genes = getBM(attributes = c("external_gene_name","ensembl_gene_id"),filters = "external_gene_name",values = rpms_allCountingWindows$V4,mart = ensembl)
colnames(all_genes) = c("V4","ensembl_gene_id")

all_genes_join  = lapply(transcriptInfo,function(x) plyr::join(all_genes,x))
all_genes_join = lapply(all_genes_join,function(x) x[!duplicated(x),])
all_genes_category = vector("list",9)
GOdata = vector("list",9)

for(i in 1:length(all_genes_join)){
  all_genes_join[[i]]$categoty = 0
  all_genes_join[[i]][!complete.cases(all_genes_join[[i]]),]$categoty <- 1
  all_genes_category[[i]] = all_genes_join[[i]]$categoty
  names(all_genes_category[[i]]) = all_genes_join[[i]]$ensembl_gene_id
  GOdata[[i]] <- new("topGOdata", ontology = "BP", allGenes = all_genes_category[[i]], geneSel = function(p) p <
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dr.eg.db",
                ID = "Ensembl")
  cat(i)
}  


resultFisher <- lapply(GOdata,function(x) runTest(x, algorithm = "classic", statistic = "fisher"))
options(scipen=999)
table_goTerms = vector("list",9)

for(i in 1:9){
  table_goTerms[[i]] = GenTable(GOdata[[i]], classicFisher = resultFisher[[i]],topNodes = 10)
}


##### i also want to compare this to the list of only zygotic transcripts.. 


##### these are zygotic transcripts above 5RPM (keep in mind!!!!!)

zygoticGenes = read.table("///Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)


zygotic_goTerms = plyr::join(all_genes,zygoticGenes)
zygotic_goTerms$categoty = 0
zygotic_goTerms[!complete.cases(zygotic_goTerms),]$categoty <- 1
zygotic_goTerms_category  = zygotic_goTerms$categoty
names(zygotic_goTerms_category) = zygotic_goTerms$ensembl_gene_id
GOdata_zygotic <- new("topGOdata", ontology = "BP", allGenes = zygotic_goTerms_category, geneSel = function(p) p <
                     0.01, description = "Test", annot = annFUN.org, mapping = "org.Dr.eg.db",
                   ID = "Ensembl")
fisher_zygotic = runTest(GOdata_zygotic, algorithm = "classic", statistic = "fisher")
Goterms_zygotic = GenTable(GOdata_zygotic, classicFisher = fisher_zygotic,topNodes = 10)
Goterms_zygotic$classicFisher = as.numeric(gsub("< ", "", Goterms_zygotic$classicFisher))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/GOterms_significantIn2replicated.pdf")
Goterms_zygotic %>% mutate(log20P = -log10(as.numeric(classicFisher))  ) %>%
  ggplot() + geom_bar(aes(x=Term,y=log20P),stat="identity",col="black",fill="black") +  
  coord_flip()  + theme_ameres(type = "barplot") + ylab("-log2(pVal)") + ggtitle("ZygoticTranscripts")

### also plotting the go terms of genes that have >3TC reads in 2 replicates


for(i in 1:length(table_goTerms)){
  table_goTerms[[i]]$classicFisher = as.numeric(gsub("<", "", table_goTerms[[i]]$classicFisher))
  p = table_goTerms[[i]] %>% mutate(log20P = -log10(as.numeric(classicFisher))  ) %>%
    ggplot() + geom_bar(aes(x=Term,y=log20P),stat="identity",col="black",fill="black") +  
    coord_flip()  + theme_ameres(type = "barplot") + ylab("-log2(pVal)") + ggtitle(paste0("TP",i))
  print(p)
}

dev.off()


# 
# ###### i also want to know how many of such cases exist in the background... 
# 
# allFiles = list.files("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",pattern = "freqAllMutations_allTranscripts*")
# 
# files_untreated = allFiles[grep("Unt",allFiles)]
#  
# files_untreated_path  = paste0("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/numberOfGenes_allTranscripts/",files_untreated)
#  
#   for(i in 1:length(files_untreated_path)){
#     load(files_untreated_path[i])  
#     fractionOfreadsWithGreaterThan2 = as.data.frame.matrix(allMutations$T_C) %>% dplyr::select(.,c(4:ncol(.))) %>% dplyr::mutate(rowSums = rowSums(.))
#     totalReads = as.data.frame.matrix(allMutations$T_C)  %>% mutate(sumReads = rowSums(.))
#     confidece_int = t.test(fractionOfreadsWithGreaterThan2$rowSums/totalReads$sumReads,conf.level = 0.99,alternative = "less")$conf.int
#     singificantInBG=  which((fractionOfreadsWithGreaterThan2$rowSums/totalReads$sumReads) > confidece_int[2])
#     singificantInBG = data.frame(rowname = names(singificantInBG))
#     singificantInBG = plyr::join(singificantInBG,rpms_allCountingWindows)  %>% distinct()
#     BG_goTerms = plyr::join(all_genes,singificantInBG)
#     BG_goTerms$categoty = 0
#     BG_goTerms[!complete.cases(BG_goTerms),]$categoty <- 1
#     BG_goTerms_category  = BG_goTerms$categoty
#     names(BG_goTerms_category) = BG_goTerms$ensembl_gene_id
#     GOdata_untreated <- new("topGOdata", ontology = "BP", allGenes = BG_goTerms_category, geneSel = function(p) p <
#                             0.01, description = "Test", annot = annFUN.org, mapping = "org.Dr.eg.db",
#                           ID = "Ensembl")
#     fisher_untreated = runTest(GOdata_untreated, algorithm = "classic", statistic = "fisher")
#     Goterms_untreated = GenTable(GOdata_untreated, classicFisher = fisher_untreated,topNodes = 10)
#     Goterms_untreated$classicFisher = as.numeric(gsub("< ", "", Goterms_untreated$classicFisher))
#     
#     
#   }
# 
#   

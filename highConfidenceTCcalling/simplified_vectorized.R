
library(Rsamtools)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)

dir.create("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/")


vectorizedsample = function(vcfFile,BamFile,countingWindows,barcode,condition,multiVCF){
  
  vcfFile = read.table(vcfFile)
  VcfFie_multi = read.table(multiVCF,header=T)
  countingWindows_zygotic_plus = countingWindows_zygotic %>% filter(V6 == "+")
  countingWindows_zygotic_minus = countingWindows_zygotic %>% filter(V6 == "-")
  countingWindows_zygotic_plus = makeGRangesFromDataFrame(df = countingWindows_zygotic_plus,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
  countingWindows_zygotic_minus = makeGRangesFromDataFrame(df = countingWindows_zygotic_minus,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
  
  bamFile <- BamFile(file = BamFile)
  open(bamFile)
  params <- ScanBamParam(which  = countingWindows_zygotic_plus, what = scanBamWhat(),tag=c("RA","MP"),flag=scanBamFlag(isMinusStrand=FALSE)) #### set the parameters required
  
  aln <- scanBam(bamFile, param = params) ### scan the bam
  
  params <- ScanBamParam(which  = countingWindows_zygotic_minus, what = scanBamWhat(),tag=c("RA","MP"),flag=scanBamFlag(isMinusStrand=T)) #### set the parameters required
  bamFile <- BamFile(file = BamFile)  #### assign and open bam file
  aln_minus <- scanBam(bamFile, param = params) ### scan the bam
  cat("this is okay, aln_plus and aln_minus are created")
  aln = c(aln,aln_minus)

    numberReads = melt(lapply(aln,function(x) length(x$pos)))
  
  aln[which(numberReads<1)] <- NULL #### just excluding genes for which there are no reads...
  vcfFile_mutn = paste(vcfFile$V1,(vcfFile$V2+1),sep="_")
  vcfFile_muth_full = paste(VcfFie_multi$Chrom,(VcfFie_multi$Position+1),sep="_")
  vcfFile_mutn = c(vcfFile_mutn,vcfFile_muth_full)
  vcfFile_mutn = vcfFile_mutn[!duplicated(vcfFile_mutn)]
  #cat("length=",vcfFile_muth_full)
  #### some counting windows have reads but have 0 mis matches... I want to add an NA there... 
  
  RAtag_num = melt(lapply(aln,function(x) length(x$tag$RA))) ### i need the RA tag here to add the exact number of MP tags 
  MPtag_num = melt(lapply(aln,function(x) length(x$tag$MP)))
  RA_MP_num = cbind.data.frame(RAtag_num,MPtag_num)
  noMPtags = which(RA_MP_num[,3] == 0)
  
  for(MPtags in 1:length(noMPtags)){
    
    aln[[noMPtags[MPtags]]]$tag$MP = rep(NA, RA_MP_num[noMPtags[MPtags],1])
  }
  cat("this has worked")
  
  aln_chr = as.data.frame(unlist(lapply(aln,function(x) x$rname)))
  aln_RA = as.data.frame(unlist(lapply(aln,function(x) x$tag$RA)))
  aln_MP = as.data.frame(unlist(lapply(aln,function(x) x$tag$MP)))
  aln_position = as.data.frame(unlist(lapply(aln,function(x) x$pos)))
  aln_strand = as.data.frame(unlist(lapply(aln,function(x) x$strand)))
  aln_phred = as.data.frame(unlist(lapply(aln,function(x) as.data.frame(x$qual))))
  alignmentMatrix = cbind.data.frame(aln_chr,aln_RA,aln_MP,aln_position,aln_strand,aln_phred)
  colnames(alignmentMatrix) = c("chr","RA","MP","startPos","strand","qual")
  alignmentMatrix$MP = as.character(alignmentMatrix$MP)
  alignmentMatrix$id = names(unlist(lapply(aln,function(x) x$rname)))
  
  ###### i want to remove the reads that have no mutations... 
  ###### i will add these back at a later stage... 
  
  readsWith0Mutations = alignmentMatrix[!complete.cases(alignmentMatrix),]
  alignmentMatrix =  alignmentMatrix[complete.cases(alignmentMatrix),]
  
  
  
  library(stringr)
  library(tidyr)
  maxCount =  max(str_count(alignmentMatrix$MP,",") + 1)
  MPtags_separated  = separate(data = alignmentMatrix,col = MP,sep = "\\," ,into = paste0("mutNum",c(1:maxCount)))  %>% select(matches("mutNum|start"))
  
  MPtags_separated = MPtags_separated %>%
    mutate_at(.funs = funs(mutId = as.numeric(str_extract(.,pattern = "[^:]+"))), .vars = vars(mutNum1:paste0("mutNum",maxCount)))
  
  
  #str_extract(string = MPtags_separated$mutNum1,pattern = "[:^]+")
  
  
  ##################### identifying low quality bases
  
  ##################### the read position is used to get rid of low quality bases
  
  
  MPtags_separated_readPos = MPtags_separated %>% mutate_at(.funs = funs(readPos =  as.numeric(str_sub(str_extract(.,"\\:([0-9]+)\\:"), 2, -2))),.vars= vars(mutNum1:paste0("mutNum",maxCount)) ) 
  phredScores = alignmentMatrix %>% select(qual) 
  MPtags_separated_readPos = bind_cols(phredScores,MPtags_separated_readPos)
  #### i need to now substr the phredScores
  
  MPtags_separated_readPos = MPtags_separated_readPos %>% mutate_at(.funs = funs(baseMutated = str_sub(qual,.,.)),.vars = vars(mutNum1_readPos:paste0("mutNum",maxCount,"_readPos")))  
  
  
  convertToNumeric = function(col){
    intScore=mapply(function(x) (utf8ToInt(x)-33) < 27, col)
    return(intScore)
  }
  
  PhredScores_numeric = as.data.frame(t(apply(MPtags_separated_readPos %>% select(contains("baseMutated")),1, function(x) convertToNumeric(x)) ))
  MPtags_separated_readPos = MPtags_separated_readPos %>% select(-contains("baseMutated")) %>% select(-contains("qual")) 
  MPtags_separated_readPos = bind_cols(MPtags_separated_readPos,PhredScores_numeric)
  
  ########################################################## So now I have the per position quality for each of the mutations... 
  
  
  
  ######### i also need to check if there is a SNP or not...for this
  ### extract thr ref position,
  ### add this to the start positon
  ### check if this is a SNP as detected in the VCF file...
  
  MPtags_separated = bind_cols(MPtags_separated,alignmentMatrix %>% select("chr"))
  
  SNPCheck = MPtags_separated %>% mutate_at(.funs = funs(refPos = paste0(chr,"_",as.numeric(str_extract(string = .,pattern = "([0-9]+)$")) + startPos ) %in% vcfFile_mutn),.vars = vars(mutNum1:paste0("mutNum",maxCount)))
  
  ####### so now I know if every mutation is a SNP or is a low qulaity base 
  
  SNPCheck_binary = SNPCheck %>% select(contains("refPos"))
  MPtags_binary = MPtags_separated_readPos %>% select(contains("baseMutated"))
  removeBase = as.matrix(SNPCheck_binary+ MPtags_binary )
  removeBase[which(removeBase >1)] <-1
  removeBase[which(removeBase ==0)] <- NA ### this base does not have to be subtracted --> true mutations
  mutationTYpe = as.matrix( MPtags_separated %>% select(contains ("mutId")) )
  removeBase = as.data.frame(removeBase)
  #removeBase = removeBase == 1
  
  idOfBaseRemoval = as.data.frame(removeBase * mutationTYpe)
  idOfBaseRemoval = idOfBaseRemoval + 1 
  
  RATag_split = separate(data = alignmentMatrix %>% select("RA"),col = "RA",sep = "\\,", into =   c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N"))%>% mutate_all(.funs = funs(as.numeric))
  
  RAtag_mutRemovo = bind_cols(RATag_split,idOfBaseRemoval)
  
  library(purrr)
  
  ## now i go over each of the mutation types and change the type... i need to do this for the n mutations in each read... 
  
  
  
  
  nMutRemove = RAtag_mutRemovo %>% mutate(n_A_A = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 1)))) %>%
    mutate(n_A_C = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 2)))) %>%
    mutate(n_A_G = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 3)))) %>%
    mutate(n_A_T = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 4)))) %>%
    mutate(n_C_A = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 6)))) %>%
    mutate(n_C_C = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 7)))) %>%
    mutate(n_C_G = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 8)))) %>%
    mutate(n_C_T = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 9)))) %>%
    mutate(n_G_A = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 11)))) %>%
    mutate(n_G_C = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 12)))) %>%
    mutate(n_G_G = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 13)))) %>%
    mutate(n_G_T = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 14)))) %>%
    mutate(n_T_A = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 16)))) %>%
    mutate(n_T_C = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 17)))) %>%
    mutate(n_T_G = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 18)))) %>% 
    mutate(n_T_T = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 18)))) %>%  select(contains("n_",ignore.case = F))
  
  RAtags_final = RAtag_mutRemovo %>% select(c(A_A,A_C,A_G,A_T,C_A,C_C,C_G,C_T,G_A,G_C,G_G,G_T,T_A,T_C,T_G,T_T))
  RAtags_final = RAtags_final - nMutRemove
  RAtags_final$id = alignmentMatrix$id
  RAtags_final$strand = alignmentMatrix$strand
  
  #### i also need to add back the reads that do not have any mutations... 
  
  readsWith0Mutations_subset = separate(data = readsWith0Mutations,col = "RA",sep = "\\,", into =   c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N"))  %>%   select(c(A_A,A_C,A_G,A_T,C_A,C_C,C_G,C_T,G_A,G_C,G_G,G_T,T_A,T_C,T_G,T_T)) %>% mutate_all(.funs = funs(as.numeric))
  readsWith0Mutations_subset$id = readsWith0Mutations$id
  readsWith0Mutations_subset$strand = readsWith0Mutations$strand
  
  finalMatrix = bind_rows(RAtags_final,readsWith0Mutations_subset)
  write.table(finalMatrix,paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/freqAllMutations_allTranscripts_totalSNPremoval",barcode,condition,".txt"),sep="\t",quote = F)
}




### from the RPMS of the data 
#zygoticGenes = read.table("/////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)

rpms_allCountingWindows = read.table("//////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)

#countingWindows_zygotic = rpms_allCountingWindows[rpms_allCountingWindows$V4 %in% zygoticGenes$external_gene_name,]
#rpms_allCountingWindows_mt = rpms_allCountingWindows[which(rpms_allCountingWindows$V1 == "chrM"),]

##### change this
rpms_allCountingWindows = rpms_allCountingWindows %>% mutate( sums = rowSums(rpms_allCountingWindows %>% select(-contains("V")))  ) %>% filter(sums != 0) %>% select(contains("V"))



# countingWindows_zygotic = rpms_allCountingWindows %>% mutate(mean_tp9 = rowMeans(data.frame(Inj_R1_TP9,Inj_R2_TP9,Inj_R3_TP9))) %>% filter(mean_tp9 >= 5) %>% select(V1:V6)
# rpms_allCountingWindows_mt = rpms_allCountingWindows_mt %>% select(1:6)

countingWindows_zygotic = rpms_allCountingWindows

countingWindows_zygotic_gr = makeGRangesFromDataFrame(df = countingWindows_zygotic,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)



##### what are the barcodes for which i want this info: 

sampleInfo = read.table("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",stringsAsFactors = F)
sampleInfo = sampleInfo %>% filter(grepl(paste(c("Inj_"), collapse="|"), V3))

barcodesOfInterest = sampleInfo$V2
conditionsOfInterest = sampleInfo$V3
infoSplit = tidyr::separate(data = sampleInfo,col = V3,sep="_",into=c('condn','replicate','timepoint'))
#countingWindows_zygotic_gr = countingWindows_zygotic_gr[1:100] 

for(i in 1:length(barcodesOfInterest)){
  vcfFile_ofInterest = paste0("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/snp/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf")
  BamFile_ofInterest = paste0("/////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam")
  infoSplit_interest =  paste0("/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/varSample_",infoSplit$replicate[i],".vcf")
  vectorizedsample(vcfFile = vcfFile_ofInterest,BamFile = BamFile_ofInterest,countingWindows = countingWindows_zygotic,barcode = barcodesOfInterest[i],condition = conditionsOfInterest[i],multiVCF = infoSplit_interest)
}


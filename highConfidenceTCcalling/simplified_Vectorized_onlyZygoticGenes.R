
library(Rsamtools)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)

################################################################################

dir.create("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/") 


vectorizedsample = function(vcfFile,BamFile,countingWindows,barcode,condition){ 
  
  #### vcfFile - path to the vcf file.
  #### BamFile - path to the bam file. 
  #### countingWindows - gRanges object of the regions you want to inspect
  #### barcode info - to know the identity of the sample as i run this currently in a for loop
  #### condition - only required for the name while writing out the table...
  
  vcfFile = read.table(vcfFile)  #### read in the vcf file
  bamFile <- BamFile(file = BamFile)  #### assign and open bam file
  open(bamFile)
  params <- ScanBamParam(which  = countingWindows, what = scanBamWhat(),tag=c("RA","MP")) #### set the parameters required
  
  aln <- scanBam(bamFile, param = params) ### scan the bam
  ######some of the counting windows do not have any reads, i want to know which they are and i want to remove these...################################################################################
  
  numberReads = melt(lapply(aln,function(x) length(x$pos))) 
  aln[which(numberReads<1)] <- NULL #### just excluding genes for which there are no reads...
  
  ################################################################################################################################################################
  
  
  
  #### some counting windows have reads but have 0 mis matches...I still want to consider these counting windows at this point..,so adding an NA for thesse #### 
  
  RAtag_num = melt(lapply(aln,function(x) length(x$tag$RA)))
  MPtag_num = melt(lapply(aln,function(x) length(x$tag$MP)))
  RA_MP_num = cbind.data.frame(RAtag_num,MPtag_num)
  noMPtags = which(RA_MP_num[,3] == 0)
  
  for(MPtags in 1:length(noMPtags)){
    
    aln[[noMPtags[MPtags]]]$tag$MP = rep(NA, RA_MP_num[noMPtags[MPtags],1]) #### just adding one NA would work as well... 
  }
  
  ####### getting the relavant information from the alignment file as a table...
  
  aln_chr = as.data.frame(unlist(lapply(aln,function(x) x$rname)))
  aln_RA = as.data.frame(unlist(lapply(aln,function(x) x$tag$RA)))
  aln_MP = as.data.frame(unlist(lapply(aln,function(x) x$tag$MP)))
  aln_position = as.data.frame(unlist(lapply(aln,function(x) x$pos)))
  aln_strand = as.data.frame(unlist(lapply(aln,function(x) x$strand)))
  aln_phred = as.data.frame(unlist(lapply(aln,function(x) as.data.frame(x$qual))))
  alignmentMatrix = cbind.data.frame(aln_chr,aln_RA,aln_MP,aln_position,aln_strand,aln_phred)
  colnames(alignmentMatrix) = c("chr","RA","MP","startPos","strand","qual")
  alignmentMatrix$MP = as.character(alignmentMatrix$MP)
  alignmentMatrix$id = row.names(alignmentMatrix)
  
  ###### i want to remove the reads that have no mutations... - and assign it to readsWith0Mutations
  ###### i will add these back at a later stage...
  
  readsWith0Mutations = alignmentMatrix[!complete.cases(alignmentMatrix),]
  alignmentMatrix =  alignmentMatrix[complete.cases(alignmentMatrix),]
  

  library(stringr)
  library(tidyr)
  ########################## spearating the number of mutations that each read contains using the MP tag #############################################################
  
  
  maxCount =  max(str_count(alignmentMatrix$MP,",") + 1) #### getting the maximum number of mutations... str_count conunts the comma, so the number of mutations is 1+number of commas
  
  MPtags_separated  = separate(data = alignmentMatrix,col = MP,sep = "\\," ,into = paste0("mutNum",c(1:maxCount)))  %>% select(matches("mutNum|start")) ###### separating the different mutations per read into columns
                                                                                                                                                        ###### selecting and working with only few columns  
  
  MPtags_separated = MPtags_separated %>%
    mutate_at(.funs = funs(mutId = as.numeric(str_extract(.,pattern = "[^:]+"))), .vars = vars(mutNum1:paste0("mutNum",maxCount))) ######### getting the id of the mutation - position 1 of the RA tag.
  
  ############################################################# ############################################################# #############################################################
  
  
  ##################### identifying low quality bases #######################################
  ##### for this i need to know the read position, and get the same position grom the quality score #######
  
  MPtags_separated_readPos = MPtags_separated %>% mutate_at(.funs = funs(readPos =  as.numeric(str_sub(str_extract(.,"\\:([0-9]+)\\:"), 2, -2))),.vars= vars(mutNum1:paste0("mutNum",maxCount)) ) ### getting the read position for each of the mutations
  phredScores = alignmentMatrix %>% select(qual)  ##### getting the phredScores from the aligbment file...
  MPtags_separated_readPos = bind_cols(phredScores,MPtags_separated_readPos) ##### adding the phred scores to the MP tag info 
  
  #### i need to now substr the phredScores
  
  MPtags_separated_readPos = MPtags_separated_readPos %>% mutate_at(.funs = funs(baseMutated = str_sub(qual,.,.)),.vars = vars(mutNum1_readPos:paste0("mutNum",maxCount,"_readPos")))  ### subsetting the qulaity vector according to the position of the SNP.
  
  
  convertToNumeric = function(col){
    intScore=mapply(function(x) (utf8ToInt(x)-33) < 27, col) ### this gives me info if the the bases is of good quality or not
    return(intScore)
  }
  
  PhredScores_numeric = as.data.frame(t(apply(MPtags_separated_readPos %>% select(contains("baseMutated")),1, function(x) convertToNumeric(x)) )) ### converting the qual score into numeric .. the function utf8ToInt only works on single values and not a character vector!
  MPtags_separated_readPos = MPtags_separated_readPos %>% select(-contains("baseMutated")) %>% select(-contains("qual"))  ### just removing un-necessary previous columnes
  MPtags_separated_readPos = bind_cols(MPtags_separated_readPos,PhredScores_numeric)
  
  #### So now I have the per position quality for each of the mutations... 
  
  ########################################################## i also need to check if there is a SNP or not...for this  ########################################################## 
  
  ### extract thr ref position,
  ### add this to the start positon
  ### check if this is a SNP as detected in the VCF file...
  
  vcfFile_mutn = paste(vcfFile$V1,(vcfFile$V2+1),sep="_") ### the vcf file is 0 based... 
  
  ### the read position is 1 based... 
  
  ## https://www.biostars.org/p/51504/ -  if you access a BAM file via a tool like samtools that turns BAM into SAM then it will be turned into a 1 based format.
  
  MPtags_separated = bind_cols(MPtags_separated,alignmentMatrix %>% select("chr")) ### adding the chromosome information from the the sample...
  
  ### adding start position to the read position and  checking if this is in the vcf file..
  
  SNPCheck = MPtags_separated %>% mutate_at(.funs = funs(refPos = paste0(chr,"_",as.numeric(str_extract(string = .,pattern = "([0-9]+)$")) + startPos ) %in% vcfFile_mutn),.vars = vars(mutNum1:paste0("mutNum",maxCount)))
  
  ####### so now I know if every mutation is a SNP or is a low qulaity base   ##########################################################   ########################################################## 
  
  #### there is now a probalem that a base may be a SNP and that it could be a low quality base... 
  
  SNPCheck_binary = SNPCheck %>% select(contains("refPos")) ### only getting the position of thie bases with SNPS
  MPtags_binary = MPtags_separated_readPos %>% select(contains("baseMutated")) ### getting the bases that have low quality
  removeBase = as.matrix(SNPCheck_binary+ MPtags_binary ) ### just adding these --, if a base has both SNP and lowBq --> sum = 2, if only lowBq or SNP --> sum = 1, of mutation --> sum = 1
  removeBase[which(removeBase >1)] <-1 ###juts assigning everything that is more than 1 to 1 ... so this is an indication of a base that should not be considered..
  removeBase[which(removeBase ==0)] <- NA ### setting the bases that do not have to be considered to NA --> i.e the good bases  
  mutationTYpe = as.matrix( MPtags_separated %>% select(contains ("mutId")) ) ### the information of the type of mutations -- > i want to figure out which of these have to be removed
  removeBase = as.data.frame(removeBase)
  #removeBase = removeBase == 1
  
  idOfBaseRemoval = as.data.frame(removeBase * mutationTYpe) #### multiplying the information of the bases to be removed (1s or NA) by the mutation type gives me the identity of SNPs that have to be removed
  idOfBaseRemoval = idOfBaseRemoval + 1  #### the MP tag position start from 0-24, but the RA tag is a vector, accessed from 1-25...so adding 1
  
  ############ splitting the RA tag....  ############   ############   ############   ############   ############   ############ 
 
  RATag_split = separate(data = alignmentMatrix %>% select("RA"),col = "RA",sep = "\\,", into =   c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N"))%>% mutate_all(.funs = funs(as.numeric))
  RAtag_mutRemovo = bind_cols(RATag_split,idOfBaseRemoval)
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
    mutate(n_T_T = apply(RAtag_mutRemovo %>% select(contains("mutNum")),1, function(x) length(which(x == 19)))) %>%  select(contains("n_",ignore.case = F)) ### checking per row how many of each mutation there are
                                                                                                                  #### selecting only the number of table..
                                                                                                                  #### gives me the number and type of 'bad bases' per read 
  RAtags_final = RAtag_mutRemovo %>% select(c(A_A,A_C,A_G,A_T,C_A,C_C,C_G,C_T,G_A,G_C,G_G,G_T,T_A,T_C,T_G,T_T)) ### seleting relavant columns
  RAtags_final = RAtags_final - nMutRemove ##### subrtacting the number of bad bases per read from the number of mutatiosn (in the RA tag)
  RAtags_final$id = alignmentMatrix$id ### adding other in0
  RAtags_final$strand = alignmentMatrix$strand
  
  ##############################################################################################################################################################################################
  
  
  
  
  #### i also need to add back the reads that do not have any mutations... ##############################################################################################################
  
  readsWith0Mutations_subset = separate(data = readsWith0Mutations,col = "RA",sep = "\\,", into =   c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N"))  %>%   select(c(A_A,A_C,A_G,A_T,C_A,C_C,C_G,C_T,G_A,G_C,G_G,G_T,T_A,T_C,T_G,T_T)) %>% mutate_all(.funs = funs(as.numeric))
  readsWith0Mutations_subset$id = readsWith0Mutations$id
  readsWith0Mutations_subset$strand = readsWith0Mutations$strand
  
  finalMatrix = bind_rows(RAtags_final,readsWith0Mutations_subset)
  write.table(finalMatrix,paste0("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/freqAllMutations_zygoticTranscripts",barcode,condition,".txt"),sep="\t",quote = F)

  #### end of function...
  
  }




### from the RPMS of the data 
zygoticGenes = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)

rpms_allCountingWindows = read.table("//////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)

countingWindows_zygotic = rpms_allCountingWindows[rpms_allCountingWindows$V4 %in% zygoticGenes$external_gene_name,]
rpms_allCountingWindows_mt = rpms_allCountingWindows[which(rpms_allCountingWindows$V1 == "chrM"),]
#countingWindows_zygotic = rpms_allCountingWindows %>% mutate(mean_tp9 = rowMeans(data.frame(Inj_R1_TP9,Inj_R2_TP9,Inj_R3_TP9))) %>% filter(mean_tp9 >= 5) %>% select(V1:V6)
rpms_allCountingWindows_mt = rpms_allCountingWindows_mt %>% select(1:6)
countingWindows_zygotic =countingWindows_zygotic %>%  select(1:6)
countingWindows_zygotic = rbind.data.frame(countingWindows_zygotic,rpms_allCountingWindows_mt)

countingWindows_zygotic_gr = makeGRangesFromDataFrame(df = countingWindows_zygotic,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)



##### what are the barcodes for which i want this info: 

sampleInfo = read.table("/////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",stringsAsFactors = F)
sampleInfo = sampleInfo %>% filter(grepl(paste(c("Inj_","Unt"), collapse="|"), V3))

barcodesOfInterest = sampleInfo$V2
conditionsOfInterest = sampleInfo$V3

#countingWindows_zygotic_gr = countingWindows_zygotic_gr[1:100] 

for(i in 1:length(barcodesOfInterest)){
  vcfFile_ofInterest = paste0("/////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/snp/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf")
  BamFile_ofInterest = paste0("//////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam")
  vectorizedsample(vcfFile = vcfFile_ofInterest,BamFile = BamFile_ofInterest,countingWindows = countingWindows_zygotic_gr,barcode = barcodesOfInterest[i],condition = conditionsOfInterest[i])
}


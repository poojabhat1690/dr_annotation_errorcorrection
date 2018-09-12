
library(Rsamtools)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)



getMutationsAllTypes = function(vcfFile,BamFile,countingWindows,barcode,condition){
  vcfFile = read.table(vcfFile)
  
  bamFile <- BamFile(file = BamFile)
  open(bamFile)
  params <- ScanBamParam(which  = countingWindows, what = scanBamWhat(),tag=c("RA","MP"))
  
  ### the reads that map to these intervals on both strands are considered, i need to remove this... 
  
  
  aln <- scanBam(bamFile, param = params)
  aln_copy = aln
  numberReads = melt(lapply(aln,function(x) length(x$pos)))
  
  aln[which(numberReads<1)] <- NULL
  vcfFile_mutn = paste(vcfFile$V1,(vcfFile$V2+1),sep="_")
  
  
  #### some counting windows have reads but have 0 mis matches... I want to add an NA there... 
  
  RAtag_num = melt(lapply(aln,function(x) length(x$tag$RA)))
  MPtag_num = melt(lapply(aln,function(x) length(x$tag$MP)))
  
  RA_MP_num = cbind.data.frame(RAtag_num,MPtag_num)
  noMPtags = which(RA_MP_num[,3] == 0)
  
  for(MPtags in 1:length(noMPtags)){
    
    aln[[noMPtags[MPtags]]]$tag$MP = rep(NA, RA_MP_num[noMPtags[MPtags],1])
  }
  cat("this has worked")
  ######## now i can convert the whole file into a dataframe and process it easier...
  aln_chr = as.data.frame(unlist(lapply(aln,function(x) x$rname)))
  aln_RA = as.data.frame(unlist(lapply(aln,function(x) x$tag$RA)))
  aln_MP = as.data.frame(unlist(lapply(aln,function(x) x$tag$MP)))
  aln_position = as.data.frame(unlist(lapply(aln,function(x) x$pos)))
  aln_strand = as.data.frame(unlist(lapply(aln,function(x) x$strand)))
  aln_phred = as.data.frame(unlist(lapply(aln,function(x) as.data.frame(x$qual))))
  alignmentMatrix = cbind.data.frame(aln_chr,aln_RA,aln_MP,aln_position,aln_strand,aln_phred)
  colnames(alignmentMatrix) = c("chr","RA","MP","startPos","strand","qual")
  alignmentMatrix$MP = as.character(alignmentMatrix$MP)
  
  cat ("initial alignment file made")
  alignmentMatrix$splitMP = mapply(function(x) strsplit(x,",",T), alignmentMatrix$MP)
  alignmentMatrix$finalMutationSplit =  mapply(function(x) strsplit(x,":",T), alignmentMatrix$splitMP)
  
  alignmentMatrix$refPosition = (mapply(function(x) lapply(x,function(y) as.numeric(y[3])), alignmentMatrix$finalMutationSplit))
  alignmentMatrix$readPosition = (mapply(function(x) lapply(x,function(y) as.numeric(y[2])), alignmentMatrix$finalMutationSplit))
  alignmentMatrix$mutationType = (mapply(function(x) lapply(x,function(y) as.numeric(y[1])), alignmentMatrix$finalMutationSplit))
  alignmentMatrix = alignmentMatrix %>% mutate (readCorrected = apply(alignmentMatrix,1,function(x)   list(paste0( x$chr,"_",unlist(x$readPosition) + x$startPos)  )) ) %>%  mutate (refCorrected = apply(alignmentMatrix,1,function(x)   list(paste0(x$chr,"_", unlist(x$refPosition) + x$startPos)  ))) 
  
  
  vcfFile$mutn = paste0(vcfFile$V1,"_",vcfFile$V2+1)
  #vcfFile %>% filter(as.character(V1) == as.character(x$chr)) %>% select(mutn)
  
  ###### 
  # alignmentMatrix = alignmentMatrix[1:1000,]
  
  alignmentMatrix =  alignmentMatrix %>% mutate(isSNP = apply(alignmentMatrix,1, function(x) lapply(x$refCorrected, function(y) y %in%  vcfFile_mutn  == T))) 
  alignmentMatrix = alignmentMatrix %>% mutate(snpTYPE = apply(alignmentMatrix,1,function(y) lapply(y$isSNP, function(x) y$mutationType[which(x==T)]) ))
  alignmentMatrix = alignmentMatrix %>%  mutate(whichMutnisSNP = apply(alignmentMatrix,1,function(y) lapply(y$isSNP, function(x) which(x==T)) ))
  alignmentMatrix = alignmentMatrix %>% mutate(bqOfMismatches = apply(alignmentMatrix,1,function(y) lapply(y$readPosition,function(x) utf8ToInt(substr(y$qual,x,x))-33 < 27) )) 
  alignmentMatrix = alignmentMatrix %>% mutate(lowBQTYPE = apply(alignmentMatrix,1,function(y) y$mutationType[which(unlist(y$bqOfMismatches) == T) ] )) 
  alignmentMatrix = alignmentMatrix %>% mutate(whichMutnIsBq = apply(alignmentMatrix,1,function(y) which(unlist(y$bqOfMismatches) == T)  )) 
  
  
  alignmentMatrix = alignmentMatrix %>% mutate(bqidentity =  apply(alignmentMatrix,1,function(y) (paste0(unlist(y$whichMutnIsBq),":",unlist(y$lowBQTYPE )+1)  )  ) ) %>% mutate(SNPidentity =  apply(alignmentMatrix,1,function(y) (paste0(unlist(y$whichMutnisSNP),":",unlist(y$snpTYPE) +1)  )  ) )
  alignmentMatrix = alignmentMatrix %>% mutate(uniqueBad = apply(alignmentMatrix,1,function(y) unique(c(y$bqidentity,y$SNPidentity))))
  
  ### these are the mutations that have to be subtracted... 
  
  alignmentMatrix = alignmentMatrix %>% mutate(mutationTobeRemoved =   apply(alignmentMatrix,1,function(z) list(names(table(unlist(mapply(function(x) lapply(strsplit(x,":",T),function(y) as.numeric(y[2])),z$uniqueBad[z$uniqueBad != ":"]))))))) 
  
  alignmentMatrix = alignmentMatrix %>% mutate( numberTOBeremoved = apply(alignmentMatrix,1,function(z) as.numeric(table(unlist(mapply(function(x) lapply(strsplit(x,":",T),function(y) as.numeric(y[2])),z$uniqueBad[z$uniqueBad != ":"]))))))
  
  mutations_reads = as.data.frame(t(apply(do.call(rbind.data.frame,lapply(apply(alignmentMatrix,1,function(x) strsplit(as.character(x$RA),",",T)),function(y) y[[1]])),1,as.numeric)))
  colnames(mutations_reads) =   nucleotideOrder = c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N")
  
  cat("alignment matrix is completed")
  
  for(i in 1:nrow(alignmentMatrix)){
    rowDf = alignmentMatrix[i,]
    if(length(rowDf$mutationTobeRemoved[[1]][[1]])>0){
      mutations_reads[i, ][lapply(rowDf$mutationTobeRemoved,function(x) as.numeric(unlist(x)) )[[1]]] <-  mutations_reads[i, ][lapply(rowDf$mutationTobeRemoved,function(x) as.numeric(unlist(x)) )[[1]]] -   lapply(rowDf$numberTOBeremoved,function(x) unlist(x)) [[1]]
      
    }
  }
  
  alignmentMatrix = cbind.data.frame(mutations_reads,alignmentMatrix)
  
  save( alignmentMatrix,file = paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/mutations_data_",barcode,condition,".Rdata"))
  
  
}



rpms_allCountingWindows = read.table("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)
countingWindows_zygotic = rpms_allCountingWindows %>% mutate(mean_tp9 = rowMeans(data.frame(Inj_R1_TP9,Inj_R2_TP9,Inj_R3_TP9))) %>% filter(mean_tp9 >= 4) %>% select(V1:V6)
countingWindows_zygotic_gr = makeGRangesFromDataFrame(df = countingWindows_zygotic,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)


##### what are the barcodes for which i want this info: 

sampleInfo = read.table("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",stringsAsFactors = F)
sampleInfo = sampleInfo %>% filter(grepl(paste(c("Inj_"), collapse="|"), V3))

barcodesOfInterest = sampleInfo$V2
conditionsOfInterest = sampleInfo$V3

#countingWindows_zygotic_gr = countingWindows_zygotic_gr[1:100] 

for(i in 1:length(barcodesOfInterest)){
  vcfFile_ofInterest = paste0("///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/snp/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf")
  BamFile_ofInterest = paste0("////groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/filter/combinedFile_",barcodesOfInterest[i],".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam")
  getMutationsAllTypes(vcfFile = vcfFile_ofInterest,BamFile = BamFile_ofInterest,countingWindows = countingWindows_zygotic_gr,barcode = barcodesOfInterest[i],condition = conditionsOfInterest[i])
}







#   vcfFile = read.table(vcfFile)
#   
#   bamFile <- BamFile(file = bamFile_current)
#   open(bamFile)
#   params <- ScanBamParam(which  = countingWindows, what = scanBamWhat(),tag=c("RA","MP"))
#   
#   ### the reads that map to these intervals on both strands are considered, i need to remove this... 
#   
#   
#   aln <- scanBam(bamFile, param = params)
#   aln_copy = aln
#   numberReads = melt(lapply(aln,function(x) length(x$pos)))
#   
#   aln[which(numberReads<1)] <- NULL
#   vcfFile_mutn = paste(vcfFile$V1,(vcfFile$V2+1),sep="_")
#   
#   
#   
#   
#   
#   
#   #### some counting windows have reads but have 0 mis matches... I want to add an NA there... 
#   
#   RAtag_num = melt(lapply(aln,function(x) length(x$tag$RA)))
#   MPtag_num = melt(lapply(aln,function(x) length(x$tag$MP)))
#   
#   RA_MP_num = cbind.data.frame(RAtag_num,MPtag_num)
#   noMPtags = which(RA_MP_num[,3] == 0)
#   
#   for(MPtags in 1:length(noMPtags)){
#   
#     aln[[noMPtags[MPtags]]]$tag$MP = rep(NA, RA_MP_num[noMPtags[MPtags],1])
#   }
#   
#   ######## now i can convert the whole file into a dataframe and process it easier...
#   aln_chr = as.data.frame(unlist(lapply(aln,function(x) x$rname)))
#   aln_RA = as.data.frame(unlist(lapply(aln,function(x) x$tag$RA)))
#   aln_MP = as.data.frame(unlist(lapply(aln,function(x) x$tag$MP)))
#   aln_position = as.data.frame(unlist(lapply(aln,function(x) x$pos)))
#   aln_strand = as.data.frame(unlist(lapply(aln,function(x) x$strand)))
#   aln_phred = as.data.frame(unlist(lapply(aln,function(x) as.data.frame(x$qual))))
#   alignmentMatrix = cbind.data.frame(aln_chr,aln_RA,aln_MP,aln_position,aln_strand,aln_phred)
#   colnames(alignmentMatrix) = c("chr","RA","MP","startPos","strand","qual")
#   alignmentMatrix$MP = as.character(alignmentMatrix$MP)
#   
#   alignmentMatrix$splitMP = mapply(function(x) strsplit(x,",",T), alignmentMatrix$MP)
#   alignmentMatrix$finalMutationSplit =  mapply(function(x) strsplit(x,":",T), alignmentMatrix$splitMP)
#   
#   alignmentMatrix$refPosition = (mapply(function(x) lapply(x,function(y) as.numeric(y[3])), alignmentMatrix$finalMutationSplit))
#   alignmentMatrix$readPosition = (mapply(function(x) lapply(x,function(y) as.numeric(y[2])), alignmentMatrix$finalMutationSplit))
#   alignmentMatrix$mutationType = (mapply(function(x) lapply(x,function(y) as.numeric(y[1])), alignmentMatrix$finalMutationSplit))
#   alignmentMatrix = alignmentMatrix %>% mutate (readCorrected = apply(alignmentMatrix,1,function(x)   list(paste0( x$chr,"_",unlist(x$readPosition) + x$startPos)  )) ) %>%  mutate (refCorrected = apply(alignmentMatrix,1,function(x)   list(paste0(x$chr,"_", unlist(x$refPosition) + x$startPos)  ))) 
#   
#   
#   vcfFile$mutn = paste0(vcfFile$V1,"_",vcfFile$V2+1)
#   vcfFile %>% filter(as.character(V1) == as.character(x$chr)) %>% select(mutn)
#   
# ###### 
#  # alignmentMatrix = alignmentMatrix[1:1000,]
#    
#   alignmentMatrix =  alignmentMatrix %>% mutate(isSNP = apply(alignmentMatrix,1, function(x) lapply(x$refCorrected, function(y) y %in%  vcfFile_mutn  == T))) 
#   alignmentMatrix = alignmentMatrix %>% mutate(snpTYPE = apply(alignmentMatrix,1,function(y) lapply(y$isSNP, function(x) y$mutationType[which(x==T)]) ))
#   alignmentMatrix = alignmentMatrix %>%  mutate(whichMutnisSNP = apply(alignmentMatrix,1,function(y) lapply(y$isSNP, function(x) which(x==T)) ))
#   alignmentMatrix = alignmentMatrix %>% mutate(bqOfMismatches = apply(alignmentMatrix,1,function(y) lapply(y$readPosition,function(x) utf8ToInt(substr(y$qual,x,x))-33 < 27) )) 
#   alignmentMatrix = alignmentMatrix %>% mutate(lowBQTYPE = apply(alignmentMatrix,1,function(y) y$mutationType[which(unlist(y$bqOfMismatches) == T) ] )) 
#   alignmentMatrix = alignmentMatrix %>% mutate(whichMutnIsBq = apply(alignmentMatrix,1,function(y) which(unlist(y$bqOfMismatches) == T)  )) 
#   
#   
#   alignmentMatrix = alignmentMatrix %>% mutate(bqidentity =  apply(alignmentMatrix,1,function(y) (paste0(unlist(y$whichMutnIsBq),":",unlist(y$lowBQTYPE )+1)  )  ) ) %>% mutate(SNPidentity =  apply(alignmentMatrix,1,function(y) (paste0(unlist(y$whichMutnisSNP),":",unlist(y$snpTYPE) +1)  )  ) )
#   alignmentMatrix = alignmentMatrix %>% mutate(uniqueBad = apply(alignmentMatrix,1,function(y) unique(c(y$bqidentity,y$SNPidentity))))
#   
#   ### these are the mutations that have to be subtracted... 
#   
#   alignmentMatrix = alignmentMatrix %>% mutate(mutationTobeRemoved =   apply(alignmentMatrix,1,function(z) list(names(table(unlist(mapply(function(x) lapply(strsplit(x,":",T),function(y) as.numeric(y[2])),z$uniqueBad[z$uniqueBad != ":"]))))))) 
#   
#   alignmentMatrix = alignmentMatrix %>% mutate( numberTOBeremoved = apply(alignmentMatrix,1,function(z) as.numeric(table(unlist(mapply(function(x) lapply(strsplit(x,":",T),function(y) as.numeric(y[2])),z$uniqueBad[z$uniqueBad != ":"]))))))
# 
#   mutations_reads = as.data.frame(t(apply(do.call(rbind.data.frame,lapply(apply(alignmentMatrix,1,function(x) strsplit(as.character(x$RA),",",T)),function(y) y[[1]])),1,as.numeric)))
#   colnames(mutations_reads) =   nucleotideOrder = c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N")
# 
#   cat("alignment matrix is completed")
#   
#   for(i in 1:nrow(alignmentMatrix)){
#     rowDf = alignmentMatrix[i,]
#     if(length(rowDf$mutationTobeRemoved[[1]][[1]])>0){
#         mutations_reads[i, ][lapply(rowDf$mutationTobeRemoved,function(x) as.numeric(unlist(x)) )[[1]]] <-  mutations_reads[i, ][lapply(rowDf$mutationTobeRemoved,function(x) as.numeric(unlist(x)) )[[1]]] -   lapply(rowDf$numberTOBeremoved,function(x) unlist(x)) [[1]]
# 
#     }
#   }
#   
#   alignmentMatrix = cbind.data.frame(mutations_reads,alignmentMatrix)
# 
#   save( alignmentMatrix,file = "///groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/mutations_data.Rdata")
#   
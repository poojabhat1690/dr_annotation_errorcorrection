##### common functions for the detection of SNPs in transcripts and for determining the probability to find reads with differnet number of mutations 


### theme_ameres
## additional theme for the ggplot package
###

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


##### removing SNPs and getting high quality bases from bam files... 

  #### input - vcf file and bam file path, counting window ranges from which the bam file is subset.
  #### output - the cleaned up alignment file  - for SNPs and base quality 
  ###### hard base quality of 27 is set

removeSNPs_bq = function(vcfFile,bamFile_current,countingWindows){
  vcfFile = read.table(vcfFile)
  
  bamFile <- BamFile(file = bamFile_current)
  open(bamFile)
  params <- ScanBamParam(which  = countingWindows, what = scanBamWhat(),tag=c("RA","MP"))
  
  ### the reads that map to these intervals on both strands are considered, i need to remove this... 
  
  
  aln <- scanBam(bamFile, param = params)
  aln_copy = aln
  numberReads = melt(lapply(aln,function(x) length(x$pos)))
  aln[which(numberReads<1)] <- NULL
  
  vcfFile_mutn = paste(vcfFile$V1,(vcfFile$V2+1),sep="_")
 
  for(m in 1:length(aln)){
    
    if(length(aln[[m]]$MP)>=1){ ####### only want to go into this if there are some mutations.. so checking the MP tag 
    
      for(i in 1:length(aln[[m]]$rname)){
      
              differentMutations = strsplit(aln[[m]]$tag$MP[i],",",T)
            if(is.na(differentMutations) == F){
              for(j in 1:length(differentMutations[[1]])){ ##### for 1 read
                SingleMutationPosition = strsplit(differentMutations[[1]][[j]][1],":",T)
                mutPosition = aln[[m]]$pos[i] + (as.numeric(SingleMutationPosition[[1]][3])) ### this neeeds to be the position of the reference as the vcf will report the position wrt reference
                chr_mutposition = paste(aln[[m]]$rname[1],mutPosition,sep="_")
                isMutn = chr_mutposition %in% vcfFile_mutn
                if(isMutn == TRUE){
                  ### getting the RA tag for this read and checking the mutation in that 
                  mutNum = as.numeric(SingleMutationPosition[[1]][1])
                  RATag_read = aln[[m]]$tag$RA[i]
                  split_RAtag =  strsplit(RATag_read,",",T)[[1]]
                  split_RAtag[mutNum+1]<- as.numeric(split_RAtag[mutNum+1]) -1  ### +1 as the mutation numbers have been named from 0 and the RA tag vector starts with 1
                  RAtag_collapsed =  paste(as.character(split_RAtag),collapse=",",sep="")
                  aln[[m]]$tag$RA[i] <- RAtag_collapsed
                }
                else{
                  #   ### if the base is not a mutation, i want to check the quality of the base...
                  baseQualOf_j = utf8ToInt(as.character(aln[[m]]$qual[i][[1]][ (as.numeric(SingleMutationPosition[[1]][2]))]))-33 ### the base quality has to be of the read 
                  
                  if(baseQualOf_j < 27){
                    mutNum = as.numeric(SingleMutationPosition[[1]][1])
                    RATag_read = aln[[m]]$tag$RA[i]
                    split_RAtag =  strsplit(RATag_read,",",T)[[1]]
                    split_RAtag[mutNum+1]<- as.numeric(split_RAtag[mutNum+1]) -1   ### +1 as the mutation numbers have been named from 0 and the RA tag vector starts with 1
                    RAtag_collapsed =  paste(as.character(split_RAtag),collapse=",",sep="")
                    aln[[m]]$tag$RA[i] <- RAtag_collapsed
                    #cat(baseQualOf_j)
                    
                  } 
                  
                } 
                
              }
            }
            
          }
    }
    cat(m)
  }
  return(aln)
}


###########


### this function gets the fraction of reads with 1,2,3,4 mutations per gene/region that has been speicifed in the counting windows...

### input the aligned file that has been corrected for SNPs and basequality

getParamsFromBam = function(aln){
  RAtags = lapply(aln,function(x) mapply(function(y) strsplit(y,",",T),x$tag$RA))
  MPtags = lapply(aln,function(x) mapply(function(y) strsplit(y,",",T),x$tag$MP))
  phredScores = lapply(aln,function(x) as.data.frame(x$qual))
  split_singleNucleotides = lapply(phredScores,function(y) mapply(function(x) strsplit(x,"",T) , y$x))
  RAtags = lapply(RAtags,function(x) lapply(x,function(y) as.numeric(y) ))
  
  
  RAtags=lapply(RAtags,function(x) do.call(rbind.data.frame,x))
  nucleotideOrder = c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N")
  lengthOfColumns  = melt(lapply(RAtags,ncol))
  RAtags[which(lengthOfColumns$value != 25)] = NULL
  
  RAtags = lapply(RAtags, setNames, nm = nucleotideOrder)
  
  #### i want to integrate base quality information at this level
  
  RAtags_copy =RAtags
  RAtags = lapply(RAtags, setNames, nm = nucleotideOrder)
  
  lengthOfColumns  = melt(lapply(RAtags,nrow)) ### calculating the number of reads 
  
  RAtags = RAtags[which(lengthOfColumns$value > 10)]
  
  RAtags_readWith1mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 1)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith2mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 2)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith3mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 3)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith4mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 4)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith4mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 5)))) )#### per counting window how many reads have 1 mutation
  colnames(RAtags_readWith1mutation) = nucleotideOrder
  colnames(RAtags_readWith2mutation) = nucleotideOrder
  colnames(RAtags_readWith3mutation) = nucleotideOrder
  colnames(RAtags_readWith4mutation) = nucleotideOrder
  
  
  numberOfreads = melt(lapply(aln,function(x) length(x$seq)))
  strand_samples = melt(lapply(aln,function(x) x$strand[1]))
  numberOfreads$strand = strand_samples$value
  
  #### normalizing to the number of reads ... i.e getting the fraction of total reads that have different mutations
  library(dplyr)
  calculateFraction = function(df,numReads){
    df = as.data.frame(df)
    
    df$L1 = row.names(df)
    df = join(df,numReads)
    df_fraction = df[,c(1:25)]/df$value
    df_fraction$strand = df$strand
    df_fraction$numRead = df$value
    minusTab = df_fraction %>% dplyr::filter(strand == "-") %>%
      dplyr::select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)
    plusTab = df_fraction %>% dplyr::filter(strand == "+") %>%
      dplyr::select(A_A , G_G , C_C , T_T, A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G )
    
    df_fraction = rbind(minusTab,plusTab)
    df_fraction = df_fraction %>% dplyr::select(A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G) 
    return(df_fraction)
  }
  
  library(plyr)
  library(ggplot2)
  
  mutation_n1 = calculateFraction(df =RAtags_readWith1mutation,numReads = numberOfreads )
  mutation_n2 = calculateFraction(df =RAtags_readWith2mutation,numReads = numberOfreads )
  mutation_n3 = calculateFraction(df =RAtags_readWith3mutation,numReads = numberOfreads )
  mutation_n4 = calculateFraction(df =RAtags_readWith4mutation,numReads = numberOfreads )
  
  numberOfreadsWithMutations = list(mutation_n1,mutation_n2,mutation_n3,mutation_n4)
  names(numberOfreadsWithMutations) = c("mut_1","mut_2","mut_3","mut_4")
  return(numberOfreadsWithMutations)
}







##### getting the absolute number of mutations per gene or counting windows that has been specified... 

### input : the alignment file cleaned up foe SNPs and low quality bases... 

getNumberOfMutationsPerGene = function(aln){ ### does not caluclate the fraction of reads with mutations, but reports only the number of mutations...
  RAtags = lapply(aln,function(x) mapply(function(y) strsplit(y,",",T),x$tag$RA))
  MPtags = lapply(aln,function(x) mapply(function(y) strsplit(y,",",T),x$tag$MP))
  phredScores = lapply(aln,function(x) as.data.frame(x$qual))
  split_singleNucleotides = lapply(phredScores,function(y) mapply(function(x) strsplit(x,"",T) , y$x))
  #split_singleNucleotides = lapply(split_singleNucleotides,function(x) lapply(x,function(y) mapply(function(z) utf8ToInt(z)-33,y)))
  RAtags = lapply(RAtags,function(x) lapply(x,function(y) as.numeric(y) ))
  
  #RAtags = do.call(rbind.data.frame,mapply(function(x) strsplit(x,",",T),aln$`chr13:50307021-50307151`$tag$RA))
  
  
  RAtags=lapply(RAtags,function(x) do.call(rbind.data.frame,x))
  nucleotideOrder = c("A_A" ,"A_C" ,"A_G", "A_T", "A_N", "C_A",  "C_C" ,"C_G", "C_T", "C_N", "G_A" ,"G_C",  "G_G", "G_T" ,"G_N" ,"T_A", "T_C", "T_G",  "T_T", "T_N", "N_A" ,"N_C" ,"N_G" ,"N_T" ,"N_N")
  lengthOfColumns  = melt(lapply(RAtags,ncol))
  RAtags[which(lengthOfColumns$value != 25)] = NULL
  
  RAtags = lapply(RAtags, setNames, nm = nucleotideOrder)
  
  #### i want to integrate base quality information at this level
  
  RAtags_copy =RAtags
  
  RAtags = lapply(RAtags, setNames, nm = nucleotideOrder)
  
  lengthOfColumns  = melt(lapply(RAtags,nrow)) ### calculating the number of reads 
  
  RAtags = RAtags[which(lengthOfColumns$value > 10)]
  
  RAtags_readWith1mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 1)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith2mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 2)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith3mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 3)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith4mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 4)))) )#### per counting window how many reads have 1 mutation
  RAtags_readWith4mutation = do.call(rbind,lapply(RAtags,function(x) apply(x,2,function(x) length(which(x == 5)))) )#### per counting window how many reads have 1 mutation
  colnames(RAtags_readWith1mutation) = nucleotideOrder
  colnames(RAtags_readWith2mutation) = nucleotideOrder
  colnames(RAtags_readWith3mutation) = nucleotideOrder
  colnames(RAtags_readWith4mutation) = nucleotideOrder
  
  
  numberOfreads = melt(lapply(aln,function(x) length(x$seq)))
  strand_samples = melt(lapply(aln,function(x) x$strand[1]))
  numberOfreads$strand = strand_samples$value
  
  #### normalizing to the number of reads ... i.e getting the fraction of total reads that have different mutations
  library(dplyr)
  calculateFraction = function(df,numReads){
    df = as.data.frame(df)
    
    df$L1 = row.names(df)
    df = join(df,numReads)
    df_fraction = df[,c(1:25)]
    df_fraction$strand = df$strand
    df_fraction$numRead = df$value
    df_fraction$range = df$L1
    minusTab = df_fraction %>% dplyr::filter(strand == "-") %>%
      dplyr::select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,numRead,range,strand)
    plusTab = df_fraction %>% dplyr::filter(strand == "+") %>%
      dplyr::select(A_A , G_G , C_C , T_T, A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,numRead,range,strand)
    
    df_fraction = rbind(minusTab,plusTab)
    df_fraction = df_fraction %>% dplyr::select(A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G,range,strand,numRead) 
    return(df_fraction)
  }
  
  library(plyr)
  library(ggplot2)
  
  mutation_n1 = calculateFraction(df =RAtags_readWith1mutation,numReads = numberOfreads )
  mutation_n2 = calculateFraction(df =RAtags_readWith2mutation,numReads = numberOfreads )
  mutation_n3 = calculateFraction(df =RAtags_readWith3mutation,numReads = numberOfreads )
  mutation_n4 = calculateFraction(df =RAtags_readWith4mutation,numReads = numberOfreads )
  
  numberOfreadsWithMutations = list(mutation_n1,mutation_n2,mutation_n3,mutation_n4)
  names(numberOfreadsWithMutations) = c("mut_1","mut_2","mut_3","mut_4")
  return(numberOfreadsWithMutations)
}



### caculate the means of fracrion of conversions across genes and plot mean fraction of conversions (1-4) for the counting winfoes of interest.
### input : output from function getParamsFromBam

calcMeans = function(df  ){
  means_mutationProbability = do.call(rbind,lapply(df,function(x) colMeans(x)))
  means_mutationProbability = as.data.frame(means_mutationProbability)
  means_mutationProbability = melt(means_mutationProbability)
  means_mutationProbability$nmut = c(1:4)
  colSacles = c(RColorBrewer::brewer.pal(n = 8,name = "Dark2"), RColorBrewer::brewer.pal(n=4,name="Set1"))
  p = ggpubr::ggline(data = means_mutationProbability,y='value',x='nmut',group='variable',col='variable',size = 1) + scale_color_manual(values = colSacles) + theme_ameres(type="barplot")
  return(p)
}









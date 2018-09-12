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


rpms_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/allCountingWindows_TCConversionRate.txt",sep="\t",stringsAsFactors = F ,header = T)
RPMS_allCountingWindows = read.table("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/dataTables/RPMData_allCountingWindows.txt",sep="\t",stringsAsFactors = F ,header = T)

mir430 = rpms_allCountingWindows[c(1:3),]

mir430_unt =  mir430 %>% select(contains("Unt"))
mir430_injR1 =  mir430 %>% select(contains("Inj_R1")) - mir430_unt
mir430_injR2 =  mir430 %>% select(contains("Inj_R2")) - mir430_unt
mir430_injR3 =  mir430 %>% select(contains("Inj_R3")) - mir430_unt

allMir430 = list(mir430_injR1,mir430_injR2,mir430_injR3)
names(allMir430) = c("R1","R2","R3")
allMir430 = melt(allMir430)
allMir430$locus = paste0("locus_",c(1:3))
allMir430$time = unlist(lapply(strsplit(as.character(allMir430$variable),"_",T),function(x) x[[3]]))
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/probabilityOfTC/mir430_conversionRates.pdf")
ggpubr::ggline(data = allMir430,x = 'time',y='value',group='locus',col='locus',size=1.5) + facet_grid(L1 ~.) + scale_colour_brewer(palette = "Dark2") + theme_ameres(type = "barplot") + ylab("conversion rate")
dev.off()
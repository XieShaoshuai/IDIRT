#load mix protein data

load_mix <- function(path, col_select,i){
  #load data
  mix <- read.delim(path)
  metadata <- read.delim(col_select,sep=",")
  
  #selct master protein groups
  mix <- mix %>% filter(Master=="IsMasterProtein", Protein.FDR.Confidence.Combined=="High", 
                        Number.of.Unique.Peptides>1,Contaminant=="False")
  
  #select columns
  col <- metadata %>% filter(experiment==i)
  mix <- mix[,c('Accession','Description',col$colname)]
  
  #rename column 
  colnames(mix) <- c('Accession','Description',col$new_colname)
  
  #replace NA to 0
  mix<- replace(mix, is.na(mix), 0)
  
  #calculate ratio using normalized abundance
  mix$protein_ratio1 <- mix$Sample_1/(mix$Sample_1+mix$Control_1)
  mix$protein_ratio2 <- mix$Sample_2/(mix$Sample_2+mix$Control_2)
  mix$protein_ratio3 <- mix$Sample_3/(mix$Sample_3+mix$Control_3)
  mix$protein_ratio4 <- mix$Sample_4/(mix$Sample_4+mix$Control_4)
  mix$protein_ratio5 <- mix$Sample_5/(mix$Sample_5+mix$Control_5)
  mix$protein_ratio6 <- mix$Sample_6/(mix$Sample_6+mix$Control_6)
  
  #remove NA
  mix <- mix %>% filter(mix$protein_ratio1!="NaN")
  
  #calculate protein median ratio
  ratio <- mix %>% dplyr::select(starts_with("protein_ratio"))
  ratio <- as.matrix(ratio)
  mix$protein_median_ratio <- rowMedians(ratio, na.rm=TRUE)
  
  
  #plot correlation
  d <- mix %>% dplyr::select(contains("protein_"))
  p <-ggpairs(d, title="Ratio correlation") 
  
  ggsave(p,file=paste0("./Output/Mix/ratio_correlation_",i,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  
  #selct colname and retrun
  mix <-  mix %>% dplyr::select(Accession,Description,starts_with("protein_"))
  return(mix)
}

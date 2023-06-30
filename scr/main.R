library(mixR)
library(tidyverse)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(patchwork)
library(matrixStats)
library(GGally)
library(VennDiagram)
library(ggrepel)
library(BayesFactor)
library(matrixStats)


set.seed(123456)
#parameter setting
mix_protein_path <- "./data/Mix data/Mix_all_Proteins.txt"
idirt_protein_path <- "./data/IDIRT data/IDIRTall_Proteins.txt"

col_select <- "./data/Mix data/metadata.csv"
col_select_idirt <- "./data/IDIRT data/IDIRT_metadata swap.csv"

buffer <- c(2,3,16,19,25,29)

#-----------------------------------------------
#MIX data analysis
#-----------------------------------------------

#load source
source("./scr/Mix/1_load_data.R")
source("./scr/Mix/5_check_overlap.R")
source("./scr/Mix/6_get_cutoff.R")
source("./scr/idirt/1_load_idirt.R")
source("./scr/idirt/bayesian.R")
source("./scr/idirt/bayesian_global.R")


#load MIX proteome protein and peptide data
mix <- load_mix(mix_protein_path,col_select)

#check the overlap between Mix and IDIRT data
check_ovelap(mix,idirt_protein_path)

#establish cut-off, method1: just use the distribution of mix protein ratio as a cutoff
ratio <- mix %>% dplyr::filter(protein_median_ratio > 0 & protein_median_ratio<1)
ratio <- ratio$protein_median_ratio
#stimulate peak distribtion
x_binned <- mixR::bin(ratio, brks = seq(0, 1, length = 20))
mod <- mixR::mixfit(x_binned, ncomp = 1,init.method = "kmeans")
cutoff <- mod[["mu"]] + 2*mod[["sd"]]
cutoff_list <- list("cutoff"=cutoff,"sd"=mod[["sd"]], mu = mod[["mu"]])
cutoff_list
plot(mod)



#--------------------------------------------------------------
#IDIRT data analysis
#--------------------------------------------------------------


#detach("package:mixR", unload = TRUE)

#library(dplyr)
data_all <- as.data.frame(as.matrix("NA",1,0))
colnames(data_all) <- "Accession"
for(k in buffer){
  #load IDIRT data
  idirt  <- load_idirt(idirt_protein_path, col_select_idirt, k, mix)
  
  idirt_raw <- idirt
  idirt <- idirt[rowSums((idirt %>% select(contains("idirt_")))>0)>=3,] #proteins found >3 replicates for Bayesian test
  idirt <- idirt[rowSums(!is.na(idirt %>% select(contains("idirt_"))))>=3,]
  
  idirt_bayesian <- bayesian2(idirt,mix,cutoff_list,k)
  
  idirt_bayesian <- idirt_bayesian %>% select(Accession,bf)
  
  data <- left_join(idirt_raw,idirt_bayesian,by="Accession")
  write.table(data,file=paste0("Buffer ",k,".csv"),sep=",")
  
  data_temp <- data %>% select(Accession,median_ratio,bf)
  colnames(data_temp) <- c("Accession",paste0("Buffer",k,"_Idirt ratio"),paste0("Buffer",k,"_Bayes Factor"))
  
  data_all <- full_join(data_temp,data_all, by="Accession")
  
}

write.table(data_all,"all_data.csv",sep=",")

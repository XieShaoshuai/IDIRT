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
mix_protein_path <- "./data/Mix data/Mix_90min_Proteins.txt"
idirt_protein_path <- "./data/IDIRT data/IDIRT_Proteins.txt"

col_select <- "./data/Mix data/metadata.csv"
col_select_idirt <- "./data/IDIRT data/IDIRT_metadata.csv"

buffer <- c(2,3,16,19,25,29)

#-----------------------------------------------
#MIX data analysis
#-----------------------------------------------

#load source
source("./scr/Mix/1_load_data.R")
source("./scr/Mix/5_check_overlap.R")
source("./scr/Mix/6_get_cutoff.R")

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
source("./scr/idirt/1_load_idirt.R")
source("./scr/idirt/bayesian.R")
source("./scr/idirt/bayesian_global.R")

#detach("package:mixR", unload = TRUE)

#library(dplyr)
data <- as.data.frame(as.matrix("NA",1,0))
colnames(data) <- "Accession"
for(k in buffer){
  #load IDIRT data
  data1  <- load_idirt(idirt_protein_path, col_select_idirt, k, mix)
  
  data3 <- bayesian2(data1,mix,cutoff_list,k)
  
  data3_temp <- data3 %>% select(Accession,bf)
  data1 <- left_join(data1,data3_temp,by="Accession")
  write.table(data1,file=paste0("Buffer ",k,".csv"),sep=",")
  
  data_temp <- data3 %>% select(Accession,Description,bf,median_ratio)
  colnames(data_temp) <- c("Accession",paste0("Buffer",k,"_Bayes Factor"),paste0("Buffer",k,"_Idirt ratio"))
  
  data <- full_join(data_temp,data, by="Accession")
  
}

dd <- mix %>% select(Accession,Description,protein_median_ratio)

data <- left_join(data,dd,by="Accession")
write.table(data,"all_data.csv",sep=",")

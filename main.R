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
library(ggrepel)


set.seed(123456)
#parameter setting
mix_protein_path <- "./data/Mix data/Mix_all_Proteins.txt"
idirt_protein_path <- "./data/IDIRT data/IDIRTall_Proteins.txt"

col_select <- "./data/Mix data/metadata.csv"
col_select_idirt <- "./data/IDIRT data/IDIRT_metadata swap.csv"

buffer <- c(2,3,16,19,25,29)
experiment <- c(1,2)  # 1: IDIRT Heavy for BG1 CDH1-HA, #2 IDIRT swap, Heavy for BG1 wt

#-----------------------------------------------
#MIX data analysis
#-----------------------------------------------

#load source
source("./scr/Mix/1_load_data.R")
source("./scr/Mix/6_get_cutoff.R")
source("./scr/idirt/1_load_idirt.R")
source("./scr/idirt/bayesian_global.R")
source("./scr/idirt/stochiometry.R")
source("./scr/idirt/get_interactors.R")


#load MIX proteome protein and peptide data
mix1 <- load_mix(mix_protein_path,col_select,1)
mix2 <- load_mix(mix_protein_path,col_select,2)


#establish cut-off, method1: just use the distribution of mix protein ratio as a cutoff
cutoff_1 <- get_cutoff(mix1,1)
cutoff_2 <- get_cutoff(mix2,2)


#--------------------------------------------------------------
#IDIRT data analysis
#--------------------------------------------------------------
#get interactors in each condition
con2 <- get_interactor(2)
con3 <- get_interactor(3)
con16 <- get_interactor(16)
con19 <- get_interactor(19)
con25 <- get_interactor(25) 
con29 <- get_interactor(29)


#combine all interactors
list2 <- con2 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list3 <- con3 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list16 <- con16 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list19 <- con19 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list25 <- con25 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list29 <- con29 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
interactors <- unique(c(list2$Accession,list3$Accession,list16$Accession,list19$Accession,list25$Accession,list29$Accession))

write.table(interactors,"interactors.txt",row.names = FALSE)


list2 <- con2 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
list3 <- con3 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
list16 <- con16 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
list19 <- con19 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
list25 <- con25 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
list29 <- con29 %>% filter(group=="Stable interactor") %>% 
  select(Gene,group)
interactors <- unique(c(list2$Gene,list3$Gene,list16$Gene,list19$Gene,list25$Gene,list29$Gene))

write.table(interactors,"stable interactors.txt",row.names = FALSE)

list2 <- con2 %>% filter(group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list3 <- con3 %>% filter(group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list16 <- con16 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list19 <- con19 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list25 <- con25 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
list29 <- con29 %>% filter(group=="Stable interactor"|group=="Unstable interactor"|group=="Bait") %>% 
  select(Accession,Gene,group)
interactors <- unique(c(list2$Gene,list3$Gene,list16$Gene,list19$Gene,list25$Gene,list29$Gene))

write.table(interactors,"interactors2.txt",row.names = FALSE)


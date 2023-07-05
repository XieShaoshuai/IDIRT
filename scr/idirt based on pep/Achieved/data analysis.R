library(ggplot2)
library(dplyr)
library(bayestestR)
library(rstanarm)

set.seed(123456)

# peptides level comapre
source("./scr/idirt based on pep/pep_ratio.R")



# read file
psm <- read.delim("./data/IDIRT data/IDIRTall_PSM_PSMs.txt")

#change Heavy and light to Sample,Control----
psm <- psm %>%
  mutate(Quan.Channel = if_else((grepl("Heavy", Quan.Channel) & grepl("IDIRT_buffer", Spectrum.File))|
                                  (grepl("Light", Quan.Channel) & grepl("IDIRT_SWAP", Spectrum.File)), "Sample", "Control"))

#remove Contaminant and reomve blank in protein accession
psm <- psm %>% filter(Master.Protein.Accessions!="", Contaminant=="False") 


#idirt data analysis

get_ratio <- function(psm,buffer)  #buffer: buffer2_
{
  
}


#IDIRT and IDIRT swap analysis seperately
get_ratio_0 <- function(idirt,buffer){
  filname <- paste0(idirt,"_",buffer)  #get the 4 replicates in buffer
  df <- psm %>% filter(grepl(filname,Spectrum.File))
  protein_accession <- unique(df$Master.Protein.Accessions) #get unique protein accession
  data_container <-  data.frame(accession=protein_accession) #set a container to save the ratio for each protein_accession
  for (acc in protein_accession){
    df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
    df_temp <- replace(df_temp,is.na(df_temp),0)  #replace NA to 0
    df_temp <- df_temp %>%                        #group same sequence by file name and quan channel          
      group_by(Sequence,Quan.Channel,Spectrum.File) %>%
      summarise(C_sum = sum(Precursor.Abundance, na.rm = TRUE))
    #filter protein with unique sequence(peptide)>2 in total 4 replicates
    if(length(unique(df_temp$Sequence))>2){
      #analysis separately by replicate
      replicate <- unique(df_temp$Spectrum.File)
      
      ratio <- ""
      for(rep in replicate){
        df_temp2 <- df_temp %>%  filter(grepl(rep,Spectrum.File))
        peptide <- data.frame(Sequence=unique(df_temp2$Sequence))
        
        #extract sample channel
        peptide_S <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Sample"),by="Sequence")  #add heavy channel
        peptide_S$Quan.Channel <- "Sample"
        peptide_S$C_sum[is.na(peptide_S$C_sum)] <-0
        
        #Extract control channel
        peptide_C <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Control"),by="Sequence")  #add light channel
        peptide_C$Quan.Channel <- "Control"
        peptide_C$C_sum[is.na(peptide_C$C_sum)] <-0
        peptide <- rbind(peptide_S,peptide_C)
        
        peptide <- peptide %>% group_by(Sequence)%>%
          summarise(ratio = C_sum[Quan.Channel == "Sample"] /(C_sum[Quan.Channel == "Control"]+C_sum[Quan.Channel == "Sample"]),
                    group ="F1")
        
        #remove NaN
        peptide <- peptide[peptide$ratio!="NaN",]
        
        #combine all ratio
        ratio <- c(ratio,peptide$ratio) 
      }
      
      
    }
    #Generate posterior data from prior data(ratio) 
    alpha_prior_A <- 2
    beta_prior_A <- 2
    
    # Observed data
    # Replace these with your actual data
    data_A <- as.numeric(ratio[-1])
    data_A <- na.omit(data_A)
    
    # Number of trials and number of successes
    n_trials_A <- length(data_A)
    successes_A <- sum(data_A)
    
    # Posterior distribution parameters
    alpha_posterior_A <- alpha_prior_A + successes_A
    beta_posterior_A <- beta_prior_A + n_trials_A - successes_A
    
    # Generate samples from the posterior distributions
    samples_A <- rbeta(1000, shape1 = alpha_posterior_A, shape2 = beta_posterior_A)
    
    #get mu and sd
    mu <- mean(samples_A)
    sd <- sd(samples_A)
    
    data_container$mu[data_container$accession==acc]<-mu
    data_container$sd[data_container$accession==acc]<-sd
  }
  return(data_container)
}
#1 extract buffer
df <- psm %>% filter(grepl("IDIRT_buffer2_",Spectrum.File))




protein_accession <- unique(df$Master.Protein.Accessions)
data_container <- data.frame(accession=unique(df$Master.Protein.Accessions))


for (acc in protein_accession){
  #extract protein
  df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
  
  #replace NA to 0
  df_temp <- replace(df_temp,is.na(df_temp),0)
  #sum intensity of the same peptide
  df_temp <- df_temp %>% 
    group_by(Sequence,Quan.Channel,Spectrum.File) %>%
    summarise(C_sum = sum(Precursor.Abundance, na.rm = TRUE))
  
  #if unique peptide number is less than 2 in all replicate, skip
  if(length(unique(df_temp$Sequence))>2){
    #analysis separately by replicate
    replicate <- unique(df_temp$Spectrum.File)
    
    ratio <- ""
    for(rep in replicate){
      df_temp2 <- df_temp %>%  filter(grepl(rep,Spectrum.File))
      peptide <- data.frame(Sequence=unique(df_temp2$Sequence))
      
      #extract sample channel
      peptide_S <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Sample"),by="Sequence")  #add heavy channel
      peptide_S$Quan.Channel <- "Sample"
      peptide_S$C_sum[is.na(peptide_S$C_sum)] <-0
      
      #Extract control channel
      peptide_C <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Control"),by="Sequence")  #add light channel
      peptide_C$Quan.Channel <- "Control"
      peptide_C$C_sum[is.na(peptide_C$C_sum)] <-0
      peptide <- rbind(peptide_S,peptide_C)
      
      peptide <- peptide %>% group_by(Sequence)%>%
        summarise(ratio = C_sum[Quan.Channel == "Sample"] /(C_sum[Quan.Channel == "Control"]+C_sum[Quan.Channel == "Sample"]),
                  group ="F1")
      
      #remove NaN
      peptide <- peptide[peptide$ratio!="NaN",]
      
      #combine all ratio
      ratio <- c(ratio,peptide$ratio) 
    }
    
    #Generate Bayesian posterior data
    
    alpha_prior_A <- 2
    beta_prior_A <- 2
    
    # Observed data
    # Replace these with your actual data
    data_A <- as.numeric(ratio[-1])
    data_A <- na.omit(data_A)
    
    # Number of trials and number of successes
    n_trials_A <- length(data_A)
    successes_A <- sum(data_A)
    
    # Posterior distribution parameters
    alpha_posterior_A <- alpha_prior_A + successes_A
    beta_posterior_A <- beta_prior_A + n_trials_A - successes_A
    
    # Generate samples from the posterior distributions
    samples_A <- rbeta(1000, shape1 = alpha_posterior_A, shape2 = beta_posterior_A)
    
    #get mu and sd
    mu <- mean(samples_A)
    sd <- sd(samples_A)
    
    data_container$mu[data_container$accession==acc]<-mu
    data_container$sd[data_container$accession==acc]<-sd
    
  }
  return(data_container)
}

data_container <- na.omit(data_container)
colnames(data_container) <- c("accession","idirt_mu","sd")

data_container0 <- data_container

d <- full_join(data_container,data_container0,by="accession")

d0 <- d

d[is.na(d)]<-0.5

ggplot(data_container,aes(x=idirt_mu))+
  geom_histogram(bins = 20)+
  theme_bw()


con <- Contaminant %>% select(`Uniprot ID`)
colnames(con) <- "accession"
con$contaminant <-"yes"


d1 <- left_join(d,con,by="accession")
d1[is.na(d1)]<-"no"

d1 <- d1[d1$contaminant!="yes",]


ggplot(d1,aes(x=idirt_mu,y=mu))+
  geom_point()+
  theme_bw()

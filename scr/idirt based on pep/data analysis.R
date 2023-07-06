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
  mutate(Quan.Channel0 = if_else((grepl("Heavy", Quan.Channel) & grepl("IDIRT_buffer", Spectrum.File))|
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
    ratio <- ""
    df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
    df_temp <- replace(df_temp,is.na(df_temp),0)  #replace NA to 0
    df_temp <- df_temp %>%                        #group same sequence by file name and quan channel          
      group_by(Sequence,Quan.Channel0,Spectrum.File) %>%
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
        peptide_S <- left_join(peptide, df_temp2 %>% filter(Quan.Channel0=="Sample"),by="Sequence")  #add heavy channel
        peptide_S$Quan.Channel0 <- "Sample"
        peptide_S$C_sum[is.na(peptide_S$C_sum)] <-0
        
        #Extract control channel
        peptide_C <- left_join(peptide, df_temp2 %>% filter(Quan.Channel0=="Control"),by="Sequence")  #add light channel
        peptide_C$Quan.Channel0 <- "Control"
        peptide_C$C_sum[is.na(peptide_C$C_sum)] <-0
        peptide <- rbind(peptide_S,peptide_C)
        
        peptide <- peptide %>% group_by(Sequence)%>%
          summarise(ratio = C_sum[Quan.Channel0 == "Sample"] /(C_sum[Quan.Channel0 == "Control"]+C_sum[Quan.Channel0 == "Sample"]),
                    group =rep)
        
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



d <- get_ratio_0("IDIRT","buffer2_")



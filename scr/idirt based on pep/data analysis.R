library(bayestestR)
library(rstanarm)

set.seed(123456)

# peptides level comapre
source("./scr/idirt based on pep/pep_ratio.R")



# read file
psm <- read.delim("./data/IDIRT data/IDIRTall_PSM_PSMs.txt")

#change Heavy and light to Sample,Control
psm <- psm %>%
  mutate(Quan.Channel = if_else((grepl("Heavy", Quan.Channel) & grepl("IDIRT_buffer", Spectrum.File))|
                                  (grepl("Light", Quan.Channel) & grepl("IDIRT_SWAP", Spectrum.File)), "Sample", "Control"))

#select unique peptide and remove Contaminant and reomve blank in protein accession
psm <- psm %>% filter(!grepl("^.+;.+$", Master.Protein.Accessions, Contaminant=="False")) %>% 
  filter(Master.Protein.Accessions!="")



#idirt data analysis

#1 extract buffer
df <- psm %>% filter(grepl("IDIRT_SWAP_buffer2_",Spectrum.File))

"
#check peptide numbers in protein
df0 <- df %>% select(Master.Protein.Accessions,Sequence,File.ID)%>% 
 group_by(Master.Protein.Accessions,File.ID) %>%
  summarise(C_sum = length(unique(Sequence)))


ggplot(df0,aes(x=File.ID,y=C_sum,fill=File.ID))+
  geom_point()+
  theme_bw()
xx <- df %>% group_by(Master.Protein.Accessions,Sequence) %>% summarise(C_sum = sum(Intensity, na.rm = TRUE))
length(unique)
"


protein_accession <- unique(df$Master.Protein.Accessions)
data_container <- data.frame(accession=unique(df$Master.Protein.Accessions))




for (acc in protein_accession){
  #extract protein
  df_temp <- df %>% filter(Master.Protein.Accessions==acc)
  
  #replace NA to 0
  df_temp <- replace(df_temp,"NA",0)
  #sum intensity of the same peptide
  df_temp <- df_temp %>% 
    group_by(Sequence,Quan.Channel,Spectrum.File) %>%
    summarise(C_sum = sum(Precursor.Abundance, na.rm = TRUE))
  
  #if unique peptide number is less than 2 in all replicate, skip
  if(length(unique(df_temp$Sequence))>1){
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
}

data_container <- na.omit(data_container)
colnames(data_container) <- c("accession","idirt_mu","sd")

data_container0 <- data_container

d <- full_join(data_container,data_container0,by="accession")

d[is.na(d)]<-0.5

ggplot(d,aes(x=idirt_mu,y=mu))+
  geom_point()+
  theme_bw()


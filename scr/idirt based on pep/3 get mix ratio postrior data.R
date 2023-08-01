#get mix ratio

psm <- read.delim("./data/Mix data/mix_PSM_PSMs new_PSMs.txt")

#remove Contaminant and reomve blank in protein accession
psm <- psm %>% filter(Contaminant=="False") 


#change Heavy and light to Sample,Control, #replace NA to 0
#import!!! remove petides not used for quantification
psm <- psm %>% filter(Quan.Channel!="")

psm <- psm %>%
  mutate(group = if_else((grepl("Heavy", Quan.Channel) & grepl("Mix_rep", Spectrum.File))|
                           (grepl("Light", Quan.Channel) & grepl("Mix_swap", Spectrum.File)), "Sample", "Control"))

psm <- replace(psm,is.na(psm),0)

protein <- unique(psm$Master.Protein.Accessions)


get_mix_ratio <- function(mix){  #IDIRT_buffer2_
  #get mix
  
  df <- psm %>% filter(grepl(mix,Spectrum.File))
  df_container <- data.frame(protein)
  #acc <-"P06576"
  for (acc in protein){
    #acc <- "P06732"
    df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
    if(nrow(df_temp)>0){
      #group peptide by sequence,group and file name
      df_temp <- df_temp %>% group_by(Sequence,Spectrum.File,group) %>% 
        summarise(Intensity=sum(Intensity))
      
      #fill missing peptides----
      peptide <- unique(df_temp$Sequence)
      len <- length(peptide)
      file <- unique(df$Spectrum.File)
      df_filled <- data.frame(Sequence=rep(peptide,12),
                              Spectrum.File = rep(rep(file,each=len),2),
                              group=c(rep("Sample",len*6),rep("Control",len*6)))
      
      df_filled <- left_join(df_filled,df_temp,by=c("Sequence","Spectrum.File","group"))
      #repalce NA to 0
      df_filled <- replace(df_filled,is.na(df_filled),0)
      
      #get ratio and remove NaN
      ratio <- df_filled %>% group_by(Sequence,Spectrum.File) %>% 
        summarise(ratio = Intensity[group == "Sample"] /(Intensity[group == "Control"]+Intensity[group =="Sample"]))
      ratio <- ratio[ratio$ratio!="NaN",]
      
      #generate posterior ratio-----
      raw_mu <- median(ratio$ratio)  #the median of row ratio
      raw_sd <- sd(ratio$ratio)
      
      #Generate posterior data from prior data(ratio) 
      alpha_prior <- 2
      beta_prior <- 2
      # Observed data
      # Replace these with your actual data
      # Number of trials and number of successes
      n_trials <- length(ratio$ratio)
      successes <- sum(ratio$ratio)
      
      # Posterior distribution parameters
      alpha_posterior <- alpha_prior + successes
      beta_posterior <- beta_prior + n_trials - successes
      
      # Generate samples from the posterior distributions
      post_ratio <- rbeta(1000, shape1 = alpha_posterior, shape2 = beta_posterior)
      
      #get mu and sd
      mu <- median(post_ratio)
      sd <- sd(post_ratio)
      df_container$ratio[grepl(acc,df_container$protein)]<-mu
      df_container$sd[grepl(acc,df_container$protein)]<-sd
    }
    
  }
  return(df_container)
}


mix <- get_mix_ratio("Mix_rep")
mix_swap <- get_mix_ratio("Mix_swap_rep")


post_mix <- list(mix=mix,
                 mix_swap=mix_swap)

saveRDS(post_mix,"./data/Mix data/post_mix.rsd")

median(mix_swap$ratio)+2*sd(mix_swap$ratio)

ggplot(mix_swap,aes(x=ratio))+
  geom_density()

qqnorm(mix_swap$ratio)

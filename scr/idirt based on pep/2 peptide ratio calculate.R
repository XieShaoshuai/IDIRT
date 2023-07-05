library(ggplot2)
library(dplyr)
library(bayestestR)
library(rstanarm)

set.seed(123456)

# read file
psm <- read.delim("./data/IDIRT data/IDIRTall_PSM_PSMs.txt")
pro <- read.delim("./data/IDIRT data/IDIRTall_Proteins.txt")

#remove Contaminant and reomve blank in protein accession
psm <- psm %>% filter(Contaminant=="False") 

#get potein accession from pro
pro <- pro %>% filter(Protein.FDR.Confidence.Combined == "High",Master == "IsMasterProtein",
                      Number.of.Unique.Peptides >1, Contaminant=="False")

protein <- pro$Accession


#change Heavy and light to Sample,Control, #replace NA to 0----
#import!!! remove petides not used for quantification
psm <- psm %>% filter(Quan.Channel!="")
psm <- psm %>%
  mutate(group = if_else((grepl("Heavy", Quan.Channel) & grepl("IDIRT_buffer", Spectrum.File))|
                                   (grepl("Light", Quan.Channel) & grepl("IDIRT_SWAP", Spectrum.File)), "Sample", "Control"))
psm <- replace(psm,is.na(psm),0)




get_ratio <- function(buffer,protein){  #IDIRT_buffer2_
  #get buffer2
  df <- psm %>% filter(grepl(buffer,Spectrum.File))
  df_container <- data.frame(protein)
  
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
      df_filled <- data.frame(Sequence=rep(peptide,8),
                              Spectrum.File = rep(rep(file,each=len),2),
                              group=c(rep("Sample",len*4),rep("Control",len*4)))
      
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

#get protein ratio
d <- get_ratio("IDIRT_buffer2_",protein)
d2 <- get_ratio("IDIRT_SWAP_buffer2_",protein)

dd <- full_join(d,d2,by="protein")

dd2 <- replace(dd,is.na(dd),0)


ggplot(dd2,aes(x=ratio.x,y=ratio.y,label=protein))+
  geom_point()+
  geom_text_repel(size = 5)+
  xlab("I-DIRT:H/(H+L)")+
  ylab("I-DIRT swap:L/(H+L)")+
  theme_bw()



ggplot(data,aes(median_ratio_1,median_ratio_2,color=group,label=Gene))+
  geom_point(aes(size=average_stochiometry),alpha=0.7)+
  scale_color_manual(values=c("#4daf4a", "#a6d854", "grey50","grey70","lightblue"))+
  geom_text_repel(size = 2.1)+
  xlab("I-DIRT:H/(H+L)")+
  ylab("I-DIRT swap:L/(H+L)")+
  theme_bw()


dd2 %>% filter(ratio.x <0.2, ratio.y>0.8) %>% filter(ratio.x>0)




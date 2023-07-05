pep_ratio <- function(idirt,p){
  df <- psm %>% filter(File.ID==idirt) %>% 
    filter(Master.Protein.Accessions==p)
  
  #sum intensity of each peptide
  df <- df %>%
    group_by(Sequence,Quan.Channel,File.ID) %>%
    summarise(C_sum = sum(Intensity, na.rm = TRUE))
  
  #get the quantified peptides
  peptide <- data.frame(Sequence=unique(df_sum$Sequence))
  
  
  #add channel information
  peptide_H <- left_join(peptide, df_sum %>% filter(Quan.Channel=="Heavy"),by="Sequence")  #add heavy channel
  peptide_H$Quan.Channel <- "Heavy"
  peptide_H$File.ID=idirt
  peptide_H$C_sum[is.na(peptide_H$C_sum)] <-0
  
  peptide_L <- left_join(peptide, df_sum %>% filter(Quan.Channel=="Light"),by="Sequence")  #add light channel
  peptide_L$Quan.Channel <- "Light"
  peptide_L$File.ID=idirt
  peptide_L$C_sum[is.na(peptide_L$C_sum)] <-0
  
  #combine H and L channel
  peptide <- rbind(peptide_H,peptide_L)
  return(peptide)
}

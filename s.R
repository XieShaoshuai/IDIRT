Q02413

df <- psm %>% filter(grepl("buffer2_",Spectrum.File))

df_temp <- df %>% filter(Master.Protein.Accessions=="P12830")
df_temp <- replace(df_temp,"NA",0)
#sum intensity of the same peptide
df_temp <- df_temp %>% 
  group_by(Sequence,Quan.Channel,Spectrum.File) %>%
  summarise(C_sum = sum(Precursor.Abundance, na.rm = TRUE))


replicate <- unique(df_temp$Spectrum.File)

ratio <- data.frame(Sequence="",Spectrum.File="",ratio="")
for(rep in replicate){
  df_temp2 <- df_temp %>%  filter(grepl(rep,Spectrum.File))
  peptide <- data.frame(Sequence=unique(df_temp2$Sequence))
  
  #extract sample channel
  peptide_S <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Sample"),by="Sequence")  #add heavy channel
  peptide_S$Quan.Channel <- "Sample"
  peptide_S$Spectrum.File <- rep
  peptide_S$C_sum[is.na(peptide_S$C_sum)] <-0
  
  #Extract control channel
  peptide_C <- left_join(peptide, df_temp2 %>% filter(Quan.Channel=="Control"),by="Sequence")  #add light channel
  peptide_C$Quan.Channel <- "Control"
  peptide_C$Spectrum.File <- rep
  peptide_C$C_sum[is.na(peptide_C$C_sum)] <-0
  peptide <- rbind(peptide_S,peptide_C)
  
  peptide <- peptide %>% group_by(Sequence,Spectrum.File)%>%
    summarise(ratio = C_sum[Quan.Channel == "Sample"] /(C_sum[Quan.Channel == "Control"]+C_sum[Quan.Channel == "Sample"]))
  peptide <- peptide[peptide$ratio!="NaN",]
  #combine all ratio
  ratio <- rbind(ratio,peptide) 
}

ratio <- ratio[-1,]
ratio$ratio <- as.numeric(ratio$ratio)

ggplot(ratio[-1,], aes(x=Spectrum.File, y=ratio, fill=Spectrum.File)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")

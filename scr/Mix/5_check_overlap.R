#load mix protein data

check_ovelap <- function(mix, path){
  df <- read.delim(path)
  
  df <- df %>% filter(Master=="IsMasterProtein")
  
  idirt_number <- df$Accession
  mix_number <- mix$Accession
  
  
  venn.diagram(
    x = list(idirt_number, mix_number),
    category.names = c("IDIRT" , "MIX"),
    filename = './Output/Mix/venn_Mix_Idirt.png',
    output=TRUE
  )
  
  
  #check ratio in MIX of which protein is found in IDIRT
  idirt <-df[,c("Accession","Description")]
  mix_ratio <- mix[,c("Accession","protein_median_ratio")]
  ratio <- left_join(idirt,mix_ratio,by="Accession")
  
  p<-ggplot(ratio,aes(x=protein_median_ratio))+
    geom_density(fill="#118ab2")+
    ggtitle("Protein ratio in mix of which found in IDRT")+
    theme_bw()
  
  ggsave(p,file="./Output/Mix/ratio of IDIRT protein in Mix.pdf",
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600,)
}

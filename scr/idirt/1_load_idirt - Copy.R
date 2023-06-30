#idirt data analysis

load_idirt <- function(path, col_select,buffer,mix){
  buffer0 <- buffer
  
  df <- read.delim(path)
  df <- df %>%  filter(Protein.FDR.Confidence.Combined == "High"&Master == "IsMasterProtein" &
                         Number.of.Unique.Peptides >1)
  
  #which colname will be selcted for ratio calculates
  metadata <- read.delim(col_select,sep=",")
  col <- metadata %>% filter(select=="Yes")
  df <- df[,col$colname]
  colnames(df) <- col$new_colname
  
  #calculate ratio
  buffer <- as.character(na.omit(col$new_colname[col$Condition==buffer]))
  
  df <- df %>% select(Accession,all_of(buffer)) %>% 
    filter(rowSums(is.na(df[,buffer]))<length(buffer))
  df <- df%>% 
    replace(is.na(df),0)
  
  light_value <- df %>% select(starts_with("Light"))              #extract light channel and heavy channel 
  
  heavy_value <- df %>% select(starts_with("heavy"))    
  
  ratio <- heavy_value/(light_value+heavy_value)                  #calculate protein ratio (H/H+L)
  colnames(ratio) <- paste0("idirt_rep",c(1:ncol(ratio)))
  
  #remove protein abundance information
  #df <- df %>% select(Accession)
  
  #add ratio
  df <- cbind(df,ratio)
  df$median_ratio <- rowMedians(as.matrix(ratio), na.rm=TRUE)
  df$median_light <- rowMedians(as.matrix(light_value), na.rm=TRUE)
  df$median_heavy <- rowMedians(as.matrix(heavy_value), na.rm=TRUE)
  
  #light abundance vs heavy abundance
  
  dd <- df %>% select(starts_with(c("light","heavy")))
  
  
  p <-ggpairs(log2(dd+1), title="Ratio correlation")
  ggsave(p, file=paste0("./Output/IDIRT/Light vs Heavy correlation/correlation_in_Buffer_",buffer0,".pdf"))
  
  
  #ratio missing check
  dd <- df %>% select(starts_with("idirt"))
  dd$missing <- rowSums(is.na(dd))
  
  p<-ggplot(dd)+
    geom_bar(aes(x=missing))+
    ggtitle(paste0("missing ratio in Buffer",buffer))+
    theme_bw()
  
  ggsave(p, file=paste0("./Output/IDIRT/Missing value/missing_ratio_in_Buffer_",buffer0,".pdf"))
  
  #proteins were found in at least 2 replicates
  df <- df %>% filter(rowSums(is.na(df))<3)
  
  
  #venn
  mix_number <- mix$Accession
  idirt_number <- df$Accession
  venn.diagram(
    x = list(idirt_number, mix_number),
    category.names = c("IDIRT" , "MIX"),
    filename = paste0('./Output/IDIRT/Venn mix and idirt/venn_Mix_Idirt',buffer0, ".png"),
    output=TRUE
  )
  
  #check the ratio in Mix
  mix_ratio <- mix %>% select(Accession,median_ratio)
  colnames(mix_ratio) <- c("Accession","Mix_ratio")
  dd <- left_join(df,mix_ratio,by="Accession")
  
  dd <- dd %>% 
    replace(is.na(dd),-0.15) %>% 
    select(Mix_ratio,median_ratio)
  
  dd$group <- "Low"
  dd$group[dd$median_ratio>0.71] <- "High"
  
  p1 <-ggparcoord(dd,
             columns = 1:2,scale="globalminmax", groupColumn = 3
  )+
    geom_hline(yintercept = 0.71,color="blue",linetype = "dashed")+
    theme_bw()
  
  ggsave(p1, file=paste0("./Output/IDIRT/Mix vs idirt/ratio_coordinate_",buffer0,".pdf"))
  
  p2 <- ggplot(dd)+
    geom_histogram(aes(x=median_ratio),col="grey20",fill="red",alpha=0.4,binwidth = 0.025)+
    geom_histogram(aes(x=Mix_ratio),col="grey20",fill="blue",alpha=0.4,binwidth = 0.025)+
    theme_bw()
  ggsave(p2, file=paste0("./Output/IDIRT/ratio_distribution_",buffer0,".pdf"))
  
  p3 <- ggplot(dd)+
    geom_histogram(aes(x=median_ratio),col="grey20",fill="red",alpha=0.4,binwidth = 0.025)+
    theme_bw()
  ggsave(p3, file=paste0("./Output/IDIRT/ratio_distribution2_",buffer0,".pdf"))
  
  
  
  #plot abundance and ratio
  
  # transform the format
  d <- df %>%  select(Accession,median_light,median_heavy,median_ratio) %>%
    filter(median_heavy>0)
  
  light <- d %>% select(median_light) %>% 
    filter(median_light>0)
  
  quartile <- quantile(as.matrix(light))[2]
  
  p4 <-ggplot(data= d %>% filter(median_heavy>0),
         aes(x=reorder(Accession, -median_heavy), y =log2(median_heavy+1)))+
    geom_point(color= "#69b3a2",alpha=0.5)+
    geom_point(data = d %>% filter(median_ratio>=0.71),
               aes(x=Accession,y =log2(median_heavy+1)),color="red")+
    geom_text_repel(data =d %>% filter(Accession=="P12830"),
              aes(label=Accession))+
    geom_hline(yintercept=log2(quartile+1),linetype="dashed", color="red")+
    labs(x="", y = "lgo2(Abundance of Heavy labeled proteins)")+
    theme_bw()+
    theme(panel.grid=element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())
  ggsave(p4, file=paste0("./Output/IDIRT/Abundance and ratio/ratio_distribution2_",buffer0,".pdf"))
  
  return(df)
}

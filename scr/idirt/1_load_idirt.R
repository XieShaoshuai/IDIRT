#idirt data analysis

load_idirt <- function(path, col_select,i,mix,k){
  
  buffer0 <- i
  df <- read.delim(path)
  #selct master protein groups
  df <- df %>%  dplyr::filter(Protein.FDR.Confidence.Combined == "High"&Master == "IsMasterProtein" &
                                Number.of.Unique.Peptides >1 & 
                                Contaminant=="False")
  
  
  #select columns
  metadata <- read.delim(col_select,sep=",")
  col <- metadata %>% dplyr::filter(selectcol=="Yes",experiment==k)
  df <- df[,c('Accession','Description',"Number.of.Peptides",col$colname)]
  
  #rename columns
  colnames(df) <- c('Accession','Description',"Number.of.Peptides",col$new_colname)
  
  
  #select columns contains specific buffer
  i <- as.character(na.omit(col$new_colname[col$Condition==i]))
  df <- df %>% dplyr::select(Accession,Description,Number.of.Peptides,all_of(i))
  
  #sel proteins that NA found in >=1 replicates of sample
  df <- df[(rowSums(df %>% dplyr::select(starts_with("Sample"))>0)>=3),]
  df<- replace(df, is.na(df), 0)
  
  #calculate ratio
  Control_a <- paste0("Control.a.buffer",buffer0)
  Control_b <- paste0("Control.b.buffer",buffer0)
  Control_c <- paste0("Control.c.buffer",buffer0)
  Control_d <- paste0("Control.d.buffer",buffer0)
  
  Sample_a <- paste0("Sample.a.buffer",buffer0)
  Sample_b <- paste0("Sample.b.buffer",buffer0)
  Sample_c <- paste0("Sample.c.buffer",buffer0)
  Sample_d <- paste0("Sample.d.buffer",buffer0)
  
  
  df[,paste0("idirt_",k,"_a")] <- df[,Sample_a]/(df[,Sample_a]+df[,Control_a])
  df[,paste0("idirt_",k,"_b")] <- df[,Sample_b]/(df[,Sample_b]+df[,Control_b])
  df[,paste0("idirt_",k,"_c")] <- df[,Sample_c]/(df[,Sample_c]+df[,Control_c])
  df[,paste0("idirt_",k,"_d")] <- df[,Sample_d]/(df[,Sample_d]+df[,Control_d])
  
  
  #calculate median ratio
  ratio <- df %>% dplyr::select(contains("idirt_"))
  df$median_ratio <- matrixStats::rowMedians(as.matrix(ratio), na.rm=TRUE)
  
  #calculate median abundance stoichiometry
  df$median_Sample <- matrixStats::rowMedians(as.matrix(df %>% dplyr::select(contains("Sample"))), na.rm=TRUE)
  df$median_Control <- matrixStats::rowMedians(as.matrix(df %>% dplyr::select(contains("Control"))), na.rm=TRUE)
  
  control_median <- median(df$median_Control)
  cdh1 <- df %>% filter (Accession=="P12830")
  cdh1 <- (cdh1$median_Sample)/cdh1$Number.of.Peptides
  df$abundance_stoichiometry <- (df$median_Sample)/df$Number.of.Peptides/cdh1
  
  #remove zero in sample
  df <- df %>% filter(median_Sample>0)

  #check the ratio in Mix
  mix_ratio <- mix %>% dplyr::select(Accession,protein_median_ratio)
  colnames(mix_ratio) <- c("Accession","Mix_ratio")
  dd <- left_join(df,mix_ratio,by="Accession")
  
  dd <- dd %>% 
    replace(is.na(dd),0.5) %>% 
    dplyr::select(Mix_ratio,median_ratio)
  
  dd$group <- "Low"
  dd$group[dd$median_ratio>0.71] <- "High"
  
  p1 <-ggparcoord(dd,
             columns = 1:2,scale="globalminmax", groupColumn = 3
  )+
    theme_bw()
  
  ggsave(p1, file=paste0("./Output/IDIRT/ratio_coordinate_",k,"_Buffer",buffer0,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  p2 <- ggplot(dd)+
    geom_histogram(aes(x=median_ratio),col="grey20",fill="red",alpha=0.4,binwidth = 0.05)+
    geom_histogram(aes(x=Mix_ratio),col="grey20",fill="blue",alpha=0.4,binwidth = 0.05)+
    theme_bw()
  ggsave(p2, file=paste0("./Output/IDIRT/ratio_distribution_",k,"_Buffer_",buffer0,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  p3 <- ggplot(dd)+
    geom_histogram(aes(x=median_ratio),col="grey20",fill="red",alpha=0.4,binwidth = 0.05)+
    theme_bw()
  ggsave(p3, file=paste0("./Output/IDIRT/only_idirt_distribution_",k,"Buffer_",buffer0,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  return(df)
}

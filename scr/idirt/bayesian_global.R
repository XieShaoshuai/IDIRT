
bayesian2 <- function(idirt,mix,cutoff_list,k,i){
  idirt_raw <- idirt
  idirt <- idirt[rowSums((idirt %>% dplyr::select(contains("idirt_")))>0)>=3,] #proteins found >3 replicates for Bayesian test
  idirt <- idirt[rowSums(!is.na(idirt %>% dplyr::select(contains("idirt_"))))>=3,]
  df_idirt <- idirt %>% dplyr::select(Accession,Description,abundance_stoichiometry,dplyr::starts_with("idirt"))
  df_mix <- mix %>% dplyr::select(Accession,dplyr::starts_with("protein_ratio"))
  #df <- left_join(df_idirt,df_mix, by="Accession")
  
  df_raw <-df
  #generate random data
  #lower <- cutoff_list$mu - 1.5*cutoff_list$sd
  #upper <- cutoff_list$mu + 1.5*cutoff_list$sd
  
  prior_data <- prior_data <- rnorm(2000,cutoff_list$mu,cutoff_list$sd)

  
  df_clean <- df_idirt
  

  for(k in c(1:nrow(df_clean))){
    posterior_data <- na.omit(as.numeric(df_clean[k,4:7]))
    bf <- ttestBF(x = prior_data, y = posterior_data)
    df_clean$bf[k] <- bf@bayesFactor$bf
  }

  df_clean$median_ratio <- rowMedians(as.matrix(df_clean %>% dplyr::select(starts_with("idirt"))), na.rm=TRUE)  
  
  
  
  p1 <-ggplot(df_clean)+
    geom_histogram(aes(x=median_ratio),fill="red",color="black",alpha=0.5,binwidth = 0.05)+
    geom_histogram(data=df_clean %>% dplyr::filter(bf>3 & median_ratio>0.5), aes(x=median_ratio),fill="#B0CDFE",color="black" ,binwidth = 0.05)+
    ggtitle(paste0("IDIRT_",i,"_Buffer",k))+
    theme_bw()
  ggsave(p1, file=paste0("./Output/Bayesian/ratio_",i,"_Buffer",k,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  #rename colname
  if(i==1){
    colnames(df_clean)=c("Accession","Description","abundance_stoichiometry_1",",idirt_1_a","idirt_1_b","idirt_1_c","idirt_1_d","bf_1","median_ratio_1")
  }else{
    colnames(df_clean)=c("Accession","Description","abundance_stoichiometry_2","idirt_2_a","idirt_2_b","idirt_2_c","idirt_2_d","bf_2","median_ratio_2")
  }
  
  return(df_clean)
}

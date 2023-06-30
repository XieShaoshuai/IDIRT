
bayesian2 <- function(IDIRT,mix,cutoff_list,i){
  df_idirt <- IDIRT %>% dplyr::select(Accession,Description,dplyr::starts_with("idirt"))
  df_mix <- mix %>% dplyr::select(Accession,dplyr::starts_with("protein_ratio"))
  #df <- left_join(df_idirt,df_mix, by="Accession")
  
  df_raw <-df
  #generate random data
  lower <- cutoff_list$mu - 2*cutoff_list$sd
  upper <- cutoff_list$mu + 2*cutoff_list$sd
  
  prior_data <- prior_data <- rnorm(2000,cutoff_list$mu,cutoff_list$sd)

  
  df_clean <- df_idirt
  

  for(k in c(1:nrow(df_clean))){
    posterior_data <- na.omit(as.numeric(df_clean[k,3:6]))
    bf <- ttestBF(x = prior_data, y = posterior_data)
    df_clean$bf[k] <- bf@bayesFactor$bf
  }

  df_clean$median_ratio <- rowMedians(as.matrix(df_clean %>% select(starts_with("idirt"))), na.rm=TRUE)  
  
  
  
  p1 <-ggplot(df_clean)+
    geom_histogram(aes(x=median_ratio),fill="red",color="black",alpha=0.5,binwidth = 0.05)+
    geom_histogram(data=df_clean %>% dplyr::filter(bf>3 & median_ratio>0.5), aes(x=median_ratio),fill="#B0CDFE",color="black" ,binwidth = 0.05)+
    geom_vline(xintercept = 0.72, linetype = "dashed", color="red")+
    theme_bw()
  ggsave(p1, file=paste0("./Output/Bayesian/ratio_histogram_global_",i,".pdf"),
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600)
  
  
  return(df_clean)
}

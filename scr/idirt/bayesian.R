
bayesian <- function(IDIRT,mix,cutoff_list,i){
  df_idirt <- IDIRT %>% dplyr::select(Accession,dplyr::starts_with("Buffer"))
  df_mix <- mix %>% dplyr::select(Accession,dplyr::starts_with("ratio"))
  df <- left_join(df_idirt,df_mix, by="Accession")
  
  df_raw <-df
  #generate random data
  lower <- cutoff_list$mu - 2*cutoff_list$sd
  upper <- cutoff_list$mu + 2*cutoff_list$sd
  
  N <-6
  K <- sum(is.na(df$ratio_2))
  
  df[is.na(df$ratio_2),c("ratio_1","ratio_2","ratio_3",
                         "ratio_4","ratio_5","ratio_6")] <- matrix(runif(K*N, min=lower,max=upper),nrow=K,ncol=N)
  
  
  
  # Select only the rows where "idirt" have values
  idirt_col <- grep("Buffer", colnames(df), value = TRUE)
  rows_idirt <- rowSums(is.na(df[,idirt_col])) == length(idirt_col)
  
  df_clean <- df[!rows_idirt,]
  
  #index_col <-rowMedians(as.matrix(df_clean %>% select(starts_with("Buffer"))), na.rm=TRUE)> 1.0*rowMedians(as.matrix(df_clean %>% select(starts_with("ratio"))), na.rm=TRUE)
  #df_clean <- df_clean[index_col,]
  
  df_tidy <- df_clean %>% 
    pivot_longer(cols = -Accession, names_to = "Experiment", values_to = "Ratio") %>% 
    # add a condition
    mutate(Condition = if_else(grepl("Buffer", Experiment), "IDIRT", "MIX")) %>% 
    # Remove NA values
    na.omit()
  
  df_tidy$Ratio <- log2(df_tidy$Ratio+0.5)

  
  ## Bayessian t-test
  d <- df_tidy %>% 
    nest(data=-Accession) %>% 
    dplyr::mutate(Btest = map(data, ~ttestBF(
      formula = Ratio ~ Condition,
      data = .x)@bayesFactor$bf)) %>% 
    tidyr::unnest(Btest) 
  #%>% dplyr::select(-data)
  
  
  d <- left_join(d,df_clean,by="Accession")
  d$idirt_ratio <- rowMedians(as.matrix(d %>% dplyr::select(starts_with("Buffer"))), na.rm=TRUE)
  d$mix_ratio <- rowMedians(as.matrix(d %>% dplyr::select(starts_with("ratio"))), na.rm=TRUE)
  
  d$color <-"BF <3"
  
  d$color[(d$idirt_ratio/(1-d$idirt_ratio))>2*(d$mix_ratio/(1-d$mix_ratio)) & d$Btest>3] <- "BF>3, log2(FC)>1"
  d$color[d$idirt_ratio<=1.5*d$mix_ratio & d$Btest>3] <- "BF>3"
  
  
  
  p1 <-ggplot(d)+
    geom_histogram(aes(x=idirt_ratio,y = ..count..),fill="red",color="black",alpha=0.5)+
    geom_histogram(data=d %>% dplyr::filter(Btest>3 & idirt_ratio>mix_ratio), aes(x=idirt_ratio),fill="#B0CDFE",color="black")+
    geom_vline(xintercept = 0.68,linetype="dashed",color="red")+
    theme_bw()
  ggsave(p1, file=paste0("./Output/Bayesian/ratio_histogram_BF3_",i,".pdf"))
  
  p2 <- ggplot(d)+
    geom_point(aes(x=idirt_ratio,y = mix_ratio,col=color,size=3), alpha=0.5)+
    geom_vline(xintercept = 0.68, linetype = "dashed", color="red")+
    ylim(0,1)+
    xlim(0,1)+
    theme_bw()
    
  ggsave(p2, file=paste0("./Output/Bayesian/ratio_scatter_BF3_",i,".pdf"))
  
  return(d)
}

cal_ratio <- function(df){
  
  #replace NA with 0
  df<- replace(df, is.na(df), 0)
  
  df <- df %>% dplyr::filter(Grouped_ratio>0.1 & Grouped_ratio<10)
  
  #calculate the median of protein_ratio
  ratio <- df %>% dplyr::select(dplyr::starts_with("ratio"))              #extract ratio of replicate

  ratio$median_ratio <- rowMedians(as.matrix(ratio), na.rm=TRUE)
  
  ratio <- ratio/(ratio+1)               #calculate protein ratio (H/H+L)
  
  df <- df %>% dplyr::select(Accession,Description) #remove protein abundance information
  
  #add ratio
  df <- cbind(df,ratio)
  

  return(df)
}
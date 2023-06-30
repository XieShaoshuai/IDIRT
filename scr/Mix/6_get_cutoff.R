get_cutoff <- function(mix0,i){
  
  ratio <- mix0 %>% dplyr::filter(protein_median_ratio > 0 & protein_median_ratio<1)
  ratio <- ratio$protein_median_ratio
  #stimulate peak distribtion
  x_binned <- mixR::bin(ratio, brks = seq(0, 1, length = 20))
  mod <- mixR::mixfit(x_binned, ncomp = 1,init.method = "kmeans")
  cutoff <- mod[["mu"]] + 2*mod[["sd"]]
  cutoff_list <- list("cutoff"=cutoff,"sd"=mod[["sd"]], mu = mod[["mu"]])
  
  p <- plot(mod)
  ggsave(p,file=paste0("./Output/Mix/mix_distribution_",i,".pdf"),
         width = 5000,
         height = 6000,
         units = "px",
         dpi = 600)
  return(cutoff_list)
}

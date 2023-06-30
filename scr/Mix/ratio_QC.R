ratio_QC <- function(df){
  
  d <- df %>% dplyr::select(contains("ratio"))
  p <-ggpairs(d, title="Ratio correlation") 
  
  ggsave(p,file="./Output/Mix/ratio_correlation.pdf",
         width = 6000,
         height = 6000,
         units = "px",
         dpi = 600,)
  

}

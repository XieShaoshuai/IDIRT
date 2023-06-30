check_missing <- function(df){

  ## the missing values in light or heavy channel
  light_missing <- table(rowSums(is.na(df %>% select(starts_with("Light")))))
  
  heavy_missing <- table(rowSums(is.na(df %>% select(starts_with("Heavy")))))
  
  if (!dir.exists("./Output/Mix/")) dir.create("./Output/Mix/", recursive = TRUE)
  pdf(file="./Output/Mix/0_light missing.pdf")
  barplot(light_missing, 
          col="blue",
          xlab="Missing number in replicate", 
          ylab="Count", 
          main="Light"
  )
  dev.off()
  
  pdf(file="./Output/Mix/0_heavy missing.pdf")
  barplot(heavy_missing,
          col="red",
          xlab="Missing number in replicate",
          ylab="Count",
          main="Heavy")
  dev.off()
  
  cat("The details missing in light (0,1,2,3,4,5,6):", light_missing)
  cat("\n")
  cat("The details missing in heavy (0,1,2,3,4,5,6):", heavy_missing)
  cat("\n")
  cat("batplots were saved in folder ./Output/Mix/")
}

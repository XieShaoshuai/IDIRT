get_interactor <- function(k){
  idirt1  <- load_idirt(idirt_protein_path, col_select_idirt, k, mix1,1)
  idirt2  <- load_idirt(idirt_protein_path, col_select_idirt, k, mix2,2)
  
  
  idirt_1_bayesian <- bayesian2(idirt1,mix1,cutoff_1,k,1)
  idirt_2_bayesian <- bayesian2(idirt2,mix2,cutoff_2,k,2)
  
  data <- full_join(idirt_1_bayesian,idirt_2_bayesian,by=c("Accession","Description"))
  
  data <- stochiometry(data,k)
  
  write.table(data,file=paste0("Buffer ",k,".csv"),sep=",")
  return(data)
}
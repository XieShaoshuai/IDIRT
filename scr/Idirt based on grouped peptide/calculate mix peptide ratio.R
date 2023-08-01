set.seed(84735)

pep <- read.delim("./data/Mix data/Mix_all_PeptideGroups.txt")
pro <- read.delim("./data/Mix data/Mix_all_Proteins.txt")
f_info <- read.delim("./data/Mix data/peptide new name.csv",sep=",")

#select colnames----
sel <- f_info %>% filter(selectcol=="Yes") %>% select(colname,new_colname)
pep <- pep %>% select(sel$colname)
colnames(pep) <- sel$new_colname
rm(sel,f_info)

#get potein accession from pro
pro <- pro %>% filter(Protein.FDR.Confidence.Combined == "High",Master == "IsMasterProtein",
                      Number.of.Unique.Peptides >1, Contaminant=="False")



#data analysis----
#repalce NA with0
pep <- replace(pep,is.na(pep),0)
pep <- pep %>% filter(Contaminant=="False")


pep_ratio <- function(idirt,pep){
  protein <- data.frame(Accession=pro$Accession)  #get protein accession
  df <- pep %>% select(Annotated.Sequence,Master.Protein.Accessions, contains(idirt))
  #reomve 0 in all replicate
  df <- df[rowSums(df[,3:14])>0,]
  for(acc in protein$Accession){
    df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
    if(nrow(df_temp)>1){  #check if protein is found at least 2 peptides in this buffer
      df_temp <- df_temp %>% group_by(Annotated.Sequence,Master.Protein.Accessions) %>%  #group the same peptide
        summarise(across(where(is.numeric),sum),.groups="drop")
      for(rep in c("rep1","rep2","rep3","rep4")){  #check ratio in each replicate
        df_rep <- df_temp %>% select(Annotated.Sequence,contains(rep))
        sample <- df_rep %>% select(contains("Sample"))
        control <- df_rep %>% select(contains("Control"))
        ratio <- sample[,1]/(sample[,1]+control[,1])
        ratio <- data.frame(ratio=ratio[ratio[,1]!="NaN",])
        
        #Observed data
        #Number of trials and number of successes
        
        successes <- sum(ratio[,1])
        n_trials <- length(ratio[,1])
        #Observed data
        #plot_beta_binomial(alpha = 2, beta = 2, y = successes, n = n_trials)
        y <- summarize_beta_binomial(alpha = 2, beta = 2, y = successes, n = n_trials)
        #Posterior distribution parameters
        
        #add ratio to protein----
        protein[protein$Accession==acc,rep] <- y$mean[2]
      }
    }
  }
  
  protein$mean <- apply(protein[, c("rep1", "rep2", "rep3","rep4")], 1, function(x) mean(x, na.rm = TRUE))
  
  return(protein)
}

#get ratio
mix <- pep_ratio("Mix",pep)
mix_swap <- pep_ratio("Swap",pep)


mix_data <- list(mix=mix,mix_swap=mix_swap)

saveRDS(mix_data,file="./data/Mix data/mix_07.rsd")



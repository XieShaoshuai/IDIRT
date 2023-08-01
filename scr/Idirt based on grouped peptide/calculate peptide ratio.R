#load data----
library(bayesrules)
set.seed(84735)
pep <- read.delim("./data/IDIRT data/IDIRT_pep_PeptideGroups.txt")
pro <- read.delim("./data/IDIRT data/IDIRTall_Proteins.txt")
f_info <- read.delim("./data/IDIRT data/peptide new name.csv",sep=",")

#select colnames----
sel <- f_info %>% filter(selectcol=="Yes") %>% select(colname,new_colname)
pep <- pep %>% select(sel$colname)
colnames(pep) <- sel$new_colname
rm(sel,f_info)

#get potein accession from pro
pro <- pro %>% filter(Protein.FDR.Confidence.Combined == "High",Master == "IsMasterProtein",
                      Number.of.Unique.Peptides >1)



#data analysis----
#repalce NA with0
pep <- replace(pep,is.na(pep),0)
pep <- pep %>% filter(Contaminant=="False")


pep_ratio <- function(idirt,pep){
  protein <- data.frame(Accession=pro$Accession)  #get protein accession
  df <- pep %>% select(Annotated.Sequence,Master.Protein.Accessions, contains(idirt))
  #reomve 0 in all replicate
  df <- df[rowSums(df[,3:10])>0,]
  for(acc in protein$Accession){
    df_temp <- df %>% filter(grepl(acc,Master.Protein.Accessions))
    if(nrow(df_temp)>1){  #check if protein is found at least 2 peptides in this buffer
      df_temp <- df_temp %>% group_by(Annotated.Sequence,Master.Protein.Accessions) %>%  #group the same peptide
        summarise(across(where(is.numeric),sum),.groups="drop")
      for(rep in c("a.","b.","c.","d.")){  #check ratio in each replicate
        df_rep <- df_temp %>% select(Annotated.Sequence,contains(rep))
        sample <- df_rep %>% select(contains("Sample."))
        control <- df_rep %>% select(contains("Control"))
        ratio <- sample[,1]/(sample[,1]+control[,1])
        ratio <- data.frame(ratio=ratio[ratio[,1]!="NaN",])
        
        #generate posterior ratio using beta bayesian----
     
        successes <- sum(ratio[,1])
        n_trials <- length(ratio[,1])
        #Observed data
        #plot_beta_binomial(alpha = 2, beta = 2, y = successes, n = n_trials)
        y <- summarize_beta_binomial(alpha = 2, beta = 2, y = successes, n = n_trials)
        #Posterior distribution parameters
     
        posterior <- rbeta(100,y$alpha[2],y$beta[2])
        #add ratio to protein----
        protein[protein$Accession==acc,rep] <- median(posterior)
        protein[protein$Accession==acc,paste0(rep,"raw")] <- mean(ratio[,1]) 
      }
    }
  }
  
  protein$mean<- apply(protein[, c("a.", "b.", "c.","d.")], 1, function(x) mean(x, na.rm = TRUE))
  
  return(protein)
}

#get ratio
IDIRT_buffer2 <- pep_ratio("IDIRT_buffer02",pep)
IDIRT_buffer3 <- pep_ratio("IDIRT_buffer3",pep)
IDIRT_buffer16 <- pep_ratio("IDIRT_buffer16",pep)
IDIRT_buffer19 <- pep_ratio("IDIRT_buffer19",pep)
IDIRT_buffer25 <- pep_ratio("IDIRT_buffer25",pep)
IDIRT_buffer29 <- pep_ratio("IDIRT_buffer29",pep)

Swap_buffer2 <- pep_ratio("Swap_buffer02",pep)
Swap_buffer3 <- pep_ratio("Swap_buffer3",pep)
Swap_buffer16 <- pep_ratio("Swap_buffer16",pep)
Swap_buffer19 <- pep_ratio("Swap_buffer19",pep)
Swap_buffer25 <- pep_ratio("Swap_buffer25",pep)
Swap_buffer29 <- pep_ratio("Swap_buffer29",pep)

#assign gene name----
#load gene name
gene <- read.delim("./data/IDIRT data/Gene name.csv",sep=",")
IDIRT_buffer2 <- left_join(IDIRT_buffer2,gene,by="Accession")
IDIRT_buffer3 <- left_join(IDIRT_buffer3,gene,by="Accession")
IDIRT_buffer16 <- left_join(IDIRT_buffer16,gene,by="Accession")
IDIRT_buffer19 <- left_join(IDIRT_buffer19,gene,by="Accession")
IDIRT_buffer25 <- left_join(IDIRT_buffer25,gene,by="Accession")
IDIRT_buffer29 <- left_join(IDIRT_buffer29,gene,by="Accession")


#generate list
IDIRT_data <- list("buffer2"=IDIRT_buffer2,"buffer2_swap"=Swap_buffer2,
                   "buffer3"=IDIRT_buffer3,"buffer3_swap"=Swap_buffer3,
                   "buffer16"=IDIRT_buffer16,"buffer16_swap"=Swap_buffer16,
                   "buffer19"=IDIRT_buffer19,"buffer19_swap"=Swap_buffer19,
                   "buffer25"=IDIRT_buffer25,"buffer25_swap"=Swap_buffer25,
                   "buffer29"=IDIRT_buffer29,"buffer29_swap"=Swap_buffer29)



saveRDS(IDIRT_data,file="./data/IDIRT data/IDIRT_peptide_ratio_07.rsd")


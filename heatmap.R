library('ComplexHeatmap')
library(dendextend)
#interactors analysis
gProfile <- read.csv("./gProfile.csv")
colnames(gProfile) <- c("Accession","Gene")
gProfile <- gProfile[!duplicated(gProfile$Accession),]

df <- data
df <- df %>% replace(is.na(df),0)

df <- left_join(df,gProfile,by="Accession")


df$hit2 <- df$`Buffer2_Idirt ratio`>0.5 & df$`Buffer2_Bayes Factor`>3
df$hit3 <- df$`Buffer3_Idirt ratio`>0.5 & df$`Buffer3_Bayes Factor`>3
df$hit16 <- df$`Buffer16_Idirt ratio`>0.5 & df$`Buffer16_Bayes Factor`>3
df$hit19 <- df$`Buffer19_Idirt ratio`>0.5 & df$`Buffer19_Bayes Factor`>3
df$hit25 <- df$`Buffer25_Idirt ratio`>0.5 & df$`Buffer25_Bayes Factor`>3
df$hit29 <- df$`Buffer29_Idirt ratio`>0.5 & df$`Buffer29_Bayes Factor`>3

df$hit <- rowSums(df[,c("hit2","hit3","hit16","hit19","hit25","hit29")])

df_6 <- df %>% filter(hit>=6) %>% 
  dplyr::select(Accession,Gene,contains("ratio"))

df_6 <- df_6[!duplicated(df_6$Gene),]

rownames(df_6) <- df_6$Gene

df_6 <-df_6 %>% 
  dplyr::select(contains("ratio"))


# Clusterisation using 3 variables
dend <-df_6 %>% 
  dist(method="maximum") %>% 
  hclust() %>% 
  as.dendrogram()
dend <-df_6 %>% 
  dist() %>% 
  hclust(method = "average") %>% 
  as.dendrogram() %>% 
  set("branches_k_color", value = c("skyblue", "orange"), k = 2) %>% 
  plot(horiz=FALSE, axes=TRUE)


'plot(data=df_6[c("KRT1",
                   "CTNBB1",
                   "AZGP1",
                   "SLC7A5",
                   "SLC3A2")])






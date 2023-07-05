library(ggplot2)
library(dplyr)

#load data----
pro <- read.delim("./data/IDIRT data/IDIRTall_Proteins.txt")
pep <- read.delim("./data/IDIRT data/IDIRT_pep_PeptideGroups.txt")

#data clean
pro <- pro %>% filter(Protein.FDR.Confidence.Combined == "High",Master == "IsMasterProtein",
                      Number.of.Unique.Peptides >1, Contaminant=="False")

pro <- replace(pro,is.na(pro),0)  #replace NA with0
pep <- replace(pep,is.na(pep),0)

#compare protein abundance with peptide abundance in F1(IDIRT buffe2), F25(IDIRT swap buffer2)
f <- "F1"
f.light <- paste0("Abundance.",f,".Light.Sample")
f.heavy <- paste0("Abundance.",f,".Heavy.Sample")

protein <- pro %>% select(Accession,contains(c(f.light,f.heavy))) %>% 
  filter(Accession!="")
colnames(protein) <- c("Accession","Protein.Light","Protein.Heavy")
protein <- protein[rowSums(protein[,2:3])>0,]  #reomove protein not found in file

#sum peptide total abundance
peptide0 <- pep %>% select(Master.Protein.Accessions,Annotated.Sequence,f.light,f.heavy) %>% 
  filter(Master.Protein.Accessions!="")

colnames(peptide0) <- c("Accession","sequence","Peptide.Light","Peptide.Heavy")
peptide <- peptide0%>% group_by(Accession) %>% 
  summarise(Abundance.peptide.Light=sum(Peptide.Light,na.rm=TRUE),
         Abundance.peptide.Heavy=sum(Peptide.Heavy,na.rm=TRUE))

#bind protein ana peptide data
data <- left_join(protein,peptide,by="Accession")

data[,2:5] <- log2(data[,2:5]+1)

#plot
ggplot(data)+
  geom_point(aes(x=Protein.Light,y=Abundance.peptide.Light,color="blue",alpha=0.5))+
  geom_point(aes(x=Protein.Heavy,y=Abundance.peptide.Heavy,color="red",alpha=0.5))+
  theme_bw()


#check P42704
d <- peptide0 %>% filter(Accession=="P42704")
d$ratio <- d$Peptide.Heavy/(d$Peptide.Heavy+d$Peptide.Light)
d <- d %>% filter(ratio!="NaN")

protein_ratio <- sum(d$Peptide.Heavy)/(sum(d$Peptide.Heavy+d$Peptide.Light))

#P42704 check
p <- ggplot(d)+
  geom_segment(aes(x=sequence, xend=sequence, 
                   y=log2(Peptide.Light+1), yend=log2(Peptide.Heavy+1)), 
               color="grey",linetype="dotdash",size=0.6 )+
  geom_point(aes(x=sequence,y=log2(Peptide.Light+1)),color="blue",size=3,alpha=0.5)+
  geom_point(aes(x=sequence,y=log2(Peptide.Heavy+1)),color="red",size=3,alpha=0.5)+
  ggtitle("LRPPRC in IDIRT")+
  ylab("log2(Abundance)")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90,size=5.5,hjust=1))

  



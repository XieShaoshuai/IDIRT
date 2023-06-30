
stochiometry <- function(data){
  #idirt data analysis
  
  
  data <- na.omit(data)
  #gene name convert
  gene <- read.delim("./data/IDIRT data/Gene name.csv",sep=",")
  data<- left_join(data,gene, by="Accession")
  
  
  
  
  data$average_stochiometry <- (data$abundance_stoichiometry_1+data$abundance_stoichiometry_2)/2
  
  data$group <- "Background"
  
  cutoff1 <- cutoff_1$mu-cutoff_1$sd
  cutoff2 <- cutoff_2$mu-cutoff_2$sd
  
  strong_interaction <- data$bf_1>3 & data$median_ratio_1 >0.5 & data$median_ratio_2 >0.5 & data$bf_2>3
  random_hit1 <- data$bf_1>3 & data$median_ratio_1 <0.5 & data$median_ratio_2>0.5 & data$bf_2>3
  random_hit2 <- data$bf_1>3 & data$median_ratio_1 >0.5 & data$median_ratio_2<0.5 & data$bf_2>3
  bg <- data$bf_1<3 & data$bf_2 <3
  
  weak_interaction <- (data$bf_1>3 & data$median_ratio_1 >0.5 & data$median_ratio_2 > cutoff2 & data$bf_2<3)|
    (data$bf_1<3 & data$median_ratio_1 >cutoff1 & data$median_ratio_2 > 0.5 & data$bf_2>3)
  
  
  data$group[strong_interaction] <- "Stable interactor"
  data$group[random_hit1] <- "Random hit"
  data$group[random_hit2] <- "Random hit"
  data$group[weak_interaction] <- "Unstable interactor"
  
  data$label<-""
  data$label[data$group=="Stable interactor"] <- data$Gene[data$group=="Stable interactor"]
  data$group[data$Accession=="P12830"] <- "Bait"
  
  
  data$group <- factor(data$group, levels=c('Stable interactor', 'Unstable interactor', 'Random hit', 'Background',"Bait"))
  
  
  p<- ggplot(data,aes(median_ratio_1,median_ratio_2,color=group,label=label))+
    geom_point(aes(size=average_stochiometry),alpha=0.7)+
    scale_color_manual(values=c("#4daf4a", "#a6d854", "grey50","grey70","lightblue"))+
    geom_text_repel(size = 2)+
    xlab("I-DIRT:H/(H+L)")+
    ylab("I-DIRT swap:L/(H+L)")+
    theme_bw()
  ggsave(p, file=paste0("./Output/IDIRT/stochiomerty_",i,"_Buffer",k,".pdf"),
         width = 6,
         height = 4,
         units = "in",
         dpi = 600)
  
  p1<- ggplot(data,aes(median_ratio_1,median_ratio_2,color=group,label=Gene))+
    geom_point(aes(size=average_stochiometry),alpha=0.7)+
    scale_color_manual(values=c("#4daf4a", "#a6d854", "grey50","grey70","lightblue"))+
    geom_text_repel(size = 2.1)+
    xlab("I-DIRT:H/(H+L)")+
    ylab("I-DIRT swap:L/(H+L)")+
    theme_bw()
  ggsave(p1, file=paste0("./Output/IDIRT/stochiomerty_all_",i,"_Buffer",k,".pdf"),
         width = 7,
         height = 5,
         units = "in",
         dpi = 600)
}




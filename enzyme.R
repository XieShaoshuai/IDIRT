
library(dplyr)

data$Group <- factor(data$Group, levels=c('H2O', '0.5CSA', '1CSA', '2CSA'))


ggplot(data,aes(x=Enzyme,y=Value,fill=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
  


variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)

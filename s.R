


ggplot(data)+
  geom_point(aes(x=median_ratio_1,y=median_ratio_2, color= bf_1 > 3 & bf_2 >3 & median_ratio_1>0.5 & median_ratio_2>0.5))+
  theme_bw()


data$group <- "Backgroud"

strong_index <- data$bf_1>3 & data$bf_2>3 &data$median_ratio_1>0.5 & data$median_ratio_2>0.5
noise <- data$bf_1>3 & data$bf_2>3 &data$median_ratio_1<0.5 & data$median_ratio_2<0.5
random1 <- data$bf_1>3 & data$median_ratio_1>0.5 & data$median_ratio_2<0.5
random2 <- data$bf_2>3 & data$median_ratio_2>0.5 & data$median_ratio_1<0.5
weak_index1 <- data$bf_1>3 & data$median_ratio_1>0.5 & data$bf_2<=3
weak_index2 <- data$bf_2>3 & data$median_ratio_2>0.5 &  data$bf_1<=3

data$group[strong_index] <- "Strong interaction"
data$group[noise] <- "noise"
data$group[random1] <- "random interaction"
data$group[random2] <- "random interaction"
data$group[weak_index1] <- "weak interaction"
data$group[weak_index2] <- "weak interaction"

ggplot(data)+
  geom_point(aes(x=median_ratio_1,y=median_ratio_2, color=group,size=2,alpha=0.5))+
  theme_bw()

library(ellipse)

library(ggplot2)

idirt <- readRDS("./data/IDIRT data/IDIRT_peptide_ratio_07.rsd")
mix <- readRDS("./data/Mix data/mix_07.rsd")


mix_all <- data.frame(x=mix$mix$mean,y=mix$mix_swap$mean)
mix_all <- na.omit(mix_all)

#plot mix ellipse with 95% confidence
d <- buffer2 %>% filter(mean.x>cutoff1_1,buffer2$mean.y>cutoff2_1) %>% 
  select(mean.x,mean.y)

center <- unname(colMeans(mix_all))
# Calculate the covariance matrix
cov_matrix <- cov(mix_all)
# Calculate the eigenvalues and eigenvectors of the covariance matrix
eig <- eigen(cov_matrix)

ellipse <- ellipse(cov_matrix, centre = center,level=0.95)
ellipse(cov_matrix, centre = center, level = 0.95)
# Convert the ellipse to a data frame for plotting
ellipse_df <- data.frame(ellipse)

#check the cutoff
cutoff1 <- mean(mix_all$x)+ 1*sd(mix_all$x)
cutoff2 <- mean(mix_all$y)+ 1*sd(mix_all$y)
cutoff1_1 <- mean(mix_all$x) - 1*sd(mix_all$x)
cutoff2_1 <- mean(mix_all$y)- 1*sd(mix_all$y)


buffer2 <- full_join(idirt$buffer2[,c("Accession","Gene","median")],
                     idirt$buffer2_swap[,c("Accession","median")], by="Accession")

buffer3 <- full_join(idirt$buffer3[,c(1,11,10)],
                     idirt$buffer3_swap[,c(1,10)], by="Accession")

buffer16 <- full_join(idirt$buffer16[,c(1,11,10)],
                      idirt$buffer16_swap[,c(1,10)], by="Accession")

buffer19 <- full_join(idirt$buffer19[,c(1,11,10)],
                      idirt$buffer19_swap[,c(1,10)], by="Accession")

buffer25 <- full_join(idirt$buffer25[,c(1,11,10)],
                      idirt$buffer25_swap[,c(1,10)], by="Accession")

buffer29 <- full_join(idirt$buffer29[,c(1,11,10)],
                      idirt$buffer29_swap[,c(1,10)], by="Accession")


for (i in nrow(buffer2)){
  point <- c(x = buffer2[i,3], y = buffer2[i,4])
  # Calculate the difference between the point and the center
  diff <- matrix(point - center, ncol = 1)
  # Calculate the Mahalanobis distance
  mahalanobis_dist <- t(diff) %*% solve(cov_matrix) %*% diff
  
  # Check if the point is inside the ellipse
  inside <- mahalanobis_dist <= qchisq(0.95, 2)
  buffer2$inside[i]<-inside
}


weak <- buffer2$mean.x>cutoff1_1& buffer2$mean.y>cutoff2_1&buffer2$inside=="FALSE"
strong <- buffer2$mean.x>cutoff1 & buffer2$mean.y>cutoff2 & buffer2$inside=="FALSE"
buffer2$group<-"Noise"
buffer2$group[weak] <- "Weak interaction"
buffer2$group[strong] <- "Strong interaction"


buffer3$group <- "Noise"
buffer3$group[buffer3$mean.x>cutoff1&buffer3$mean.y>cutoff2] <- "Strong interaction"
buffer3$group[buffer3$mean.x>cutoff1&buffer3$mean.y>cutoff2_1&
                buffer3$mean.y<cutoff2] <- "Weak interaction"
buffer3$group[buffer3$mean.x<cutoff1&buffer3$mean.x>cutoff1_1&
                buffer3$mean.y>cutoff2] <- "Weak interaction"

buffer16$group <- "Noise"
buffer16$group[buffer16$mean.x>cutoff1&buffer16$mean.y>cutoff2] <- "Strong interaction"
buffer16$group[buffer16$mean.x>cutoff1&buffer16$mean.y>cutoff2_1&
                buffer16$mean.y<cutoff2] <- "Weak interaction"
buffer16$group[buffer16$mean.x<cutoff1&buffer16$mean.x>cutoff1_1&
                buffer16$mean.y>cutoff2] <- "Weak interaction"


buffer19$group <- "Noise"
buffer19$group[buffer19$mean.x>cutoff1&buffer19$mean.y>cutoff2] <- "Strong interaction"
buffer19$group[buffer19$mean.x>cutoff1&buffer19$mean.y>cutoff2_1&
                buffer19$mean.y<cutoff2] <- "Weak interaction"
buffer19$group[buffer19$mean.x<cutoff1&buffer19$mean.x>cutoff1_1&
                buffer19$mean.y>cutoff2] <- "Weak interaction"

buffer25$group <- "Noise"
buffer25$group[buffer25$mean.x>cutoff1&buffer25$mean.y>cutoff2] <- "Strong interaction"
buffer25$group[buffer25$mean.x>cutoff1&buffer25$mean.y>cutoff2_1&
                buffer25$mean.y<cutoff2] <- "Weak interaction"
buffer25$group[buffer25$mean.x<cutoff1&buffer25$mean.x>cutoff1_1&
                buffer25$mean.y>cutoff2] <- "Weak interaction"

buffer29$group <- "Noise"
buffer29$group[buffer29$mean.x>cutoff1&buffer29$mean.y>cutoff2] <- "Strong interaction"
buffer29$group[buffer29$mean.x>cutoff1&buffer29$mean.y>cutoff2_1&
                buffer29$mean.y<cutoff2] <- "Weak interaction"
buffer29$group[buffer29$mean.x<cutoff1&buffer29$mean.x>cutoff1_1&
                buffer29$mean.y>cutoff2] <- "Weak interaction"



ratio_all <- list(buffer2,buffer3,buffer16,buffer19,buffer25,buffer29)
saveRDS(ratio_all,"./data/IDIRT data/ratio_all_07.rsd")

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 16))


data <- buffer2
title <- "Condition 2"


# Check if the point is inside the ellipse

ggplot()+
  geom_point(data = data %>% filter (group!="Noise"),aes(x=log2(mean.x+0.5),y=log2(mean.y+0.5),fill=group,color=group),alpha=0.5,size=3,shape=16)+
  geom_point(data = data %>% filter (group =="Noise"),aes(x=log2(mean.x+0.5),y=log2(mean.y+0.5),fill=group,color=group),alpha=0.2,size=2,shape=16)+
  #geom_path(data = ellipse_df, aes(x, y), color = 'red')+
  scale_color_manual(values = c("grey40", "#35A165", "#0073C2FF"))+
  geom_text_repel(data=data %>% filter(group=="Strong interaction"),aes(x=mean.x,y=mean.y,label=Gene),size = 2,alpha=0.7,color="black")+
  xlim(-1,1)+
  ylim(-1,1)+
  xlab("log2(ratio+0.5)")+
  ylab("log2(ratio+0.5)")+
  ggtitle(title)+
  theme_bw()+
  My_Theme


pdf(" ")
  ylim(0,1)+
  geom_path(data = ellipse_df, aes(x, y), color = 'red')

dataEllipse(data[,1],data[,2],levels=0.95,col="black",xlim=c(-1,1),ylim=c(-4,2),xlab="DNMT3A",ylab="PGLYRP2",lwd=1,
            grid=F,pch=20,lty=2)
dev.off()



ggplot(mix_all, aes(x, y)) +
  geom_point() +
  xlim(0,1)+
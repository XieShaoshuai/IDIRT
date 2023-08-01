#plot data

#read IDIRT data

idirt <- readRDS("./data/IDIRT data/IDIRT_ratio.rsd")
mix <- readRDS("./data/Mix data/post_mix.rsd")




mix_cutoff1 <- mean(mix$mix$ratio) + 1.5*sd(mix$mix$ratio)
mix_cutoff2 <- mean(mix$mix$ratio) - 1.5*sd(mix$mix$ratio)

swap_cutoff1 <- mean(mix$mix_swap$ratio) + 1.5*sd(mix$mix_swap$ratio)
swap_cutoff2 <- mean(mix$mix_swap$ratio) - 1.5*sd(mix$mix_swap$ratio)

buffer16 <- full_join(idirt$buffer16,idirt$buffer16_swap,by="protein")

buffer16_filter <- buffer16 %>% filter(!(ratio.x<mix_cutoff2 & ratio.y>swap_cutoff1))


ggplot(buffer29_filter,aes(x=ratio.x,y=ratio.y,label=Gene.x))+
  geom_point(alpha=0.7)+
  geom_text_repel(size = 4)+
  xlim(0,1)+
  ylim(0,1)+
  xlab("I-DIRT:H/(H+L)")+
  ylab("I-DIRT swap:L/(H+L)")+
  theme_bw()

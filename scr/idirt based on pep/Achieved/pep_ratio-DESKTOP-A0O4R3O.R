# peptides level comapre
source("./scr/idirt based on pep/pep_ratio.R")




# read psm file
psm <- read.delim("./data/IDIRT data/IDIRTall_PSM_PSMs.txt")

p <- "P42704"

d1 <- pep_ratio("F1",p) %>% group_by(Sequence)%>%
  summarise(ratio = C_sum[Quan.Channel == "Heavy"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),
            group ="F1")

d2 <-  pep_ratio("F2",p)%>% group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Heavy"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),
            group ="F2")

d3 <- pep_ratio("F3",p)%>% group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Heavy"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F3")

d4 <- pep_ratio("F3",p)%>% group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Heavy"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F4")

d1_1 <- pep_ratio("F25",p)%>% group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Light"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F25")

d2_1  <- pep_ratio("F26",p) %>%group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Light"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F26")

d3_1 <- pep_ratio("F27",p) %>%group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Light"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F27")

d4_1 <- pep_ratio("F28",p) %>%group_by(Sequence) %>%
  summarise(ratio = C_sum[Quan.Channel == "Light"] /(C_sum[Quan.Channel == "Light"]+C_sum[Quan.Channel == "Heavy"]),group ="F28")



plot_df <- rbind(d1,d2,d3,d4,
                 d1_1,d2_1,d3_1,d4_1)


#calculate ratio

p1 <- ggplot(plot_df, aes(x=group, y=ratio, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(width=0.5,alpha=0.2)+
  geom_boxplot(width=0.2)+
  theme_bw()
p1


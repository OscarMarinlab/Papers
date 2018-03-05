Data<-read.csv("2017_PV_electro_virus.csv", header= TRUE)


library(ggplot2)
dodge <- position_dodge(width=0.9)
limits <- aes(ymax = Freq_M + error, ymin = Freq_M- error)
gp1a <- ggplot(Data, aes(x= condition, y=Freq_M)) + 
  geom_bar(stat="identity", position=dodge, colour = "black",  fill = c("light gray")) +
  geom_errorbar (limits, position=dodge, width = 0.2)+
  theme_classic (base_size = 18) +  labs (x = 'conditions', y = 'Fraction L5/L6 Translaminer PV') 


L5.df = matrix(c(22,17,8,52), nrow = 2)
fisher.test(L23.df)


Data<-read.csv("2016electrovsvirus.csv", header= TRUE)


library(ggplot2)
dodge <- position_dodge(width=0.9)
limits <- aes(ymax = Freq_M + error, ymin = Freq_M- error)
gp1a <- ggplot(Data, aes(x= condition, y=Freq_M)) + 
  geom_bar(stat="identity", position=dodge, colour = "black",  fill = c("light gray")) +
  geom_errorbar (limits, position=dodge, width = 0.2)+
  theme_classic (base_size = 18) +  labs (x = 'conditions', y = 'Fraction of Martinotti cells') 


L23.df = matrix(c(36,4,18,12), nrow = 2)
fisher.test(L23.df)


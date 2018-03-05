
library (plotrix) 
library (ggplot2) 
library (plyr)

DATA <-read.csv ("2016-06-23_meanL1axon.csv")


mm6<- ddply(DATA, .(Genotype, dist), summarise, m_I= mean(Imean), se = std.error(Imean)/2)


gp1 <- ggplot(mm6, aes(x = factor(dist), y = m_I, group = Genotype, colour=Genotype, fill = Genotype)) + 
  geom_errorbar(aes(ymin=m_I-se, ymax=m_I+se),
                width=.1,colour="light gray") +  # Width of the error bars
  geom_line() +
  scale_colour_manual(values=c("#555555","#FF6600"))+
  theme_classic (base_size = 12) +
  scale_x_discrete(breaks = c("1", "100", "300", "400", "600")) +
  labs (x = 'dist', y = 'mean int') 


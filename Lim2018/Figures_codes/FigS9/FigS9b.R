library (plotrix) 
library (ggplot2) 
library (plyr)

DATA <-read.csv("2016-06-15_non_MC.csv")
DATA$Genotype <- factor(DATA$Genotype, levels = c("Mafb-fl/+", "Mafb-fl/fl"))
mm3 <- DATA

gp4d <- ggplot(mm3, aes(factor(Genotype), Length_n_t, colour = Genotype)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 2, position=position_jitter(width=0.1))+
  scale_fill_manual(values=c("#555555","orange"))+
  scale_colour_manual(values=c("#555555","orange"))+
  theme_classic (base_size = 16) +
  xlab("Genotype") +
  ylab("total neurite")+
  ylim(0,50)

t2 <- t.test(Length_n_t~Genotype, data=DATA)


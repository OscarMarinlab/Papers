library (plotrix) 
library (ggplot2) 
library (plyr)

DATAall <-read.csv("L1_L3_ROI.csv")
DATAall$genotype <- factor(DATAall$genotype, levels = c("Mafb+/fl", "Mafb-fl/fl"))
mm2<- ddply(DATAall, . (genotype, area), summarise, mINT= mean(meanINT), se = std.error(meanINT))
mm3<- ddply(DATAall, . (genotype, brain, area), summarise, bmINT= mean(meanINT) )
mm4<- ddply(mm3, .(genotype, area), summarise, mmINT= mean(bmINT), se = std.error(bmINT))



stat1<-aov(meanINT ~ genotype * area, data = DATAall)
stat1t<- TukeyHSD(stat1)
q1 <- ddply(DATAall, .(genotype, area), summarise, q25 = quantile (meanINT, 0.25), q50 = quantile (meanINT, 0.5),
            q75 = quantile (meanINT, 0.75))

#boxplot
fig5h <- ggplot(DATAall, aes(x = factor(genotype), y = meanINT)) + 
  facet_grid(area ~.  , scale = "free_y") +
  geom_boxplot(aes(colour = genotype), outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper.., colour = genotype), outlier.colour  = NA) + 
  scale_fill_manual(values=c("gray","orange"))+
  scale_colour_manual(values=c("gray","orange"))+
  ylab("density")+
  theme_classic()

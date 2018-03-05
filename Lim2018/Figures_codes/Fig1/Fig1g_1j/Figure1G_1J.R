library (plotrix) 
library (ggplot2) 
library (plyr)

DATA <-read.csv("E14toE16allgeno.csv", header = TRUE)
DATA$Genotype <- factor(DATA$Genotype, levels = c("Nkx2.1", "Gad65-gfp", "SST", "Dlx1/2"))
mm1 <- ddply(DATA, .(Genotype, Age, Field, Brain.no), summarise, FMZ= MZ/(MZ+SVZ))
mm2 <- ddply(mm1, .(Genotype, Age, Brain.no), summarise, avgFMZ= mean(FMZ))
mm3 <- ddply (mm2, .(Genotype, Age), summarise, FractionMZ = mean(avgFMZ), Sem = std.error(avgFMZ))

mm3

#plot for figure 1G
mm1cdf <-ddply(mm1,.(Genotype), transform, ecdall = ecdf(FMZ)(FMZ)) 
cdfall <- ggplot(mm1cdf, aes(x=FMZ)) + 
  stat_ecdf(aes(colour=Genotype), size = 1) +
  scale_color_brewer(palette="Set1") +
  theme_classic()

cdfall

#computing p-value with 2 way ANOVA
stats1 <-aov(avgFMZ ~ Genotype*Age, data=mm2)
stats1t <-TukeyHSD(stats1)

statsdf = as.data.frame(do.call(rbind, stats1t)) #p-value for figure 1J is in datafram statsdf

#Figure 1J in bar graph
pMZ <- ggplot(mm3, aes(x=Age, y=FractionMZ, fill=Genotype)) + 
  scale_fill_brewer(palette="Set1") +
  geom_bar(position=position_dodge(0.8), stat="identity", width = 0.8) +
  geom_errorbar(aes(ymin=FractionMZ-Sem, ymax=FractionMZ+Sem),
                width=0.1,colour="black",  # Width of the error bars
                position=position_dodge(0.8)) +
  scale_y_continuous(limits = c(0,0.6)) +
  xlab("") + 
  ylab("Fraction of interneurons in MZ")+
  theme_classic() 

#figrue 1J in dotplot
gp2<- ggplot(mm2, aes(y=avgFMZ, x=Genotype, fill=Genotype))+
  facet_grid(.~Age)+
  stat_summary_bin(aes(y = avgFMZ), fun.y = "mean", geom = "bar", color = "black", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_dotplot(binaxis="y",stackdir="center",dotsize=1,  color = "black", fill = "black")+
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=14)) + 
  labs (x = 'Age/Stage', y = 'Fraction of GFP+ cells in MZ')


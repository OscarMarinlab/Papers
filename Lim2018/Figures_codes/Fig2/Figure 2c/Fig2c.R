
z1 <- read.csv('2016_zscore_all.csv', header=TRUE)
library (plotrix) 
library (ggplot2) 
library (plyr)
library(devtools)
library(plotflow)

z1$celltype <- factor(z1$celltype, levels = c("Sst.Chodl", "Sst.Cdk6", "Sst.Cbln4",  "Sst.Myh8","Sst.Tacstd2", "Sst.Th"))

mm1 <-ddply(z1, .(celltype, enriched), summarise, mmzs = mean(mean_zs), se = std.error(mean_zs))
mm1

zmz <- subset(z1, z1$enriched == "MZ")

zs1 <- ggplot(mm1, aes(x = factor(enriched), y = mmzs, fill = enriched)) + 
  facet_grid(~celltype)+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#555555","#FF6600"))+
  scale_colour_manual(values=c("#555555","#FF6600"))+
  geom_errorbar(aes(ymin=mmzs-se, ymax=mmzs+se),
                width=.2,colour="black",  # Width of the error bars
                position=position_dodge(0.8)) +
  theme_classic () +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=10)) + 
  labs (x = 'enriched', y = 'mean of zscore') 

#figure2c
gp2<- ggplot(z1, aes(x = factor(enriched), y = mean_zs, color = enriched)) + 
  facet_grid(~celltype)+
  stat_summary_bin(aes(y = mean_zs), fun.y = "mean", geom = "bar", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)+
  geom_jitter(fill = 'gray', position=position_jitter(0.2), size = 0.2)+
  scale_colour_manual(values=c("black","orange"))+
  theme_classic () +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=14)) + 
  labs (x = 'region', y = 'mean z score')




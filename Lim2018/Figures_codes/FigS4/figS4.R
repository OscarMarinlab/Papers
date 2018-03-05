library (plotrix) 
library (ggplot2) 
library (plyr)
DATA <-read.csv("E17_mafb_int.csv", header = TRUE)
DATA$Region <- factor(DATA$Region, levels = c("MZ", "CP", "SVZ"))


mm.ch1 <- DATA[DATA$Channel%in%c("tdt_int"),]
mm.ch2 <- DATA[DATA$Channel%in%c("Mafb_int"),]

#summarize mean
mm_all<- ddply(DATA, .(Channel, Region, n), summarise, mean_int= mean(Intensity), se = std.error(Intensity)) 
levels(mm_all$Channel) <- c("Mafb", "Tomato")


mma.ch1 <- mm_all[mm_all$Channel%in%c("Tomato"),]
mma.ch2 <- mm_all[mm_all$Channel%in%c("Mafb"),]
#ploting
#figSxc(right)
gp1<- ggplot(mma.ch1, aes(y=mean_int, x=Region))+
  stat_summary_bin(aes(y = mean_int), fun.y = "mean", geom = "bar", color = "black", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_dotplot(binaxis="y",stackdir="center", dotsize=3,  color = "red", fill = "red") +
  theme_classic() +
  labs (x = 'Region', y = 'Tdt_mean_intensity')
#figSxc(left)
gp2<- ggplot(mma.ch2, aes(y=mean_int, x=Region))+
  #facet_grid(Channel ~ ., scales = "free")+
  stat_summary_bin(aes(y = mean_int), fun.y = "mean", geom = "bar", color = "black", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_dotplot(binaxis="y",stackdir="center", dotsize=1,  color = "red", fill = "red")+
  #geom_jitter(width = 0.4, size = 0.1)+
  theme_classic () +
  labs (x = 'Region', y = 'Mafb_mean_intensity')

#stats
a1 <- aov(formula = mean_int ~ Region, data = mma.ch1)
stats_ch1 <- TukeyHSD(a1)
a2 <- aov(formula = mean_int ~ Region, data = mma.ch2)
stats_ch2 <- TukeyHSD(a2) # by Tukey Honest Significant Differences


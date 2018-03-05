
library(plyr)
library(ggplot2)
library (plotrix) #use for std.error

C1 <-read.csv("2016-03-28-trackcomplie.csv", header = TRUE, sep = ",")
C1$geno <- factor(C1$geno, levels = c("nkx", "dlx"))
C2<-ddply(C1, .(cell_no, region, geno), summarize, maxTime = max(T.min))
C3 <- merge(C2, C1)
C4 <- subset (C3, T.min == 0 | T.min == maxTime)
C5 <- subset (C3, T.min == maxTime)

C6a <-subset(C5, region == "soma")
C6b <-subset(C5, region == "N-tip")

tsub<- C6b[, c(1, 5, 6)]
colnames(tsub) <- c("cell_no", "Xtip", "Ytip")
C7 <- merge(C6a, tsub)

C8 <- ddply(C7, .(cell_no, geno), summarise, d_ston_X = sqrt((X.micron - Xtip)^2 + (Y.micron-Ytip)^2))
C9 <- ddply(C8, .(geno), summarise, md = mean (d_ston_X), sem = std.error(d_ston_X))
#statistic test
disttest<- t.test(d_ston_X ~ geno, data=C8)

#values for quantile
q1 <- ddply(C8, .(geno), summarise, q25 = quantile (d_ston_X, 0.25), q50 = quantile (d_ston_X, 0.5),
            q75 = quantile (d_ston_X, 0.75))


#Figure 4d - boxplot graph
gp1 <- ggplot(C8, aes(y= d_ston_X, x = geno)) +
  geom_boxplot(outlier.shape = NA, colour = "black") + 
  geom_jitter(aes(colour = geno))+
  ylab("Distance soma to time final")+
  scale_fill_manual(values=c("gray","orange"))+
  scale_colour_manual(values=c("gray","orange"))+
  theme_classic() 

#figure 4d dotplot
gp2<- ggplot(C8, aes(y= d_ston_X, x = geno)) +
  stat_summary_bin(aes(y = d_ston_X), fun.y = "mean", geom = "bar", color = "black", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  #geom_jitter(aes(colour = geno))+
  geom_dotplot(binaxis="y",stackdir="center",dotsize=0.5)+
  scale_fill_manual(values=c("dark gray","orange"))+
  scale_colour_manual(values=c("dark gray","orange"))+
  theme_classic () +
  labs (x = 'region', y = "Distance soma to time final")



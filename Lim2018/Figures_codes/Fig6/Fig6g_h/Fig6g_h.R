library (plotrix) 
library (ggplot2) 
library (plyr)
library(gridExtra)
library(grid)
data_all <-read.csv("2016-05-01_mafb_virus.csv", header = TRUE)
data_all$Genotype <- factor(data_all$Genotype, levels = c("Mafb+/fl", "Mafb-fl/fl"))
mm3 <- data_all
mm4 <-ddply(mm3, .(imagename, Brainno, Genotype), summarise, NfL1 = Length_n_L1/Length_n_t, L2to6 = Length_n_t-Length_n_L1)
mm7 <-ddply(mm4, .(Genotype), summarise, mNfL1 = mean(NfL1), se = std.error(NfL1))
mm5 <-merge(mm3, mm4)
mm6<- ddply(mm5, .(Genotype), summarise, m_L26= mean(L2to6), 
            seL26 = std.error(L2to6), mL1 = mean(Length_n_L1), seL1 = std.error(Length_n_L1))

#fig 6 box plots

#fig 6g and 6h box
fig6g <- ggplot(mm5, aes(factor(Genotype), Length_n_L1, colour = Genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 2, position=position_jitter(width=0.1))+
  scale_fill_manual(values=c("#555555","orange"))+
  scale_colour_manual(values=c("#555555","orange"))+
  theme_classic (base_size = 16) +
  xlab("Genotype") +
  ylab("total neurite length in L1") 
fig6h <- ggplot(mm5, aes(factor(Genotype), L2to6, colour = Genotype)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 2, position=position_jitter(width=0.1))+
  scale_fill_manual(values=c("#555555","orange"))+
  scale_colour_manual(values=c("#555555","orange"))+
  theme_classic (base_size = 16) +
  xlab("Genotype") +
  ylab("total neurite length in L2to6") 

#fig 6g right
ggplot(mm5, aes(x=Length_n_L1, colour=Genotype))+
  stat_ecdf(size=1, show.legend = TRUE)+
  scale_color_manual(breaks=c("Mafb+/fl","Mafb-fl/fl"),
                     values=c("#555555","#FF6600"))+
  theme_classic(base_family = "Arial", base_size = 18) +
  labs (x = 'Axon in L1', y = 'Cumulative frequency')

#fig 6h right
ggplot(mm5, aes(x=L2to6, colour=Genotype))+
  stat_ecdf(size=1, show.legend = TRUE)+
  scale_color_manual(breaks=c("Mafb+/fl","Mafb-fl/fl"),
                     values=c("#555555","#FF6600"))+
  theme_classic(base_family = "Arial", base_size = 18) +
  labs (x = 'Axon in L2to6', y = 'Cumulative frequency')

# values and statistics
q1 <- ddply(mm5, .(Genotype), summarise, q25 = quantile (Length_n_L1, 0.25), q50 = quantile (Length_n_L1, 0.5),
            q75 = quantile (Length_n_L1, 0.75))

q2 <- ddply(mm5, .(Genotype), summarise, q25 = quantile (L2to6, 0.25), q50 = quantile (L2to6, 0.5),
            q75 = quantile (L2to6, 0.75))

var.test(Length_n_L1~Genotype, data=data_all)
t4 <- t.test(L2to6~Genotype, data=mm5, var.equal = FALSE)
t2 <- t.test(Length_n_L1~Genotype, data=data_all,  var.equal = TRUE)
pt2 <-p.adjust(t2$p.value, method = "bonferroni", n = 2) #p-adjusted for t2
t3 <- t.test(Length_n_t~Genotype, data=data_all, var.equal = TRUE)



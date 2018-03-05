library (plotrix) 
library (ggplot2) 
library (plyr)
library(scales)
DATAall <-read.csv("E16_mafb_n9.csv")
DATAall$Genotype <- factor(DATAall$Genotype, levels = c("Mafb+/+::SSTcre", "Mafb-fl/fl::SSTcre"))


data1 <- DATAall[DATAall$Region%in%c("MZ", "SVZ"),]
data2 <- DATAall[DATAall$Region%in%c("All"),]
gfpcount <-ddply(data2,.(Fieldname, Genotype, Brainno), summarise, tgfp = sum(GFP_))

mm1<- ddply(data1, .(Fieldname, Genotype, Brainno, Region), summarise, GFPt= sum(GFP_)) 
mm2 <-merge (mm1, gfpcount)
data2 <- mm2[mm2$Region%in%c("MZ"),]
data2<- droplevels(data2)

data3 <- mm2[mm2$Region%in%c("SVZ"),]
data3<- droplevels(data3)

mm3<- ddply(data2, .(Fieldname, Genotype, Brainno), summarise, Mzf= (GFPt/tgfp)) 
mm4<- ddply(mm3, .(Genotype, Brainno), summarise, m_Mzf= mean(Mzf), se = std.error(Mzf))
mm5<- ddply(mm3, .(Genotype), summarise, bm_Mzf= mean(Mzf), se = std.error(Mzf)) 

mm3b<- ddply(data3, .(Fieldname, Genotype, Brainno), summarise, Svzf= (GFPt/tgfp)) 
mm4b<- ddply(mm3b, .(Genotype, Brainno), summarise, m_Svzf= mean(Svzf), se = std.error(Svzf))
mm5b<- ddply(mm3b, .(Genotype), summarise, bm_Svzf= mean(Svzf), se = std.error(Svzf)) 

t.test(m_Mzf~Genotype, data=mm4)
t.test(m_Svzf~Genotype, data=mm4b)

#boxplot
m1 <- ggplot(mm4, aes(x = factor(Genotype), y = m_Mzf, fill = Genotype)) + 
  stat_summary_bin(aes(y = m_Mzf), fun.y = "mean", geom = "bar", color = "black")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_jitter(size = 2, position=position_jitter(width=0.05), alpha =1)+
  scale_fill_manual(values=c("gray","orange"))+
  #geom_dotplot(binaxis="y",stackdir="center",dotsize=1, fill = "black")+
  #scale_colour_manual(values=c("gray","orange"))+
  ylim(0, 0.5)+
  theme_classic () 

s1 <- ggplot(mm4b, aes(x = factor(Genotype), y = m_Svzf, fill = Genotype)) + 
  stat_summary_bin(aes(y = m_Svzf), fun.y = "mean", geom = "bar", color = "black")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_jitter(size = 2, position=position_jitter(width=0.05), alpha =1)+
  scale_fill_manual(values=c("gray","orange"))+
  #geom_dotplot(binaxis="y",stackdir="center",dotsize=1, fill = "black")+
  #scale_colour_manual(values=c("gray","orange"))+
  theme_classic () 
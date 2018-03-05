library (plotrix) 
library (ggplot2) 
library (plyr)
library(scales)
#setworking directory
DATA <-read.csv("E14to16vgat.csv", header = TRUE) #reading in raw data

#data graphing order
DATA$Age <- factor(DATA$Age, levels = c("E14.5", "E15.5", "E16.5"))
#calculating density of cell per tisuee area
mmB1 <- ddply(DATA,.(imagename, Brainno, Age), summarise, d_gfp = sum(GFP)/sum(Areamm2))
#summation of total GFP cell per image
mm1<- ddply(DATA, .(imagename), summarise, GFPt= sum(GFP)) 
mm2 <-merge (mm1, DATA)#merging 2 dataframes
data2 <- mm2[mm2$Region%in%c("MZ"),] #subsetting data
mm3<- ddply(data2, .(imagename, Age, Brainno), summarise, Mzf= (GFP/GFPt)) 
mm4<- ddply(mm3, .(Age, Brainno), summarise, m_Mzf= mean(Mzf), se = std.error(Mzf))
mm5<- ddply(mm3, .(Age), summarise, bm_Mzf= mean(Mzf), se = std.error(Mzf)) 

#plotting of data

gp1 <- ggplot(mm5, aes(x = factor(Age), y = bm_Mzf)) + 
  geom_bar(stat = "identity", colour = "black", fill = "gray") + 
  geom_errorbar(aes(ymin=bm_Mzf-se, ymax=bm_Mzf+se),
                width=.2,colour="black",  # Width of the error bars
                position=position_dodge(0.8)) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=14)) + 
  labs (x = 'Age/Stage', y = 'Fraction of GFP+ cells in MZ') 

#figure 1D dotplot
gp2<- ggplot(mm4, aes(y=m_Mzf, x=Age))+
  stat_summary_bin(aes(y = m_Mzf), fun.y = "mean", geom = "bar", color = "black", fill = "white")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  geom_dotplot(binaxis="y",stackdir="center",dotsize=1)+
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=14)) + 
  labs (x = 'Age/Stage', y = 'Fraction of GFP+ cells in MZ')

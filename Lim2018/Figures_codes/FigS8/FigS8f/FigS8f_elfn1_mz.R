library (plotrix) 
library (ggplot2) 
library (plyr)
DATA1<- read.csv("elfn1_e16.csv")

DATA1$Genotype <- factor(DATA1$Genotype, levels = c("Elfn1+/-::SSTcre", "Elfn1-/-::SSTcre"))
##

data1 <- DATA1[DATA1$Region%in%c("MZ", "SVZ"),]
data2 <- DATA1[DATA1$Region%in%c("All"),]
gfpcount <-ddply(data2,.(Brainno, Genotype, Fieldname), summarise, tgfp = sum(GFP_))
data3 <-merge(data1, gfpcount)

mm1<- ddply(data3, .(Fieldname, Genotype, Region), summarise, fGFP=  GFP_/tgfp)
mm1mz <- subset(mm1, Region == "MZ")


gp <- ggplot(mm1mz, aes(x = factor(Genotype), y = fGFP), fill = Genotype) + 
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA, 
               fill = c("gray", "orange")) + 
  theme_classic () +
  ylim(0,0.8)+ ylab("Fraction of reporter+ cells")

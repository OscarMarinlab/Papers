library (plotrix) #use for std.error
library (ggplot2) #use for ggplot
library (plyr) #use for ddply
#load the data
DATA1 <- read.csv("dlx_sstpv.csv")
mm1 <- ddply(DATA1, .(imagename), summarise, tdt_t= sum(tdt), sst_tdt_t = sum(SST_tdt), pv_tdt_t = sum(PV_tdt))
mm2 <- ddply(mm1, .(imagename), summarise, sstf = sst_tdt_t/tdt_t, pvf = pv_tdt_t/tdt_t)
plot(sstf ~, data = mm2)
library(reshape)
mm3<- melt (mm2)
colnames(mm3) <- c("imagename", "celltype", "cellfr")
gp3 <- ggplot(mm3, aes(x = factor(celltype), y = cellfr), fill = celltype) + 
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA, 
               fill = c("green", "cyan")) + 
  theme_classic () +
  ylab("Fraction of reporter+ cells")

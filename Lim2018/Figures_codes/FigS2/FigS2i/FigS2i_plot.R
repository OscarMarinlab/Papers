library (plotrix) #use for std.error
library (ggplot2) #use for ggplot
library (plyr) #use for ddply
#load the data
DATA1 <- read.csv("Mzsvz_transplant.csv")
mm1 <- ddply(DATA1, .(experiment, donor), summarise, Mzf= Cell_MZ/Cell_t)

#mm2 <- ddply(mm1, .(imagename), summarise, sstf = sst_tdt_t/tdt_t, pvf = pv_tdt_t/tdt_t)

gp <- ggplot(mm1, aes(x = factor(donor), y = Mzf), fill = donor) + 
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA, 
               fill = c("gray", "orange")) + 
  theme_classic () +
  ylim(0,0.8)
  ylab("Fraction of reporter+ cells")
t.test(Mzf ~ donor, data = mm1)

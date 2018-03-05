library (plotrix) #use for std.error
library (ggplot2) #use for ggplot
library (plyr) #use for ddply
#load the data
DATA <- read.csv("2018-02-23-Mafb-Calb.csv")
DATA$Genotype <- factor(DATA$Genotype, levels = c("Mafb+/+", "Mafb+/fl", "Mafb-fl/fl"))

mm1 <- ddply(DATA, .(imagename, Layers, Genotype, Brainno, celltype), summarise, den_cell=cellcount/AreaMmsq )
mm1cal <-subset(mm1, celltype =="RCE_Calb")
mm1cal<- droplevels(mm1cal)
mm2a <- ddply(mm1, .(Layers, Genotype, celltype, Brainno), summarise, m_den = mean(den_cell), Sem = std.error(den_cell))
mm2acal<-subset(mm2a, celltype =="RCE_Calb")
mm2acal <- droplevels(mm2acal)
q1 <- ddply(mm1cal, .(Genotype, Layers), summarise, q25 = quantile (den_cell, 0.25), q50 = quantile (den_cell, 0.5),
            q75 = quantile (den_cell, 0.75))

#plot figSf
figSx <- ggplot(mm1, aes(x = Genotype, y = den_cell, colour = Genotype)) + 
  facet_grid(celltype ~ Layers, scales="free_y") +
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA) + 
  scale_fill_manual(values=c("gray","gray10","orange"))+
  scale_colour_manual(values=c("gray","gray10", "orange"))+
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=10)) + 
  ylab("density of cells")

##

#ANOVA analysis for double lable GFP+Calbindin
stat3 <- aov(m_den ~Genotype*Layers, data = mm2acal)
stat3t <-TukeyHSD(stat3)
s3<- stat3t$`Genotype:Layers`
#posthoc for each layer
subset(s3, grepl("1",row.names(s3)) & !grepl("[2-6]", row.names (s3)))
subset(s3, grepl("2",row.names(s3)) & !grepl("[1, 4-6]", row.names (s3)))
subset(s3, grepl("4",row.names(s3)) & !grepl("[1-3,5,6]", row.names (s3)))
subset(s3, grepl("5",row.names(s3)) & !grepl("[1-4,6]", row.names (s3)))
subset(s3, grepl("6",row.names(s3)) & !grepl("[1-5]", row.names (s3)))


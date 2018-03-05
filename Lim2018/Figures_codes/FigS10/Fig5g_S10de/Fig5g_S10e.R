library (plotrix) #use for std.error
library (ggplot2) #use for ggplot
library (plyr) #use for ddply
#load the data
DATA <- read.csv("Mafb_P21_n3_all_CR.csv")
DATA$Genotype <- factor(DATA$Genotype, levels = c("Mafb+/+", "Mafb-fl/+", "Mafb-fl/fl"))

mm1 <- ddply(DATA, .(imagename, Layers, Genotype, Brainno, celltype), summarise, den_cell=cellcount/AreaMmsq )
mm1GFP <-subset(mm1, celltype =="RCE" & Genotype !="Mafb-fl/+")
mm1GFP<- droplevels(mm1GFP)
mm2a <- ddply(mm1, .(Layers, Genotype, celltype, Brainno), summarise, m_den = mean(den_cell), Sem = std.error(den_cell))
mm2aGFP<-subset(mm2a, celltype =="RCE" & Genotype !="Mafb-fl/+")
mm2aGFP <- droplevels(mm2aGFP)
mm2c<-ddply(mm1, .(Layers, Genotype, celltype), summarise, m_den = mean(den_cell), Sem = std.error(den_cell))
mm2b<-ddply(mm2a, .(Layers, Genotype, celltype),  summarise, mb_den = mean(m_den), Sem = sd(m_den))
mm1b<-subset(mm1, celltype == "RCE_CR")
mm1c<-subset(mm1, celltype == "RCE")

q1 <- ddply(mm1GFP, .(Genotype, Layers), summarise, q25 = quantile (den_cell, 0.25), q50 = quantile (den_cell, 0.5),
            q75 = quantile (den_cell, 0.75))

q2 <- ddply(mm1b, .(Genotype, Layers), summarise, q25 = quantile (den_cell, 0.25), q50 = quantile (den_cell, 0.5),
            q75 = quantile (den_cell, 0.75))

q3 <- ddply(mm1c, .(Genotype, Layers), summarise, q25 = quantile (den_cell, 0.25), q50 = quantile (den_cell, 0.5),
            q75 = quantile (den_cell, 0.75))


#plot GFP/ SST figure for WT and fl/fl - figure 5g
fig5g <- ggplot(mm1GFP, aes(y= den_cell, x = Genotype, colour = Genotype))+
  facet_grid(.~Layers ) +
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA) + 
  scale_fill_manual(values=c("gray","orange"))+
  scale_colour_manual(values=c("gray","orange"))+
  ylab("density")+
  theme_classic()

##

#FigS9

#FigS9d & e
gp2 <- ggplot(mm1, aes(x = Genotype, y = den_cell, colour = Genotype)) + 
  facet_grid(celltype ~ Layers, scales="free_y") +
  geom_boxplot(outlier.colour  = NA, linetype = "dotted") +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.colour  = NA) + 
  scale_fill_manual(values=c("gray","gray10","orange"))+
  scale_colour_manual(values=c("gray","gray10", "orange"))+
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=10)) + 
  ylab("density of cells")

#ANOVA analysis for RCE
GFP <- subset(mm2a, celltype == "RCE")
stat3 <- aov(m_den ~Genotype*Layers, data = GFP)
stat3t <-TukeyHSD(stat3)
s3<- stat3t$`Genotype:Layers`
#posthoc for each layer
subset(s3, grepl("1",row.names(s3)) & !grepl("[2-6]", row.names (s3)))
subset(s3, grepl("2",row.names(s3)) & !grepl("[1, 4-6]", row.names (s3)))
subset(s3, grepl("4",row.names(s3)) & !grepl("[1-3,5,6]", row.names (s3)))
subset(s3, grepl("5",row.names(s3)) & !grepl("[1-4,6]", row.names (s3)))
subset(s3, grepl("6",row.names(s3)) & !grepl("[1-5]", row.names (s3)))

#ANOVA analysis for RCE_CR
GFP_cr <- subset(mm2a, celltype == "RCE_CR")
stat2 <- aov(m_den ~Genotype*Layers, data = GFP_cr)
stat2t <-TukeyHSD(stat2)
s2<- stat2t$`Genotype:Layers`
#posthoc for each layer
subset(s2, grepl("1",row.names(s2)) & !grepl("[2-6]", row.names (s2)))
subset(s2, grepl("2",row.names(s2)) & !grepl("[1, 4-6]", row.names (s2)))
subset(s2, grepl("4",row.names(s2)) & !grepl("[1-3,5,6]", row.names (s2)))
subset(s2, grepl("5",row.names(s2)) & !grepl("[1-4,6]", row.names (s2)))
subset(s2, grepl("6",row.names(s2)) & !grepl("[1-5]", row.names (s2)))
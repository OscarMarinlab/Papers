
library (plotrix) #use for std.error
library (ggplot2) #use for ggplot
library (plyr) #use for ddply
Het <- read.csv('Het15psCRACM_avg.csv', header = TRUE)
HO <- read.csv('HO15psCRACM_avg.csv', header = TRUE)
data <-rbind(Het, HO)
celldepth <- read.csv('cellid_depth.csv')
data <-merge(data, celldepth)
data_L2 <- subset(data,Layer=='L2')
data_L3 <-subset (data, Layer == 'L3')
data_L3_r1 = subset(data_L3, Row == 1)
data_L3$Row <- factor (data_L3$Row)

#subsetting data to check n for each group
data_L3_r1_het = subset(data_L3_r1, Genotype == "Het")
data_L3_r1_het = droplevels(data_L3_r1_het)

data_L3_r1_HO = subset(data_L3_r1, Genotype == "HO")
data_L3_r1_HO = droplevels(data_L3_r1_HO)

##stats test
t.test(depth ~ Genotype, data_L3_r1) #check patching depth is the same for both groups
a1 <- aov(MeanIPSC  ~ Genotype * Row, data_L3) #ANOVA for MeanIPSC
summary (a1)
TukeyHSD(a1)

#getting values for table
q1 <- ddply(data_L3, .(Genotype, Row), summarise, q25 = quantile (MeanIPSC, 0.25), q50 = quantile (MeanIPSC, 0.5),
            q75 = quantile (MeanIPSC, 0.75))

#box plot fig6i 
p3b <-ggplot(data_L3, aes(y=MeanIPSC, x=Genotype, color = Genotype)) +
  facet_grid(.~ Row, scales="free", space="free") +
  geom_boxplot(aes(colour = Genotype), outlier.shape = NA) +
  geom_jitter(position=position_jitter(w=0.1, h=0.1)) +
  scale_fill_manual(values=c("gray45","dark orange"))+
  scale_color_manual(values=c("gray45","dark orange")) + 
  theme_classic (base_size = 14) +
  ggtitle("L3 pyramids (IPSC) in each row all data") 


#heatmap fig6k
#heatmaping the data
data_L3_het = subset(data_L3, Genotype == "Het")
data_L3_het = droplevels(data_L3_het)

data_L3_HO = subset(data_L3, Genotype == "HO")
data_L3_HO = droplevels(data_L3_HO)
data_L3_het <- data_L3_het[order(data_L3_het$depth),]
data_L3_HO <- data_L3_HO[order(data_L3_HO$depth),]
library(reshape2)
library(gplots)
mhet<- acast(data_L3_het, Row~depth + cellid, value.var="MeanIPSC")
mhet_matrix <- data.matrix(mhet)

mKO<- acast(data_L3_HO, Row~depth + cellid, value.var="MeanIPSC")
mKO_matrix <- data.matrix(mKO)
tmatrix <-cbind(mhet_matrix, mKO_matrix)
colfunc1 <- colorRampPalette(c("white", "yellow", "orange", "red", "black" ))
hmcols <- c(colfunc1(500))

#fullheatmap of het then ho
heatmap.2(tmatrix, key = FALSE ,dendrogram="none",trace="none",scale="none", col=hmcols,
          Rowv=FALSE, Colv=FALSE)



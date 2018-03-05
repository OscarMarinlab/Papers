library(ggplot2)

M_allen <- read.csv("cell_data_markers.csv")

Mz_allen_sub <- M_allen[complete.cases(M_allen),]

row.names(Mz_allen_sub) <- Mz_allen_sub$symbol
ma3<-Mz_allen_sub[, 3:204]

m_matrix <- data.matrix(ma3)
ms <-scale (t(m_matrix))
ms2 <-t(ms)

#df[grep("^Andy",rownames(df)),]

##
x1 <- grep("Sst.Tacstd2*",  colnames(ms2))
x2 <- grep("Sst.Chodl*", colnames(ms2))
x3 <- grep("Sst.Cbln4*",  colnames(ms2))
x4 <- grep("Sst.Cdk6*",  colnames(ms2))
x5 <- grep("Sst.Myh8*",  colnames(ms2))
x7 <- grep("Sst.Th*",  colnames(ms2))
colnames(ms2)[x1] <- c("Sst.Tacstd2")
colnames(ms2)[x2] <- c("Sst.Chodl")
colnames(ms2)[x3] <- c("Sst.Cbln4")
colnames(ms2)[x4] <- c("Sst.Cdk6")
colnames(ms2)[x5] <- c("Sst.Myh8")
colnames(ms2)[x7] <- c("Sst.Th")


rm <-rownames(ms2)


Sst.Tacstd2 <- rowMeans(ms2[, x1], na.rm = FALSE, dims = 1)
Sst.Chodl<- rowMeans(ms2[, x2], na.rm = FALSE, dims = 1)
Sst.Cbln4 <- rowMeans(ms2[, x3], na.rm = FALSE, dims = 1)
Sst.Cdk6 <- rowMeans(ms2[, x4], na.rm = FALSE, dims = 1)
Sst.Myh8<-rowMeans(ms2[, x5], na.rm = FALSE, dims = 1)
Sst.Th<-rowMeans(ms2[, x7], na.rm = FALSE, dims = 1)

All_zmean <- cbind(Sst.Chodl, Sst.Cdk6, Sst.Cbln4,  Sst.Myh8, Sst.Tacstd2,Sst.Th)

All_zmean <- as.data.frame(as.table(All_zmean))
colnames(All_zmean) <-c("genes", "celltype", "mean_zs")

library (plotrix) 
library (ggplot2) 
library (plyr)
mm1 <-ddply(All_zmean, .(celltype), summarise, mmzs = mean(mean_zs), se = std.error(mean_zs))
gp1 <- ggplot(mm1, aes(x = factor(celltype), y = mmzs)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=mmzs-se, ymax=mmzs+se),
                width=.2,colour="black",  # Width of the error bars
                position=position_dodge(0.8)) +
  theme_classic (base_size = 18) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=14)) + 
  labs (x = 'Celltype', y = 'mean of MZ enriched zscore') 

gp4b <- ggplot(All_zmean, aes(factor(x = genes), y= mean_zs)) + 
  facet_grid (~ celltype) +
  #geom_boxplot(outlier.size=NA) +  
  geom_bar(stat = "identity") +
  theme_classic (base_size = 16) +
  theme(axis.text.x= element_text(angle=90, vjust=0.5, size=6)) + 
  xlab("genes in MZ") +
  ylab("mean_zscore")

color.map <- function(colnames) 
  { if (colnames=="Sst.Chodl") "bisque4"
    else if (colnames=="Sst.Cbln4") "orange"
    else if (colnames=="Sst.Cdk6")  "chocolate4"
    else if (colnames == "Sst.Tacstd2") "coral"
    else if (colnames=="Sst.Myh8") "gold"
    else if (colnames(ms2)[x2]) "tomato"}


color.map <- function(colnames) 
{ if (colnames=="Sst.Chodl") "bisque4"
  else if (colnames=="Sst.Cbln4") "orange"
  else if (colnames=="Sst.Cdk6")  "chocolate4"
  else if (colnames=="Sst.Th") "tomato"
  else if (colnames == "Sst.Tacstd2") "coral"
  else if (colnames=="Sst.Myh8") "gold"
  }
cellcolors <- unlist(lapply(colnames(ms2), color.map))

##
library(RColorBrewer)
library (ggplot2)
library(gplots)
colfunc <- colorRampPalette(c("black", "black", "white", "orange", "red"))
distancem <- dist(ms2)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
hmcol<- hmcol<-brewer.pal(11,"RdBu")
hm<- heatmap.2(ms2, dendrogram="row", Colv = FALSE, col=colfunc(300),  scale="row", trace="none",  key = TRUE,
               density.info ="none",
               ColSideColors=cellcolors)
dend <- as.hclust( hm$rowDendrogram )
library("scales")



# creates a 3 x 10 inch image
png("heatmap_r3.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 3*400,
    res = 300,   # 300 pixels per inch
    #compression = "none",
    bg= "transparent",
    pointsize = 10)        # smaller font size


heatmap.2(ms2,
          cellnote = ms2,  # same data set for cell labels
          col=colfunc(25),
          scale="row"
          Rowv=dendcompletem,
          key = TRUE,
          #main = "MZ vs SVZ", # heat map title
          #notecol="black",      # change font color of cell labels to black
          #density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(12,9),     # widens margins around plot
          col=hmcol,            # use on color palette defined earlier 
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering


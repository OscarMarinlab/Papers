library(ggplot2)

Mz_allen <- read.csv("MZ_allencelltype.csv")
mz.sub<-read.csv('mzsubfdr5.csv')
Mz_allen_sub <- Mz_allen[complete.cases(Mz_allen),]

ma <- merge(mz.sub, Mz_allen_sub, by = "symbol")
n_allen <- colnames(Mz_allen_sub)
ma2 <- ma[, n_allen]
row.names(ma2) <- ma2$symbol
ma3<-ma2[, 3:204]

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
x6 <- x7 [15:33]
colnames(ms2)[x1] <- c("Sst.Tacstd2")
colnames(ms2)[x2] <- c("Sst.Chodl")
colnames(ms2)[x3] <- c("Sst.Cbln4")
colnames(ms2)[x4] <- c("Sst.Cdk6")
colnames(ms2)[x5] <- c("Sst.Myh8")
colnames(ms2)[x6] <- c("Sst.Th")

ms4<-melt (ms2)

Sst.Tacstd2 <- rowMeans(ms2[, x1], na.rm = FALSE, dims = 1)
Sst.Chodl<- rowMeans(ms2[, x2], na.rm = FALSE, dims = 1)
Sst.Cbln4 <- rowMeans(ms2[, x3], na.rm = FALSE, dims = 1)
Sst.Cdk6 <- rowMeans(ms2[, x4], na.rm = FALSE, dims = 1)
Sst.Myh8<-rowMeans(ms2[, x5], na.rm = FALSE, dims = 1)
Sst.Th<-rowMeans(ms2[, x6], na.rm = FALSE, dims = 1)

All_zmean <- cbind(Sst.Chodl, Sst.Cdk6, Sst.Cbln4, Sst.Myh8, Sst.Tacstd2, Sst.Th)

All_zmean <- as.data.frame(as.table(All_zmean))
colnames(All_zmean) <-c("genes", "celltype", "mean_zs")
rn<- rownames(ms3)

#figure 2e-violin
ggplot(ms4, aes(x = factor(Var1), y = value, color = Var1)) + 
  facet_grid(~Var2)+
  geom_violin(scale = "width")+
  geom_jitter(position=position_jitter(0.03), size = 0.03)+
  theme_classic () +
  theme(legend.position="none")+
  theme(axis.text.x= element_text(angle=90, vjust=0.5)) + 
  labs (x = 'genes', y = ' z score') +
  scale_colour_manual(values=c("gray","gray", "gray", "gray",
                               "gray", "gray", "green3", "gray", "gray",
                               "brown", "gray", "gray", "gray", "dodgerblue1", "red", "yellow",
                               "gray","gray", "gray", "gray", "gray", "gray",
                               "gray","gray"))
 


color.map <- function(colnames) 
  { if (colnames=="Sst.Chodl") "bisque4"
    else if (colnames=="Sst.Cbln4") "orange"
    else if (colnames=="Sst.Cdk6")  "chocolate4"
    else if (colnames == "Sst.Tacstd2") "coral"
    else if (colnames=="Sst.Myh8") "gold"
    else if (colnames=="Sst.Th") "tomato"}
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

#Figure 2b, heatmap
heatmap.2(ms2,
          cellnote = ms2,  # same data set for cell labels
          col=colfunc(25),
          scale="row",
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



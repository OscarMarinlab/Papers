library (ggplot2)
library(gplots)
fdr5 <-read.csv("count_fdr5.csv", header = TRUE, sep = ",")
fdr5.sub <- subset(fdr5, logCPM > 5)
sfdr5 <- fdr5.sub[order(fdr5.sub$logFC, fdr5.sub$FDR) ,]
mz.sub <-subset(sfdr5, logFC<0)
svz.fdr5 <-subset(sfdr5, logFC>0)
row.names(m_fdr5) <- fdr5.sub$SYMBOL
m_fdr5<-m_fdr5[, 8:13]
m_matrix <- data.matrix(m_fdr5)
ms <-scale (t(m_matrix))
ms2 <-t(ms)
library(RColorBrewer)
distancem <- dist(ms2)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
hmcol<-brewer.pal(9,"YlOrRd")
heatmap.2(ms2, Rowv=dendcompletem, col=hmcol, Colv=NA, scale="none", trace="none",  key = TRUE,
          keysize = 1)

# creates a 3 x 10 inch image
png("heatmap_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,   # 300 pixels per inch
    compression = "none", 
    pointsize = 10)        # smaller font size

heatmap.2(ms2,
          cellnote = ms2,  # same data set for cell labels
          main = "MZ vs SVZ", # heat map title
          #notecol="black",      # change font color of cell labels to black
          #density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=hmcol,            # use on color palette defined earlier 
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering



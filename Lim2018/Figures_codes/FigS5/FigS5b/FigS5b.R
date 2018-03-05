svz_allen <- read.csv("allen_cell_data_svz.csv")
svz_allen_sub <- svz_allen[complete.cases(svz_allen),]
svz.fdr5 <- read.csv("svzfdr5.csv")
names(svz.fdr5)[names(svz.fdr5)=="SYMBOL"] <- "symbol"
sa <- merge(svz.fdr5, svz_allen_sub, by = "symbol")
svz_allen <- colnames(svz_allen_sub)
sa2 <- sa[, svz_allen]
row.names(sa2) <- sa2$symbol
sa3<-sa2[, 3:204]


s_matrix <- data.matrix(sa3)
ss <-scale (t(s_matrix))
ss2 <-t(ss)

#df[grep("^Andy",rownames(df)),]

x1 <- grep("Sst.Tacstd2*",  colnames(ss2))
x2 <- grep("Sst.Chodl*", colnames(ss2))
x3 <- grep("Sst.Cbln4*",  colnames(ss2))
x4 <- grep("Sst.Cdk6*",  colnames(ss2))
x5 <- grep("Sst.Myh8*",  colnames(ss2))
x7 <- grep("Sst.Th*",  colnames(ss2))
colnames(ms2)[x1] <- c("Sst.Tacstd2")
colnames(ms2)[x2] <- c("Sst.Chodl")
colnames(ms2)[x3] <- c("Sst.Cbln4")
colnames(ms2)[x4] <- c("Sst.Cdk6")
colnames(ms2)[x5] <- c("Sst.Myh8")
colnames(ms2)[x7] <- c("Sst.Th")
msz <-rbind(ms2, ss2)
msz2 <- msz


Sst.Tacstd2 <- rowMeans(ss2[, x1], na.rm = FALSE, dims = 1)
Sst.Chodl<- rowMeans(ss2[, x2], na.rm = FALSE, dims = 1)
Sst.Cbln4 <- rowMeans(ss2[, x3], na.rm = FALSE, dims = 1)
Sst.Cdk6 <- rowMeans(ss2[, x4], na.rm = FALSE, dims = 1)
Sst.Myh8<-rowMeans(ss2[, x5], na.rm = FALSE, dims = 1)
Sst.Th<-rowMeans(ss2[, x7], na.rm = FALSE, dims = 1)

color.map <- function(colnames) 
{ if (colnames=="Sst.Chodl") "bisque4"
  else if (colnames=="Sst.Cbln4") "orange"
  else if (colnames=="Sst.Cdk6")  "chocolate4"
  else if (colnames == "Sst.Tacstd2") "coral"
  else if (colnames=="Sst.Myh8") "gold"
  else if (colnames=="Sst.Th") "tomato"}
cellcolors <- unlist(lapply(colnames(ms2), color.map))
rc[1:24] <-c('gray')
rc[25:79] <-c('red')

color.map2 <- function(rownames) 
{ if (rownames=="MZenriched") "green"
  else "blue"
}
rowcolors <- unlist(lapply(colnames(msz2), color.map2))
##
library(RColorBrewer)
library (ggplot2)
library(gplots)
colfunc <- colorRampPalette(c("black", "black", "white", "orange", "red"))
sdistancem <- dist(ss2)
shclust_completem <- hclust(sdistancem, method = "complete")
sdendcompletem <- as.dendrogram(shclust_completem)
hmcol<- hmcol<-brewer.pal(11,"RdBu")
heatmap.2(ss2, dendrogram="row", Colv = FALSE, col=colfunc(300),  scale="row", trace="none",  key = TRUE,
               density.info ="none",
               ColSideColors=cellcolors)
dend <- as.hclust( hm$rowDendrogram )

##allgenes

adistancem <- dist(msz)
ahclust_completem <- hclust(adistancem, method = "complete")
adendcompletem <- as.dendrogram(ahclust_completem)
heatmap.2(msz, dendrogram="row", Colv = FALSE, col=colfunc(300),  scale="row", trace="none",  key = TRUE,
          density.info ="none",
          ColSideColors=cellcolors,
          RowSideColors = rc)
dend <- as.hclust( hm$rowDendrogram )
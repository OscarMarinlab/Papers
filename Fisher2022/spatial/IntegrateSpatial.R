#Integrate:
#This script takes all the processed samples and merges them together. We check the name of the file for tags indciating well
#We then filter and cluster the data


## ----receive input from nextflow----------------------------------------------

#Get parameters
args <- commandArgs(TRUE)
path<-args[1]
timepoint <- names(read.delim(args[2]))


print(paste0(args[4],' files passed to integration'))
nsamples<-as.numeric(args[4])
print(nsamples)
set.seed(123)
print(timepoint)

## ----libraries----------------------------

#Load standard package list
source(paste0(path,'/bin/auxiliary/DirectoryChecker.R'))
source(paste0(path,'/bin/auxiliary/PackageLoader.R'))
source(paste0(path,'/bin/auxiliary/SetPlottingParameters.R'))


#Get objects 
samples<-args[5:(4+nsamples)]
samples<-gtools::mixedsort(samples)
samples <- samples[grep(paste0('^', timepoint), samples)]



ncor=detectCores()
print(paste(ncor,'cores detected'))
plan("multiprocess", workers = ncor -1 )
options(future.globals.maxSize = 15000 * 1024^2)


#collect prefixes to integrate timepoints separately
# timepoints<-unlist(str_split(samples[1],'_'))[1]
# timepoint <- unique(timepoints)



#load canoncial markers
canonmarkers<-read.csv(paste0(path,'/input/supplementary/CanonicalMarkers.csv'))
canonmarkers[,1] <- toupper(canonmarkers[,1])
markers<-list()
markers[['Martinotti']]<-canonmarkers[grep('^Martinotti',canonmarkers[,2]),1]
markers[['Non-Martinotti']]<-canonmarkers[grep('Non-Martinotti',canonmarkers[,2]),1]
markers[['LRP']]<-canonmarkers[grep('projecting',canonmarkers[,2]),1]
markers[['Stressed']]<-canonmarkers[grep('Stress',canonmarkers[,2]),1]

#load markers
load(args[3])


# ---- files ----------------------------

print(samples)
samples <- mixedsort(samples)
objectList<-vector()
for(i in 1:length(samples)){
  x<-load(samples[i])
  name<- paste0('obj',as.character(i))
  assign(name,get(x[1]))
  objectList<-c(objectList,get(x))
  rm(x)
}


# objectList<-lapply(objectList, FUN=function(obj){
#   obj<-subset(obj, cells=WhichCells(obj,expression=Region!='CPu' & GeneralRegion != '' & GeneralLayers != ''))
# })


# ---- Integrate --------------------------------------------------

anch<-FindIntegrationAnchors(objectList)

integrated<-IntegrateData(anch)

integratedFull <- integrated

integratedFull$GeneralLayers<-as.character(integratedFull$GeneralLayers)


# integratedFull <- subset(integratedFull,cells=WhichCells(integratedFull,expression=Region!='CPu'))
# 
# #check if region/layers all blank
# if(sum(integratedFull$GeneralLayers != '') > 0){
#   integratedFull <- subset(integratedFull,cells=WhichCells(integratedFull,expression=GeneralLayers != ''))
# }
# if(sum(integratedFull$GeneralRegion != '') > 0){
#   integratedFull <- subset(integratedFull,cells=WhichCells(integratedFull,expression=GeneralRegion != ''))
# }

#ALTERNATIVE FILTERING


#take integrated full, keep SST +ve cells, clean PAX6 +ve

plotList<-list()
props<-c()
keepCells<-colnames(integratedFull)
DefaultAssay(integratedFull)<-'RNA'
df<-data.frame(FetchData(integratedFull,slot='data', c('SST','GAD1','GAD2','LHX6','PAX6')))


for(gene in c('SST','PAX6','GAD1','GAD2','LHX6')){
  
  dat<-df[,gene]
  d<-density(dat)
  
  #approximate peaks
  peaks<-d$x[which(diff(sign(diff(d$y)))<0)+1]
  
  valleys <- d$x[which(diff(sign(diff(d$y)))>0)+1]
  
  #ignore peaks too close to th emax/min values
  peaks<-peaks[which(abs(peaks - max(dat)) > 0.1)]
  peaks<-peaks[which(abs(peaks - min(dat)) > 0.1)]
  
  #thresh<-max(peaks)
  thresh <- min(valleys)
  col=smoothBlue[3]
  if(gene == 'PAX6'){
    thresh <- min(valleys)
    col=smoothOrange[2]
    #filter cells
    keepCells<-keepCells [keepCells %in% rownames(df)[which(dat<thresh)]]
    prop<- sum(dat>thresh) / length(dat)
    props<-c(props,prop)
  }else if(gene == 'SST'){
    thresh <- min(valleys)
    #filter cells
    keepCells<-keepCells [keepCells %in% rownames(df)[which(dat>thresh)]]
    prop<- sum(dat>thresh) / length(dat)
    props<-c(props,prop)
    
    }
  else{
    #filter cells
    keepCells<-keepCells [keepCells %in% rownames(df)[which(dat>thresh)]]
    prop<- sum(dat>thresh) / length(dat)
    props<-c(props,prop)
  }
  
  
  p <- ggplot(df, aes_string(x=gene)) + 
    geom_density(data=df,color=col,
                 alpha=.3, fill=col) +
    geom_vline(data=df, 
               aes(xintercept=thresh), 
               colour=col, linetype="dashed", size=1)
  
  dpb <- ggplot_build(p)
  
  x1 <- min(which(dpb$data[[1]]$x >= -5))
  x2 <- max(which(dpb$data[[1]]$x <= thresh))
  
  p<- p +
    geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                              y=dpb$data[[1]]$y[x1:x2]),
              aes(x=x, y=y), fill="light grey") + 
    theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg)
  
  plotList[[gene]]<-plot_grid(p, ncol=1)
  
  
}


spatialPlots <- plotList

plot_grid(plotlist = plotList)
length(keepCells)



integratedFullSST <- subset(integratedFull, cells= keepCells)


DefaultAssay(integratedFullSST) <- 'integrated'
integratedFullSST<-FindVariableFeatures(integratedFullSST)
integratedFullSST<-ScaleData(integratedFullSST, features = rownames(integratedFullSST))
integratedFullSST<-RunPCA(integratedFullSST)
integratedFullSST<-FindNeighbors(integratedFullSST)
integratedFullSST<-FindClusters(integratedFullSST, resolution = 1)
integratedFullSST<-RunUMAP(integratedFullSST,dims=1:15)

save(integratedFull,integratedFullSST, file=paste0(path,'/output/spatial/',timepoint,'/integrated_',timepoint,'_full.Rdata'),compress = TRUE)


DimPlot(integratedFullSST, label=T)
DefaultAssay(integratedFullSST) <- 'RNA'
FeaturePlot(integratedFullSST, order=T,  c('SST','MEIS2','GAD1','GAD2','LHX6','PAX6','DACH1','CHODL','ERBB4','TAC1','NFIB')) & scale_color_gradientn(colors=FeatureCol)

Idents(integratedFullSST) <- 'seurat_clusters'
DefaultAssay(integratedFullSST) <- 'RNA'
DEGenes <- FindAllMarkers(integratedFullSST, only.pos=T, logfc.threshold = log(1.2), base=exp(1))

# 
# #take clusters with interneuron markers differentially expressed
GAD1clusters <- as.character(unique(DEGenes$cluster[grep('GAD1', DEGenes$gene)]))
GAD2clusters <- as.character(unique(DEGenes$cluster[grep('GAD2', DEGenes$gene)]))
LHX6clusters <- as.character(unique(DEGenes$cluster[grep('LHX6', DEGenes$gene)]))
# 
# clusters <- intersect(GAD1clusters, intersect(GAD2clusters, LHX6clusters))
# # 
# clusters <- unique(DEGenes$cluster[grep('SST', DEGenes$gene)])
clusters <- unique(c(GAD1clusters,GAD2clusters,LHX6clusters))


if(timepoint != 'E16'){
  integratedFullSST <- subset(integratedFullSST,cells=WhichCells(integratedFullSST,expression=Region!='CPu'))
  
  #check if region/layers all blank
  if(sum(integratedFullSST$GeneralLayers != '') > 0){
    integratedFullSST <- subset(integratedFullSST,cells=WhichCells(integratedFullSST,expression=GeneralLayers != ''))
  }
  if(sum(integratedFullSST$GeneralRegion != '') > 0){
    integratedFullSST <- subset(integratedFullSST,cells=WhichCells(integratedFullSST,expression=GeneralRegion != ''))
  }
  
}





if(timepoint != 'E16'){
  
  keepCells <- WhichCells(integratedFullSST, expression=seurat_clusters %in% clusters)
}else{
  keepCells <- colnames(integratedFullSST)
}

length(keepCells)



#---- record QC plots --------------------------------------------------
# pdf(file=paste0(path,'/output/spatial/',timepoint,'/filtering_plots.pdf'))
# plot_grid(plotlist = plotList)
# dev.off()


#---- filter --------------------------------------------------
integratedSST_clean<-subset(integratedFullSST,cells=keepCells)
#integratedSST_clean<-subset(integratedSST_clean, cells=WhichCells(integratedSST_clean,expression=Region!='CPu' & GeneralRegion != '' & GeneralLayers != ''))
#lostCells<-colnames(integratedSST)[colnames(integratedSST) %nin% colnames(integratedSST_clean)]




#---- correct cell names on image (for plotting reasons) --------------------------------------------------
if(length(integratedSST_clean@images) > 0){
  
  for(id in unique(integratedSST_clean$orig.ident)){
    
    cells<-names(integratedSST_clean$orig.ident[ integratedSST_clean$orig.ident == id])
    val<-unlist(lapply(str_split(cells,'_'), FUN = function(x) x[2]))[1]
    rownames(integratedSST_clean@images[[id]]@coordinates) <- paste0(rownames(integratedSST_clean@images[[id]]@coordinates), '_' , val)
    
  }
  
}


#---- process   --------------------------------------------------
DefaultAssay(integratedSST_clean) <- 'integrated'
integratedSST_clean<-FindVariableFeatures(integratedSST_clean)
integratedSST_clean<-ScaleData(integratedSST_clean, features = rownames(integratedSST_clean), assay = 'integrated')
integratedSST_clean<-RunPCA(integratedSST_clean)
integratedSST_clean<-FindNeighbors(integratedSST_clean)
integratedSST_clean<-FindClusters(integratedSST_clean, resolution = 1)
integratedSST_clean<-RunUMAP(integratedSST_clean,dims=1:20)
DefaultAssay(integratedSST_clean) <- 'RNA'


integratedSST <- integratedSST_clean


integratedSST$seurat_clusters <- as.numeric(integratedSST$seurat_clusters) 



#reannotate with gene labels
#reannotate??
Idents(integratedSST) <- 'seurat_clusters'
DefaultAssay(integratedSST)
DEGenes<- FindAllMarkers(integratedSST, only.pos = T, logfc.threshold = log(1.5), base=exp(1))


#reannotation
integratedSST$reannotation <- NA
integratedSST$putative_identity <- NA

#Quantify proportion of markers that appear in DE genes.
for(cluster in unique(Idents(integratedSST))){
  
  sub<-DEGenes[DEGenes$cluster==cluster,]
  scores<-vector('list',length=3)
  names(scores)<- c('Martinotti','Non-Martinotti','LRP')
  
  #Compute DEscore for each type
  for(type in c('Martinotti','Non-Martinotti','LRP')){
    mark<-typeMarkers[[type]]
    m<- sub[which(sub$gene %in% mark),]
    spec<-m$pct.1-m$pct.2
    names(spec)<-m$gene
    DEscore<- -log10(m$p_val_adj)
    DEscore<-unlist(mapply(DEscore,FUN=function(x) min(x,20)))
    #DEscore<-unlist(mapply(DEscore,FUN=function(x) min(x,20)))*spec
    
    if(dim(m)[1] == 0){
      scores[[type]]<-c(sum(DEscore),0)
    }else{
      scores[[type]]<-c(sum(DEscore),m$gene[grep(max(m$avg_logFC),m$avg_logFC)])
    }
    
    
  }
  
  vals<-as.numeric(unlist(lapply(scores,FUN=function(x)x[[1]])))
  names(vals)<-names(scores)
  if(length(grep(max(vals),vals))>1){
    
    #take the gene with highest logfc
    gene<-sub$gene[grep(max(sub$avg_logFC), sub$avg_logFC)]
    cellType<-NA
    
  }else{
    
    
    #which cell type?
    cellType<-names(vals)[grep(max(vals), vals)]
    canon<-markers[[cellType]]
    full<-typeMarkers[[cellType]]
    
    #Find gene identity
    #is there a canonical one?
    if(sum(sub$gene %in% canon) >= 1){
      
      marks<-sub[sub$gene %in% canon,]
      gene<-marks$gene[grep(max(marks$avg_logFC), marks$avg_logFC)]
      
    }else if(sum(sub$gene %in% full) >= 1){
      
      marks<-sub[sub$gene %in% full,]
      gene<-marks$gene[grep(max(marks$avg_logFC), marks$avg_logFC)]
      
    }else{
      #take the gene with highest logfc
      gene<-sub$gene[grep(max(sub$avg_logFC), sub$avg_logFC)]
      
    }
  }
  
  
  #has that identity already been found in another cluster?
  if(length(grep(gene,integratedSST$reannotation))>1){
    #get all the matching labels, seek maximum number on the end
    nams<-integratedSST$reannotation[grep(gene,integratedSST$reannotation)]
    nums<-lapply(nams,FUN=function(x) substring(x,nchar(x)))
    m<-max(as.numeric(nums))
    m<-m+1
    integratedSST$reannotation[which(integratedSST$seurat_clusters==cluster)]<-paste0('Sst_',gene,'_',as.character(m))
    
  }
  else{
    integratedSST$reannotation[which(integratedSST$seurat_clusters==cluster)]<-paste0('Sst_',gene,'_1')
  }
  
  integratedSST$putative_identity[which(integratedSST$seurat_clusters==cluster)]<-cellType
  
}

DimPlot(integratedSST, group.by='reannotation', label=T)
DefaultAssay(integratedSST) <- 'RNA'
FeaturePlot(integratedSST, order=T, c('SST','MEIS2','GAD1','GAD2','LHX6','RELN','CHODL','ERBB4','TAC1','PAX6','DACH1')) & scale_color_gradientn(colors=FeatureCol)


#---- save --------------------------------------------------

integratedSST$timepoint <- timepoint

#Save to pipeline directory for later reference
save(integratedSST,integratedFullSST, integratedFull, spatialPlots, file=paste0(path,'/output/spatial/',timepoint,'/integrated_',timepoint,'_SST.Rdata'),compress = TRUE)
#Pass to working directory
save(integratedSST,integratedFullSST, integratedFull, spatialPlots, file=paste0('integrated_',timepoint,'.Rdata'),compress = TRUE)










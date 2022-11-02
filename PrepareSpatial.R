
#PrepareSpatial:
#This script loads each well of the spatial samples, alongside the metadata giving spatial position. 
#It compiles everything into a seurat object ready to be integrated togather in a subsequent script.


args <- commandArgs(TRUE)
path<-args[1]

file<-args[2]

#---- Load Libraries -------------------------------------

source(paste0(path,'/bin/auxiliary/DirectoryChecker.R'))
source(paste0(path,'/bin/auxiliary/PackageLoader.R'))
source(paste0(path,'/bin/auxiliary/LaminarDepth.R'))


#check ram
RAM<-get_ram()
print(paste(RAM,'RAM passed to Rscript'))

ncor=detectCores()-1
print(paste(ncor,'cores detected'))
plan("multiprocess", workers = 1)
#options(future.globals.maxSize = 15000 * 1024^2)



# ---- Fetch File prefix for correct well metadata names -------------------------------------
print(file)
chunks<-unlist(str_split(file,'/'))
filename<-chunks[length(chunks)]

print(filename)
chunks<-unlist(str_split(filename,'_'))
chunks<-chunks[1:length(chunks)-1]
prefix<-paste(chunks,collapse = '_')
timepoint<-chunks[1]
well<-chunks[2]



#Check if directory exists and make it if it doesn't
if(!dir.exists(paste0(path,'/output/spatial/',timepoint))){
  dir.create(paste0(path,'/output/spatial/',timepoint))
}

#Check if directory exists and make it if it doesn't
if(!dir.exists(paste0(path,'/output/spatial/',timepoint,'/',well))){
  dir.create(paste0(path,'/output/spatial/',timepoint,'/',well))
}



# ---- Load Data --------------------------------------


#load the count matrix
mat<-read.csv(file, header = T)
rownames(mat)<-toupper(mat[,1])
mat<-mat[,-1]

#load the region and postion metadata
region<-read.csv(paste0(path,'/input/developmental_single_cell/spatial/',prefix,'_index_region_layer.csv'), header=T)


#set flags for when region/layer info exists
hasRegion<-'Region' %in% colnames(region)
hasLayers<-'Layers' %in% colnames(region)

#check to see if image exists
imgpath<-NA
if(file.exists(paste0(path,'/input/developmental_single_cell/spatial/',prefix,'_overlay_image.png'))){
  #set image path
  imgpath<-paste0(path,'/input/developmental_single_cell/spatial/',prefix,'_overlay_image.png')
  
  #load the image scaling factors
  sf <-as.data.frame(read_csv(paste0(path,"/input/supplementary/",timepoint,"_spatial_sf.csv"),  col_names = FALSE))
  rownames(sf)<-sf[,1]
  sf<-sf[,-1, drop=FALSE]
}


# ---- Build Seurat Object --------------------------------------

#Some cells are regional borderline and have duplciate entries, coin flip randomly allocating borderline cells
dups<-names(which(table(region$Index)>1))

for(dup in dups){
  #randomly discard one of the duplicates
  region<-region[-sample(grep(paste0('^',dup,'$'),region$Index),1),]
}

region<-region[order(region$Index),]  
rownames(region)<-colnames(mat)


#Make seurat object
obj<-CreateSeuratObject(mat,meta.data = region)
obj$orig.ident <-well

#Add realspace coordinates: note they are rotated the wrong way, flip 'em
emb<-cbind(obj$X,-obj$Y)
rownames(emb)<-colnames(obj)
colnames(emb)<-c('Coord_1','Coord_2')
obj[['coords']]<-CreateDimReducObject(embeddings=emb,key='Coord_')
obj$col_coord<-obj$X
obj$row_coord<-obj$Y
obj@meta.data<-cbind(obj@meta.data,region)

#generalise region labels
if(hasRegion){
  obj$Region<-as.character(obj$Region)
  obj$GeneralRegion<-obj$Region
  obj$GeneralRegion[grep('^V2',obj$Region)]<-'V2'
  obj$GeneralRegion[grep('^S1',obj$Region)]<-'S1'
}else{
  obj$Region <- NA
  obj$GeneralRegion <- NA
}

if(hasLayers){
  obj$Layers<-as.character(obj$Layers)
  obj$GeneralLayers<-obj$Layers
  obj$GeneralLayers[grep('2|3|4',obj$Layers)]<-'L234'
  
  #compute depth
  obj<-LaminarDepth(obj)
  
  #Assign strata
  obj$layer_strata <- NA
  obj$layer_strata[which(obj$GeneralLayers %in% c('L1','L234'))] <- 'Upper'
  obj$layer_strata[which(obj$GeneralLayers %in% c('L5','L6'))] <- 'Lower'
  
}else{
  obj$Layers <- NA
  obj$GeneralLayers <- NA
}



# ---- Add image to object --------------------------------------
if(!is.na(imgpath)){
  #read in png as array
  png<-png::readPNG(imgpath)
  
  
  resolution <- 1
  
  #compute average radius of points approximating (width +height)/2 as diameter
  #spotRad<-mean(obj$Width+obj$Height / 2 )
  spotRad<-100
  coords<-obj@meta.data [,c('row_coord','col_coord')]
  colnames(coords)<-c('imagerow','imagecol')
  
  #fix low coords
  coords[which(coords[,1]<0),1]<-0
  
  fiducial<-spotRad*2
  hires<-1
  lowres <- sf[well,]
  spot<-hires
  sf<-scalefactors(spot =spot, fiducial = fiducial, hires = hires, lowres = lowres)
  sr<-(sf$fiducial * sf$lowres)/max(dim(png))
  
  img <- new(
    Class = 'VisiumV1',
    image = png, # an array with the image data
    scale.factors = sf, # scale factors
    coordinates = coords, # data frame with rows as spots and columns as axes
    spot.radius = sr, # radius size
    key='slice_img',
    assay='RNA'
  )
  
  obj@images[[unique(obj$orig.ident)]]<-img

}



# ---- Standard Seurat processing --------------------------------------
obj<-NormalizeData(obj)
obj<-FindVariableFeatures(obj)
obj<-ScaleData(obj, features = rownames(obj))
obj<-RunPCA(obj, npcs= 20)
obj<-FindNeighbors(obj)
# obj<-FindClusters(obj, resolution = resolution)
obj<-RunUMAP(obj,dims=1:15)


#save to directory for reference
save(obj,file=paste(path,'/output/spatial/',timepoint,'/',well,'/',prefix,'_processed','.Rdata',sep= ''))

#pass down pipeline
save(obj,file=paste(prefix,'_processed','.Rdata',sep= ''))



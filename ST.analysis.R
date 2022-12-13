####process ST data by STutility framework
####devtools::install_github("jbergenstrahle/STUtility")
####load pcakges####
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(hdf5r)
library(STutility)
library(future)
library(spdep)
library(clustree)
plan('multisession',workers = 1) ###do parallel
set.seed(101)
getwd()
source('./function/colorPalettes.R')
####read P1####
df1 <- data.frame(samples = './raw data/p1/filtered_feature_bc_matrix.h5',
                 spotfiles = './raw data/p1/spatial/tissue_positions_list.csv',
                 imgs  = './raw data/p1/spatial/tissue_hires_image.png',
                 json = './raw data/p1/spatial/scalefactors_json.json')
p1 <- InputFromTable(infotable = df1, 
                     platform =  "Visium")
p1$sample <- 'primary1'
p1@project.name <- 'primary1'
p1$group <- 'primary'
p1 <- LoadImages(p1,time.resolve = F, verbose = T)
####read R1####
df2 <- data.frame(samples = './raw data/r1/filtered_feature_bc_matrix.h5',
                  spotfiles = './raw data/r1/spatial/tissue_positions_list.csv',
                  imgs  = './raw data/r1/spatial/tissue_hires_image.png',
                  json = './raw data/r1/spatial/scalefactors_json.json')
r1 <- InputFromTable(infotable = df2, 
                     platform =  "Visium")
r1
r1$sample <- 'recurrent1'
r1@project.name <- 'recurrent1'
r1$group <- 'recurrent'
r1 <- LoadImages(r1,time.resolve = F, verbose = T)
####read P2####
df3 <- data.frame(samples = './raw data/p2/filtered_feature_bc_matrix.h5',
                  spotfiles = './raw data/p2/spatial/tissue_positions_list.csv',
                  imgs  = './raw data/p2/spatial/tissue_hires_image.png',
                  json = './raw data/p2/spatial/scalefactors_json.json')
p2 <- InputFromTable(infotable = df3, 
                     platform =  "Visium")
p2
p2$sample <- 'primary2'
p2@project.name <- 'primary2'
p2$group <- 'primary'
p2 <- LoadImages(p2,time.resolve = F, verbose = T)
####read R2####
df4 <- data.frame(samples = './raw data/r2/filtered_feature_bc_matrix.h5',
                  spotfiles = './raw data/r2/spatial/tissue_positions_list.csv',
                  imgs  = './raw data/r2/spatial/tissue_hires_image.png',
                  json = './raw data/r2/spatial/scalefactors_json.json')
r2 <- InputFromTable(infotable = df4, 
                     platform =  "Visium")
r2
r2$sample <- 'recurrent2'
r2@project.name <- 'recurrent2'
r2$group <- 'recurrent'
r2 <- LoadImages(r2,time.resolve = F, verbose = T)
####quality control######
dir.create('./1.Quality control')
obj_list <- list(p1 = p1,
                 r1 = r1,
                 p2 = p2,
                 r2 = r2
                 )
pdf('./1.Quality control/qc.nCount.before.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "nCount_RNA",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
  })
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'nCount_Spatial')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'nCount_RNA',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1,max.cutoff = round(quantile(obj_list[[i]]$nCount_RNA,0.95))) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'nCount_Spatial')
})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'nCount_Spatial')
dev.off()

pdf('./1.Quality control/qc.nFeature.before.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "nFeature_RNA",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'nFeature_Spatial')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'nFeature_RNA',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1,max.cutoff = round(quantile(obj_list[[i]]$nFeature_RNA,0.95))) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'nFeature_Spatial')})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'nFeature_Spatial')
dev.off()

for(i in 1:length(obj_list)){
  mt.genes <- grep(pattern = "^MT", x = rownames(obj_list[[i]]), value = TRUE)
  obj_list[[i]]$percent.mito <- (Matrix::colSums(obj_list[[i]]@assays$RNA@counts[mt.genes, ])/Matrix::colSums(obj_list[[i]]@assays$RNA@counts))*100
}
pdf('./1.Quality control/qc.percentMT.before.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "percent.mito",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'percent.mito')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'percent.mito',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'percent.mito')})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'percent.mito')
dev.off()

pdf('./1.Quality control/scatterplot.before.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- FeatureScatter(obj_list[[i]],feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA',plot.cor = T) + NoLegend() + labs(subtitle = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2)
gg <- lapply(names(obj_list),function(i){
  plot <- FeatureScatter(obj_list[[i]],feature1 = 'nCount_RNA',feature2 = 'percent.mito',plot.cor = T) + NoLegend() + labs(subtitle = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2)
dev.off()
####subset data####
####do not remove too much spots
####nFeature > 100&percent.mito<20
for(i in names(obj_list)){
  obj_list[[i]] <- SubsetSTData(obj_list[[i]], expression = nFeature_RNA > 100 & percent.mito < 20)
  cat("Spots removed: ", ncol(get(i)) - ncol(obj_list[[i]]), "\n")
}
####observe the qc metrics after filter####
pdf('./1.Quality control/qc.nCount.after.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "nCount_RNA",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'nCount_Spatial')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'nCount_RNA',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1,max.cutoff = round(quantile(obj_list[[i]]$nCount_RNA,0.95))) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'nCount_Spatial')
})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'nCount_Spatial')
dev.off()

pdf('./1.Quality control/qc.nFeature.after.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "nFeature_RNA",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'nFeature_Spatial')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'nFeature_RNA',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1,max.cutoff = round(quantile(obj_list[[i]]$nFeature_RNA,0.95))) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'nFeature_Spatial')})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'nFeature_Spatial')
dev.off()


pdf('./1.Quality control/qc.percentMT.after.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- VlnPlot(obj_list[[i]],features = "percent.mito",pt.size = 0,group.by = 'sample') + NoLegend() + xlab("") + labs(title = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = 'percent.mito')
gg <- lapply(names(obj_list),function(i){
  plot <- ST.FeaturePlot(obj_list[[i]],features = 'percent.mito',dark.theme = T,show.sb = F,palette = 'Spectral',pt.size = 1) + theme(legend.position = "right",strip.text = element_blank())+ 
    labs(title = obj_list[[i]]@project.name,fill = 'percent.mito')})
plot_grid(plotlist = gg,ncol = 2)+ plot_annotation(title = 'percent.mito')
dev.off()

pdf('./1.Quality control/scatterplot.after.filter.pdf')
gg <- lapply(names(obj_list),function(i){
  plot <- FeatureScatter(obj_list[[i]],feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA',plot.cor = T) + NoLegend() + labs(subtitle = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2)
gg <- lapply(names(obj_list),function(i){
  plot <- FeatureScatter(obj_list[[i]],feature1 = 'nCount_RNA',feature2 = 'percent.mito',plot.cor = T) + NoLegend() + labs(subtitle = obj_list[[i]]@project.name)
})
plot_grid(plotlist = gg,ncol = 2)
dev.off()
#####归一化（推荐使用sctransform)####
for(i in 1:length(obj_list)){
  obj_list[[i]] <- SCTransform(obj_list[[i]],assay = 'RNA',return.only.var.genes = F)
  cat(obj_list[[i]]@project.name," finished","\n")
}
for(i in 1:length(obj_list)){
  DefaultAssay(obj_list[[i]]) <- 'SCT'
}
####降维，聚类####
for(i in 1:length(obj_list)){
  obj_list[[i]] <- RunPCA(obj_list[[i]],assay = 'SCT',npcs = 100,verbose = F)
}
dir.create('./2.Reduction')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
pdf('./2.Reduction/elbow.plot.pdf')
gg <- lapply(1:length(obj_list),function(i){
  plot <- do_Elbow_quantitative(obj_list[[i]],harmony = F)
  plot$fig + labs(title = paste(obj_list[[i]]@project.name,'pc',plot$pcs))
})
plot_grid(plotlist = gg)
dev.off()

for(i in 1:length(obj_list)){
  obj_list[[i]] <- FindNeighbors(obj_list[[i]],reduction = "pca",dims = 1:20,verbose = F)
  obj_list[[i]] <- FindClusters(obj_list[[i]],resolution = seq(0.1,1,0.1),verbose = F)
  obj_list[[i]] <- RunUMAP(obj_list[[i]],reduction = "pca",assay = "SCT",dims = 1:20,verbose = F)
  cat(obj_list[[i]]@project.name," finished","\n")
}

for(i in 1:length(obj_list)){
  pdf(file = paste0('./2.Reduction/',obj_list[[i]]@project.name,' res.observe.pdf'),10,10)
  {gg <- clustree(obj_list[[i]])
  print(gg)
  merge.res <- lapply(seq(0.1,1,0.1), function(x){
    p <- DimPlot(object = obj_list[[i]], reduction = 'umap',label = TRUE, group.by = paste0('SCT', "_snn_res.", x)) + NoLegend()
    print(p)
  })}
  dev.off()
  cat(obj_list[[i]]@project.name," finished","\n")
}

p1 <- obj_list[['p1']]
p2 <- obj_list[['p2']]
r1 <- obj_list[['r1']]
r2 <- obj_list[['r2']]

p1$seurat_clusters <- p1$SCT_snn_res.0.5
p2$seurat_clusters <- p2$SCT_snn_res.0.3
r1$seurat_clusters <- r1$SCT_snn_res.0.4
r2$seurat_clusters <- r2$SCT_snn_res.0.3
obj_list <- list(p1 = p1,
                 r1 = r1,
                 p2 = p2,
                 r2 = r2
)
for(i in 1:length(obj_list)){
  pdf(file = paste0('./2.Reduction/',obj_list[[i]]@project.name,' dimplot.pdf'),10,10)
   {gg1 <- DimPlot(object = obj_list[[i]], reduction = 'umap',label = TRUE, group.by = 'seurat_clusters',cols = Palettes[['circus']])
   print(gg1)
   gg2 <- ST.FeaturePlot(obj_list[[i]],dark.theme = F,features = 'seurat_clusters',pt.size = 2.7,show.sb = F,cols = Palettes[['circus']]) + theme(strip.text = element_blank()) + 
     guides(fill = guide_legend(override.aes = list(size = 8))) + ggtitle(obj_list[[i]]@project.name)
   print(gg2)
   gg3 <- ST.FeaturePlot(obj_list[[i]],dark.theme = F,features = 'seurat_clusters',pt.size = 1.8,show.sb = F,split.labels = T,ncol = 3,cols = Palettes[['circus']]) + 
      ggtitle(obj_list[[i]]@project.name) + NoLegend()
   print(gg3)
   dev.off()}
   cat(obj_list[[i]]@project.name," finished","\n")
}

save(obj_list,file = './2.Reduction/obj.list.rData')
####无监督寻找空间高变基因（CorSpatialGenes) by spdep package####
load('./2.Reduction/obj.list.rData')
for(i in 1:length(obj_list)){
  obj_list[[i]] <- MaskImages(obj_list[[i]])
  cat(obj_list[[i]]@project.name," finished","\n")
}

for(i in 1:length(obj_list)){
  spatgenes <- CorSpatialGenes(obj_list[[i]],assay = "SCT",features = VariableFeatures(obj_list[[i]])[1:1000])
  write.csv(spatgenes,file = paste0('./2.Reduction/',obj_list[[i]]@project.name,' spatially.correlated.genes.csv'))
  top6_spfeatures <- head(spatgenes$gene,n = 6)
  obj <- obj_list[[i]]
  gg <- lapply(1:length(top6_spfeatures),function(x){
    plot <- FeatureOverlay(obj,show.sb = F,features = top6_spfeatures[x],type = 'raw',add.alpha = T,palette = 'heat',pt.size = 1.8) + labs(subtitle = "")
  })
  pdf(paste0('./2.Reduction/',obj_list[[i]]@project.name,' spatially.correlated.genes.pdf'),width = 10,height = 15)
  gg <- plot_grid(plotlist = gg,ncol = 2) + plot_annotation(title = paste('Spatially correlated genes of',obj@project.name))
  print(gg)
  dev.off()
  cat(obj_list[[i]]@project.name," finished","\n")
  }
save(obj_list,file = './2.Reduction/obj.list.rData')

p1 <- obj_list[['p1']]
r1 <- obj_list[['r1']]
p2 <- obj_list[['p2']]
r2 <- obj_list[['r2']]
save(p1,r1,p2,r2,file = './2.Reduction/all.samples.before.RCTD.rData')

####run RCTD for deconvolution####
load('./2.Reduction/all.samples.before.RCTD.rData')
library(spacexr)
library(Matrix)
dir.create('./3.RCTD')

####prepare RCTD object for P1####
dir.create('./3.RCTD/p1')
ct <- as.data.frame(p1@assays$RNA@counts)
coord <- as.data.frame(p1@tools$Staffli@meta.data)
coord <- coord %>% select(original_x,original_y)
coord$original_y <- max(coord$original_y)*2 - coord$original_y

nUMI <- colSums(ct) 
nFeature <- p1$nFeature_RNA
puck <- SpatialRNA(coords = coord, counts = ct, nUMI = nUMI)
barcodes <- colnames(puck@counts)
pdf('./3.RCTD/p1/Spatial.UMI.pdf')
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.95))), 
                          title ='plot of nUMI', size=2, alpha=0.8,small_point = T) 
dev.off()
saveRDS(puck,file = './3.RCTD/p1/primary1.RCTD.obj.rds')
puck.tmp <- SpatialRNA(coords = coord,counts = ct,nUMI = nFeature)
pdf('./3.RCTD/p1/Spatial.nFeature.pdf')
plot_puck_continuous(puck.tmp, barcodes, puck.tmp@nUMI, ylimit = c(0,round(quantile(puck.tmp@nUMI,0.95))), 
                     title ='plot of nFeature', size=2, alpha=0.8,small_point = T) 
dev.off()
####prepare RCTD object for R1####
dir.create('./3.RCTD/r1')
ct <- as.data.frame(r1@assays$RNA@counts)
coord <- as.data.frame(r1@tools$Staffli@meta.data)
coord <- coord %>% select(original_x,original_y)
coord$original_y <- max(coord$original_y)*2 - coord$original_y

nUMI <- colSums(ct) 
nFeature <- r1$nFeature_RNA
puck <- SpatialRNA(coords = coord, counts = ct, nUMI = nUMI)
barcodes <- colnames(puck@counts)
pdf('./3.RCTD/r1/Spatial.UMI.pdf')
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.95))), 
                     title ='plot of nUMI', size=2, alpha=0.8,small_point = T) 
dev.off()
saveRDS(puck,file = './3.RCTD/r1/recurrent1.RCTD.obj.rds')
puck.tmp <- SpatialRNA(coords = coord,counts = ct,nUMI = nFeature)
pdf('./3.RCTD/r1/Spatial.nFeature.pdf')
plot_puck_continuous(puck.tmp, barcodes, puck.tmp@nUMI, ylimit = c(0,round(quantile(puck.tmp@nUMI,0.95))), 
                     title ='plot of nFeature', size=2, alpha=0.8,small_point = T) 
dev.off()
####prepare RCTD object for P2####
dir.create('./3.RCTD/p2')
ct <- as.data.frame(p2@assays$RNA@counts)
coord <- as.data.frame(p2@tools$Staffli@meta.data)
coord <- coord %>% select(original_x,original_y)
coord$original_y <- max(coord$original_y)*2 - coord$original_y

nUMI <- colSums(ct) 
puck <- SpatialRNA(coords = coord, counts = ct, nUMI = nUMI)
barcodes <- colnames(puck@counts)
pdf('./3.RCTD/p2/Spatial.UMI.pdf')
plot_puck_continuous(puck3, barcodes, puck3@nUMI, ylimit = c(0,round(quantile(puck3@nUMI,0.95))), 
                     title ='plot of nUMI', size=2, alpha=0.8,small_point = T) 
dev.off()
saveRDS(puck3,file = './3.RCTD/p2/primary2.RCTD.obj.rds')
nFeature <- p2$nFeature_RNA
puck.tmp <- SpatialRNA(coords = coord,counts = ct,nUMI = nFeature)
pdf('./3.RCTD/p2/Spatial.nFeature.pdf')
plot_puck_continuous(puck.tmp, barcodes, puck.tmp@nUMI, ylimit = c(0,round(quantile(puck.tmp@nUMI,0.95))), 
                     title ='plot of nFeature', size=2, alpha=0.8,small_point = T) 
dev.off()
####prepare RCTD object for R2####
dir.create('./3.RCTD/r2')
ct <- as.data.frame(r2@assays$RNA@counts)
coord <- as.data.frame(r2@tools$Staffli@meta.data)
coord <- coord %>% select(original_x,original_y)
coord$original_y <- max(coord$original_y)*2 - coord$original_y

nUMI <- colSums(ct) 
puck <- SpatialRNA(coords = coord, counts = ct, nUMI = nUMI)
barcodes <- colnames(puck@counts)
pdf('./3.RCTD/r2/Spatial.UMI.pdf')
plot_puck_continuous(puck4, barcodes, puck4@nUMI, ylimit = c(0,round(quantile(puck4@nUMI,0.95))), 
                     title ='plot of nUMI', size=2, alpha=0.8,small_point = T) 
dev.off()
saveRDS(puck4,file = './3.RCTD/r2/recurrent2.RCTD.obj.rds')
nFeature <- r2$nFeature_RNA
puck.tmp <- SpatialRNA(coords = coord,counts = ct,nUMI = nFeature)
pdf('./3.RCTD/r2/Spatial.nFeature.pdf')
plot_puck_continuous(puck.tmp, barcodes, puck.tmp@nUMI, ylimit = c(0,round(quantile(puck.tmp@nUMI,0.95))), 
                     title ='plot of nFeature', size=2, alpha=0.8,small_point = T) 
dev.off()

####RCTD for large celltypes (do not devide T celltypes)####
####prepare the single cell reference(scRef_NC)
dir.create('./3.RCTD/T.not.devide')
sc_ref <- readRDS('./3.RCTD/scRef_NC/3.Annotation/annotated.filtered.obj.rds')
sc_ref
table(sc_ref$celltype.refined)
DimPlot(sc_ref,group.by = 'celltype.refined',reduction = 'umap',label = T)
counts <- GetAssayData(sc_ref,assay = 'RNA',slot = 'counts') %>% as.data.frame()
metadata <- sc_ref@meta.data
celltype <- metadata$celltype.refined
celltype[celltype%in%c('CD4','CD8','Treg')] <- 'T cells'
celltype[celltype=='Macro/mono'] <- 'Macro_mono'
names(celltype) <- rownames(metadata)
celltype <- as.factor(celltype)
table(celltype)
nUMI <- metadata$nUMI
names(nUMI) <- rownames(metadata)
### Create the Reference object
reference <- Reference(counts, celltype, nUMI)
str(reference)
print(dim(reference@counts)) 
#number of occurences for each cell type
table(reference@cell_types) 
saveRDS(reference,file = './3.RCTD/T.not.devide/nc.reference.obj.rds')
####Run RCTD
reference <- readRDS('./3.RCTD/T.not.devide/nc.reference.obj.rds')
p1 <- readRDS('./3.RCTD/p1/primary1.RCTD.obj.rds')
r1 <- readRDS('./3.RCTD/r1/recurrent1.RCTD.obj.rds')
p2 <- readRDS('./3.RCTD/p2/primary2.RCTD.obj.rds')
r2 <- readRDS('./3.RCTD/r2/recurrent2.RCTD.obj.rds')
obj_list <- list(p1 = p1,
                 r1 = r1,
                 p2 = p2,
                 r2 = r2)
rctd_list <- lapply(names(obj_list),function(i){
  myRCTD <- create.RCTD(obj_list[[i]], reference, max_cores = 11)
  cat(names(obj_list[i])," finished","\n")
  return(myRCTD)
})
names(rctd_list) <- names(obj_list)
str(rctd_list[['p2']])

res_rctd <- lapply(names(rctd_list),function(i){
  res <- run.RCTD(RCTD = rctd_list[[i]], doublet_mode = 'full')
  cat(names(rctd_list[i])," finished","\n")
  return(res)
})

names(res_rctd) <- names(rctd_list)

p1.res <- res_rctd[['p1']]
r1.res <- res_rctd[['r1']]
p2.res <- res_rctd[['p2']]
r2.res <- res_rctd[['r2']]
save(p1.res,r1.res,p2.res,r2.res,file = './3.RCTD/T.not.devide/RCTD.res.rData')

####RCTD for the detailed celltypes(devide T celltype)####
####prepare the single cell reference(scRef_NC)
dir.create('./3.RCTD/T.devide')
sc_ref <- readRDS('./3.RCTD/scRef_NC/3.Annotation/annotated.filtered.obj.rds')
sc_ref
table(sc_ref$celltype.refined)
DimPlot(sc_ref,group.by = 'celltype.refined',reduction = 'umap',label = T)
counts <- GetAssayData(sc_ref,assay = 'RNA',slot = 'counts') %>% as.data.frame()
metadata <- sc_ref@meta.data
celltype <- metadata$celltype.refined
celltype[celltype=='Macro/mono'] <- 'Macro_mono'
names(celltype) <- rownames(metadata)
celltype <- as.factor(celltype)
table(celltype)
nUMI <- metadata$nUMI
names(nUMI) <- rownames(metadata)
### Create the Reference object
reference <- Reference(counts, celltype, nUMI)
str(reference)
print(dim(reference@counts)) 
#number of occurences for each cell type
table(reference@cell_types) 
saveRDS(reference,file = './3.RCTD/T.devide/nc.reference.obj.rds')
####Run RCTD
reference <- readRDS('./3.RCTD/T.devide/nc.reference.obj.rds')
p1 <- readRDS('./3.RCTD/p1/primary1.RCTD.obj.rds')
r1 <- readRDS('./3.RCTD/r1/recurrent1.RCTD.obj.rds')
p2 <- readRDS('./3.RCTD/p2/primary2.RCTD.obj.rds')
r2 <- readRDS('./3.RCTD/r2/recurrent2.RCTD.obj.rds')
obj_list <- list(p1 = p1,
                 r1 = r1,
                 p2 = p2,
                 r2 = r2)
rctd_list <- lapply(names(obj_list),function(i){
  myRCTD <- create.RCTD(obj_list[[i]], reference, max_cores = 11)
  cat(names(obj_list[i])," finished","\n")
  return(myRCTD)
})
names(rctd_list) <- names(obj_list)
str(rctd_list[['p2']])

res_rctd <- lapply(names(rctd_list),function(i){
  res <- run.RCTD(RCTD = rctd_list[[i]], doublet_mode = 'full')
  cat(names(rctd_list[i])," finished","\n")
  return(res)
})

names(res_rctd) <- names(rctd_list)

p1.res <- res_rctd[['p1']]
r1.res <- res_rctd[['r1']]
p2.res <- res_rctd[['p2']]
r2.res <- res_rctd[['r2']]
save(p1.res,r1.res,p2.res,r2.res,file = './3.RCTD/T.devide/RCTD.res.rData')
####merge rctd results and spatial matrix####
rm(list = ls())
source('./function/colorPalettes.R')
load('./2.Reduction/all.samples.before.RCTD.rData')
load('./3.RCTD/T.not.devide/RCTD.res.rData')
####p1
res <- p1.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(p1@meta.data),rownames(norm.res))
p1@meta.data <- cbind(p1@meta.data,norm.res)
####r1
res <- r1.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(r1@meta.data),rownames(norm.res))
r1@meta.data <- cbind(r1@meta.data,norm.res)
####p2
res <- p2.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(p2@meta.data),rownames(norm.res))
p2@meta.data <- cbind(p2@meta.data,norm.res)
####r2
res <- r2.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(r2@meta.data),rownames(norm.res))
r2@meta.data <- cbind(r2@meta.data,norm.res)

load('./3.RCTD/T.devide/RCTD.res.rData')
####p1
res <- p1.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(p1@meta.data),rownames(norm.res))
norm.res <- norm.res[,c('CD4','CD8','Treg')]
p1@meta.data <- cbind(p1@meta.data,norm.res)
####r1
res <- r1.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(r1@meta.data),rownames(norm.res))
norm.res <- norm.res[,c('CD4','CD8','Treg')]
r1@meta.data <- cbind(r1@meta.data,norm.res)
####p2
res <- p2.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(p2@meta.data),rownames(norm.res))
norm.res <- norm.res[,c('CD4','CD8','Treg')]
p2@meta.data <- cbind(p2@meta.data,norm.res)
####r2
res <- r2.res@results$weights
norm.res <- normalize_weights(res) %>% as.data.frame()
identical(rownames(r2@meta.data),rownames(norm.res))
norm.res <- norm.res[,c('CD4','CD8','Treg')]
r2@meta.data <- cbind(r2@meta.data,norm.res)

save(p1,r1,p2,r2,file = './3.RCTD/predicted.obj.rData')

####visualize RCTD results####
load('./3.RCTD/predicted.obj.rData')
####p1
celltypes <- c('Epithelial','Fibroblast','SMC','T cells','CD4','CD8','Treg','B cells','Plasma','NK','Macro_mono','DC','pDC','Mast','Endothelial')
plot_list <- lapply(celltypes,function(pro){
  gg1 <- ST.FeaturePlot(p1,features = pro,dark.theme = T,palette = 'Spectral',max.cutoff = quantile(p1@meta.data$pro,0.9),pt.size = 1.8) + theme(legend.position = "right",strip.text = element_blank()) + 
    labs(color = 'value')
  gg2 <- FeatureOverlay(p1,features = pro,palette = 'Spectral',max.cutoff = quantile(p1@meta.data$pro,0.9),pt.size = 1.8,add.alpha = T,type = 'raw') + labs(subtitle = "",color = 'value')
  gg <- plot_grid(plotlist = list(gg1,gg2),ncol = 2)
  cat(pro," finished","\n")
  return(gg)
  })
pdf('./3.RCTD/p1/predicted.celltyeps.pdf',width = 10,height = 5)
par(mfrow = c(1,2),xpd = T)
print(plot_list)
dev.off()
####r1
celltypes <- c('Epithelial','Fibroblast','SMC','T cells','CD4','CD8','Treg','B cells','Plasma','NK','Macro_mono','DC','pDC','Mast','Endothelial')
plot_list <- lapply(celltypes,function(pro){
  gg1 <- ST.FeaturePlot(r1,features = pro,dark.theme = T,palette = 'Spectral',max.cutoff = quantile(r1@meta.data$pro,0.9),pt.size = 1.8) + theme(legend.position = "right",strip.text = element_blank()) + 
    labs(color = 'value')
  gg2 <- FeatureOverlay(r1,features = pro,palette = 'Spectral',max.cutoff = quantile(r1@meta.data$pro,0.9),pt.size = 1.8,add.alpha = T,type = 'raw') + labs(subtitle = "",color = 'value')
  gg <- plot_grid(plotlist = list(gg1,gg2),ncol = 2)
  cat(pro," finished","\n")
  return(gg)
})
pdf('./3.RCTD/r1/predicted.celltyeps.pdf',width = 10,height = 5)
par(mfrow = c(1,2),xpd = T)
print(plot_list)
dev.off()
####p2
celltypes <- c('Epithelial','Fibroblast','SMC','T cells','CD4','CD8','Treg','B cells','Plasma','NK','Macro_mono','DC','pDC','Mast','Endothelial')
plot_list <- lapply(celltypes,function(pro){
  gg1 <- ST.FeaturePlot(p2,features = pro,dark.theme = T,palette = 'Spectral',max.cutoff = quantile(p2@meta.data$pro,0.9),pt.size = 1.8) + theme(legend.position = "right",strip.text = element_blank()) + 
    labs(color = 'value')
  gg2 <- FeatureOverlay(p2,features = pro,palette = 'Spectral',max.cutoff = quantile(p2@meta.data$pro,0.9),pt.size = 1.8,add.alpha = T,type = 'raw') + labs(subtitle = "",color = 'value')
  gg <- plot_grid(plotlist = list(gg1,gg2),ncol = 2)
  cat(pro," finished","\n")
  return(gg)
})
pdf('./3.RCTD/p2/predicted.celltyeps.pdf',width = 10,height = 5)
par(mfrow = c(1,2),xpd = T)
print(plot_list)
dev.off()
####r2
celltypes <- c('Epithelial','Fibroblast','SMC','T cells','CD4','CD8','Treg','B cells','Plasma','NK','Macro_mono','DC','pDC','Mast','Endothelial')
plot_list <- lapply(celltypes,function(pro){
  gg1 <- ST.FeaturePlot(r2,features = pro,dark.theme = T,palette = 'Spectral',max.cutoff = quantile(r2@meta.data$pro,0.9),pt.size = 1.8) + theme(legend.position = "right",strip.text = element_blank()) + 
    labs(color = 'value')
  gg2 <- FeatureOverlay(r2,features = pro,palette = 'Spectral',max.cutoff = quantile(r2@meta.data$pro,0.9),pt.size = 1.8,add.alpha = T,type = 'raw') + labs(subtitle = "",color = 'value')
  gg <- plot_grid(plotlist = list(gg1,gg2),ncol = 2)
  cat(pro," finished","\n")
  return(gg)
})
pdf('./3.RCTD/r2/predicted.celltyeps.pdf',width = 10,height = 5)
par(mfrow = c(1,2),xpd = T)
print(plot_list)
dev.off()

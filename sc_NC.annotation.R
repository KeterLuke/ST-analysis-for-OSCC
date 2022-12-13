#### load package (always load all the package before you run the code)####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
library(colourpicker)
library(cowplot)
plan("multisession", workers = 1) ####do parallel
options(future.globals.maxSize = 80000 * 1024^2)
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
####read samples####
files <- list.files('./3.RCTD/scRef_NC/raw data/')
files <- files[grep(pattern = 'PBL',files,invert = T)]
dir <- './3.RCTD/scRef_NC/raw data'
sceList = lapply(files,function(pro){ 
  # pro=samples[1] 
  print(pro)
  ct <- Read10X(data.dir = file.path(dir,pro))
  sce=CreateSeuratObject(counts =  ct ,
                         project =  pro,
                         min.cells = 3,
                         min.features = 200)
  
  return(sce)
})
names(sceList) <- files
sce = merge(sceList[[1]], y = sceList[2:length(sceList)],add.cell.ids <- names(sceList))
sce$sample <- sce$orig.ident
table(sce$sample)
####quality control####
setwd('./3.RCTD/scRef_NC/')
getwd()
sce$nUMI <- sce$nCount_RNA
sce$nGene <- sce$nFeature_RNA
sce$log10GenesPerUMI <- log10(sce$nGene)/log10(sce$nUMI)
## Mitochondrial gene ratio
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce$mito.ratio <- sce$percent.mt/100
## Ribosome gene ratio
sce[["percent.rb"]] <- PercentageFeatureSet(sce, pattern = "^RPL|^RPS")
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
## HBC gene ratio
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HBA|^HBB")
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
dir.create('./1.QualityControl')
sce$project <- 'NC_tumor'
Idents(sce) <- sce$project
pdf(file = "1.QualityControl/count.feature.mt.pdf",width = 10)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"),ncol = 3,pt.size = 0,group.by = 'project')
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'project')
ploT1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'project')
plot1 + ploT1

ggdensity(sce@meta.data, x = "nCount_RNA", title = "nCount")
ggdensity(sce@meta.data, x = "nFeature_RNA", title = "nFeature") + geom_vline(xintercept = 500)
ggdensity(sce@meta.data, x = "log10GenesPerUMI", title = "log10GenesPerUMI") + geom_vline(xintercept = 0.8)
ggdensity(sce@meta.data, x = "percent.rb", title = "percent.rb") +  geom_vline(xintercept = 8)
dev.off()

source('./function/do_qc_plot.R')
do_qc_plot(sce,dir = './1.QualityControl/',name = 'qcplot.before.pdf')
saveRDS(sce,file = './merged.unfiltered.obj.rds')
####subset:nFeatures 300-5000,nCount>500,percent.mt<10,complexity>0.8####
sce_obj <- subset(sce, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 500 &
                    log10GenesPerUMI > 0.8)
do_qc_plot(sce_obj,dir = './1.QualityControl/',name = 'qcplot.after.pdf')
pdf("1.QualityControl/filtered.statistics.pdf")
source('./function/colorPalettes.R')
VlnPlot(object = sce_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rb','percent.hb'), ncol = 2, pt.size = 0)
FeatureScatter(object = sce_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = sce_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
saveRDS(sce_obj,file = './1.QualityControl/filtered.obj.rds')
####normalize and scaledata####
setwd('./3.RCTD/scRef_NC/')
getwd()
sce_obj <- readRDS('./1.QualityControl/filtered.obj.rds')
sce_obj <- sce_obj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
sce_obj <- RunPCA(sce_obj,features = VariableFeatures(sce_obj))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sce_obj <- CellCycleScoring(sce_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sce_obj <- RunPCA(sce_obj, features = c(s.genes, g2m.genes),verbose = F)
DimPlot(sce_obj,group.by = 'Phase',reduction = 'pca')
sce_obj <- ScaleData(sce_obj,vars.to.regress = c('S.Score','G2M.Score'))
sce_obj <- RunPCA(sce_obj,features = VariableFeatures(sce_obj))
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(sce_obj,harmony = F)
DimPlot(sce_obj,reduction = 'pca',group.by = 'orig.ident')
sce_obj$patient <- sapply(strsplit(sce_obj$orig.ident,'_'),"[",1)
table(sce_obj$orig.ident,sce_obj$patient)
DimPlot(sce_obj,reduction = 'pca',group.by = 'patient')
####integrate(harmony)####
sce_obj <- RunHarmony(sce_obj,group.by.vars = 'orig.ident',reduction = 'pca',verbose = F)
DimPlot(sce_obj,reduction = 'harmony',group.by = 'orig.ident')
DimPlot(sce_obj,reduction = 'harmony',group.by = 'patient')
dir.create('./2.Harmony')
pdf('./2.Harmony/observe.batcheffect.pdf')
DimPlot(sce_obj,reduction = 'pca',group.by = 'orig.ident')
DimPlot(sce_obj,reduction = 'harmony',group.by = 'orig.ident')
dev.off()
do_Elbow_quantitative(sce_obj,harmony = T)

pc <- 50
sce_obj <- FindNeighbors(sce_obj,reduction = 'harmony',dims = 1:pc,verbose = F)
set.resolutions <- seq(0.1,1,0.1)
sce_obj <- FindClusters(sce_obj,resolution = set.resolutions,verbose = F)
sce_obj <- RunUMAP(sce_obj,reduction = 'harmony',dims = 1:pc,verbose = F)

pdf('./2.Harmony/Harmony.integrate.pc50.pdf')
clustree(sce_obj)
DimPlot(sce_obj,reduction = 'umap',group.by = 'orig.ident')
merge.res <- sapply(set.resolutions,function(i){
  p <- DimPlot(sce_obj,reduction = 'umap',group.by = paste0('RNA_snn_res.',i),label = T)
  print(p)
})
dev.off()
####select:pc50 res0.8####
sce_obj$seurat_clusters <- sce_obj$RNA_snn_res.0.8
table(sce_obj$seurat_clusters)
sce_obj$seurat_clusters <- factor(sce_obj$seurat_clusters,levels = as.character(0:28))
table(sce_obj$seurat_clusters)
pdf('./2.Harmony/cluster.pdf')
DimPlot(sce_obj,group.by = 'orig.ident',reduction = 'umap')
DimPlot(sce_obj,group.by = 'seurat_clusters',label = T,reduction = 'umap')
dev.off()
saveRDS(sce_obj,file = './2.Harmony/integrated.obj.rds')

####annotation####
data.merge <- readRDS('./2.Harmony/integrated.obj.rds')
dir.create('./3.Annotation')
cell.type.markers <- read.table(file = "3.Annotation/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) ##perform 0-1 standardization for all clusters per gene
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$Celltype, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)##perform 0-1 standardization for all clusters per celltype marker
annotation.colors <- Palettes$stallion2[1:length(unique(cell.type.markers$Celltype))]
names(annotation.colors) <- unique(cell.type.markers$Celltype)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Celltype, 
                                               levels = unique(cell.type.markers$Celltype)),
                                 col = list(Type = annotation.colors),show_annotation_name = F)
pdf("3.Annotation/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 1.5), c("grey", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$Celltype, levels = unique(cell.type.markers$Celltype))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, left_annotation = row.annotations,
        width = unit(10, "cm"), height = unit(16, "cm"), cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
a <- data.merge
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(T = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='T'],
                  B = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='B'],
                  cycling = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Cycling'],
                  Plasma = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Plasma'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Myeloid = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Myeloid'],
                  Epi = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Epi'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo'],
                  Fibro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Fibro'],
                  Smooth = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Smooth'])
library(GEB)
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "seurat_clusters",dot.scale = 4)
dev.off()
pdf('./3.Annotation/marker.observe.pdf')
FeaturePlot(data.merge,features = c('PTPRC','EPCAM','PECAM1','COL1A1'),order = T,reduction = 'umap',min.cutoff = 1)
FeaturePlot(data.merge,features = c("CD3D",'CD3E','CD3G','NKG7'),order = T,reduction = 'umap',min.cutoff = 1)
FeaturePlot(data.merge,features = c('MS4A1','CD19','CD79A','CD79B'),reduction = 'umap',order = T,min.cutoff = 1)
dev.off()

####find markers####
plan('multisession',workers = 11)
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.Annotation/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.Annotation/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.Annotation/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)
top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("3.Annotation/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(subset(data.merge,downsample = 50), features = unique(top5.genes$gene), size = 2,group.bar = T) + NoLegend()
dev.off()

# Change the column of the resolution if you ended up using a different one than 0.5 
cluster.ids <- sort(as.numeric(unique(as.character(data.merge@meta.data$seurat_clusters))))
# Manually check and annotate each cluster to immune and non-immune   
celltype <- c("CD4","Epithelial","CD8","Treg","CD8","B cells",
                       "NK","Endothelial","Fibroblast","Macrophage","Endothelial",
                       "Proliferating T","Monocyte","Epithelial","CD4","Doublet",
                       "cDC2","Monocyte","CD4","Plasma","CD4",
                       "pDC","Myofibroblast","Epithelial","Proliferating B","LAMP3+DC",
                       "Mast","CD4","CD4")



# Add annotation to the Seurat object 
data.merge@meta.data$celltype <- data.merge@meta.data$seurat_clusters 
data.merge@meta.data$celltype <- plyr::mapvalues(x = data.merge@meta.data$celltype, from = cluster.ids, to = celltype)
table(data.merge$celltype,useNA = 'always')

genelist <- c('CD3D','CD3E','CD3G','IL7R','CD40LG','EPCAM','KRT18','KRT19','CD8A','CD8B','IL2RA','FOXP3','TNFRSF4','MS4A1','CD19','CD79A',
              'CD79B','NKG7','KLRB1','PRF1','GNLY','CLDN5','PECAM1','RAMP2','LUM','DCN','COL1A1','COL1A2','CD14','C1QA','APOE','MKI67','TOP2A',
              'S100A8','S100A9','IL1B','PTGS2','CD1C','CD1E','CLEC10A','FCER1A','MZB1','IGHG1','IGKC','LILRA4','IRF7','ACTA2','TAGLN','MYH11',
              'RGS13','STMN1','LRMP','LAMP3','CCL19','BIRC3','IDO1','TPSAB1','TPSB2')
pdf('./3.Annotation/marker.expression.pdf',width = 15,height = 12)
DotPlot(data.merge,group.by = 'celltype',features = genelist,scale.by = 'size',cols = Palettes[['blueRed']]) + coord_flip() + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(angle = 30))
dev.off()

pdf('./3.Annotation/celltype.pdf',12,12)
DimPlot(data.merge,group.by = 'celltype',label = T,reduction = 'umap')
dev.off()
saveRDS(data.merge,file = './3.Annotation/annotated.all.obj.rds')

###filter doublet and cycling cells####
####to ensure the accuracy of ST annotation, we remove doublets and cycling cells from the reference SC data
####for their probable mixture of other cell types
celltype.filtered <- setdiff(celltype,c('Doublet','Proliferating T','Proliferating B'))
data.merge.filtered <- subset(data.merge,subset = celltype %in% celltype.filtered)
table(data.merge.filtered$celltype)
data.merge.filtered$celltype.refined <- as.character(data.merge.filtered$celltype)
index <- data.merge.filtered$celltype.refined%in%c('cDC2','LAMP3+DC')
data.merge.filtered@meta.data$celltype.refined[index] <- 'DC'
index <- data.merge.filtered$celltype.refined%in%c('Monocyte','Macrophage')
data.merge.filtered@meta.data$celltype.refined[index] <- 'Macro/mono'


data.merge.filtered@meta.data$celltype.refined[data.merge.filtered@meta.data$celltype.refined=='Myofibroblast'] <- 'SMC'
table(data.merge.filtered$celltype.refined)
table(data.merge.filtered$celltype)
pdf('./3.Annotation/celltype.filtered.pdf',12,12)
DimPlot(data.merge.filtered,group.by = 'celltype.refined',label = T,label.size = 5,reduction = 'umap',cols = Palettes[['mycols_19']])
dev.off()


saveRDS(data.merge.filtered,file = './3.Annotation/annotated.filtered.obj.rds')

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(sctransform)
library(SeuratData)

# Try aggregating all replicates

setwd('/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_08/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_08/hg19.exonic+intronic/alignment/Puck_200727_08.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts) <- toupper(x = raw_counts[, 1])
rownames(raw_counts) <- make.names(toupper(x=raw_counts[,1]), unique = TRUE)
colnames(x = raw_counts) <- toupper(x = colnames(x = raw_counts))
raw_counts <- raw_counts[, -1]

positions <- read.csv( '../barcode_matching/Puck_200727_08_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes <- read.csv('../barcode_matching/Puck_200727_08_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions) <- barcodes$V1
positions <- data.frame(positions[, 2:3])
colnames(positions) = c('xcoord','ycoord')

### REPLICATE 2
setwd('/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_09/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts2 = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_09/hg19.exonic+intronic/alignment/Puck_200727_09.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts2) <- toupper(x = raw_counts2[, 1])
rownames(raw_counts2) <- make.names(toupper(x=raw_counts2[,1]), unique = TRUE)
colnames(x = raw_counts2) <- toupper(x = colnames(x = raw_counts2))
raw_counts2 <- raw_counts2[, -1]

positions2 <- read.csv( '../barcode_matching/Puck_200727_09_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes2 <- read.csv('../barcode_matching/Puck_200727_09_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions2) <- barcodes2$V1
positions2 <- data.frame(positions2[, 2:3])
colnames(positions2) = c('xcoord','ycoord')

### REPLICATE 3
setwd('/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_10/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts3 = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_10/hg19.exonic+intronic/alignment/Puck_200727_10.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts3) <- toupper(x = raw_counts3[, 1])
rownames(raw_counts3) <- make.names(toupper(x=raw_counts3[,1]), unique = TRUE)
colnames(x = raw_counts3) <- toupper(x = colnames(x = raw_counts3))
raw_counts3 <- raw_counts3[, -1]

positions3 <- read.csv( '../barcode_matching/Puck_200727_10_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes3 <- read.csv('../barcode_matching/Puck_200727_10_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions3) <- barcodes3$V1
positions3 <- data.frame(positions3[, 2:3])
colnames(positions3) = c('xcoord','ycoord')
###

# Merge seurat objects
rcc8 <- CreateSeuratObject(counts = raw_counts, project = "SlideSeq")
rcc9 <- CreateSeuratObject(counts = raw_counts2, project = "SlideSeq")
rcc10 <- CreateSeuratObject(counts = raw_counts3, project = "SlideSeq")

mydata_s <- merge(rcc8, y = c(rcc9, rcc10), add.cell.ids = c("rep8", "rep9", "rep10"), project = "RCC_post")
mydata_s[["percent.mt"]] <- PercentageFeatureSet(mydata_s, pattern = "^MT-")
mydata_s <- SCTransform(mydata_s, verbose = FALSE)
mydata_s <- RunPCA(mydata_s)
mydata_s <- RunUMAP(mydata_s, dims = 1:30)
mydata_s <- FindNeighbors(mydata_s, dims = 1:30)
mydata_s <- FindClusters(mydata_s, resolution = 0.3, verbose = FALSE)
plot5 <- DimPlot(mydata_s, reduction = "umap", label = TRUE)
plot5

mydata_s.markers <- FindAllMarkers(mydata_s, min.pct = 0.25, logfc.threshold = 0.25)
mydata_s.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top20 <- mydata_s.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Seurat::DoHeatmap(mydata_s, features = top20$gene) + NoLegend()

FeaturePlot(mydata_s, features = c("TRBC2",'TRAC','CD8A','CD3D','CD4'))

SpatialDimPlot(mydata_s, cells.highlight = CellsByIdentities(object = mydata_s, idents = c(0,
                                                                                           1,2,3,4,5,6,7,8,9,10,11,12, 13,14,15)), facet.highlight = TRUE,
               cols.highlight=c("#DE2D26", "#EBECF0"),stroke=0)

plot7 <- SpatialFeaturePlot(mydata_s, features = c('TRBC2','TRAC')) + theme(legend.position = "right")

cluster_assignments_to_save <- data.frame(rownames(mydata_s@meta.data))
cluster_assignments_to_save['cluster_id'] <- mydata_s@meta.data$seurat_clusters
colnames(cluster_assignments_to_save) <- c('barcode','cluster')
write.csv(cluster_assignments_to_save,'/Volumes/broad_thechenlab/Fei/TCR/Wu_200909/Aggregated_clusterassignments.csv',row.names=FALSE)

index <- c(paste('rep8_',colnames(raw_counts),sep=''),
           paste('rep9_',colnames(raw_counts2),sep=''),
           paste('rep10_',colnames(raw_counts3),sep=''))
replicate_status <- data.frame(c(replicate(length(colnames(raw_counts)),'rep8') ,replicate(length(colnames(raw_counts2)),'rep9') ,
                                 replicate(length(colnames(raw_counts3)),'rep10')))
replicate_status

save(mydata_s,file='/Volumes/broad_thechenlab/Fei/TCR/Wu_200909/vectorized_figures/rctd_all_pucks_merged_seurat_object.Rdata')

rownames(replicate_status) <- index
mydata_s$replicate <- replicate_status

# Next, switch the identity class of all cells to reflect replicate ID
Idents(mydata_s) <- "replicate"
DimPlot(mydata_s, reduction = "umap")

Idents(mydata_s) <- mydata_s@meta.data$seurat_clusters
DimPlot(mydata_s, reduction = "umap")

# plot each replicate seprately
umap_tx = mydata_s@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = mydata_s$replicate)

ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point(fill=NA,alpha=0.1,size=1) + 
  scale_color_manual(values=c("rep9" = "gray", "rep10"="gray","rep8" = "darkblue"))
#

#####
##### PRE-TREATMENT AGGREGATED REPLICATES
#Try aggregating all replicates

setwd('/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_12/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_12/hg19.exonic+intronic/alignment/Puck_200727_12.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts) <- toupper(x = raw_counts[, 1])
rownames(raw_counts) <- make.names(toupper(x=raw_counts[,1]), unique = TRUE)
colnames(x = raw_counts) <- toupper(x = colnames(x = raw_counts))
raw_counts <- raw_counts[, -1]


positions <- read.csv( '../barcode_matching/Puck_200727_12_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes <- read.csv('../barcode_matching/Puck_200727_12_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions) <- barcodes$V1
positions <- data.frame(positions[, 2:3])
colnames(positions) = c('xcoord','ycoord')


### REPLICATE 2
setwd('/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_13/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts2 = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-09-02_Puck_200727_13/hg19.exonic+intronic/alignment/Puck_200727_13.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts2) <- toupper(x = raw_counts2[, 1])
rownames(raw_counts2) <- make.names(toupper(x=raw_counts2[,1]), unique = TRUE)
colnames(x = raw_counts2) <- toupper(x = colnames(x = raw_counts2))
raw_counts2 <- raw_counts2[, -1]


positions2 <- read.csv( '../barcode_matching/Puck_200727_13_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes2 <- read.csv('../barcode_matching/Puck_200727_13_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions2) <- barcodes2$V1
positions2 <- data.frame(positions2[, 2:3])
colnames(positions2) = c('xcoord','ycoord')


# Merge seurat objects
rcc12 <- CreateSeuratObject(counts = raw_counts, project = "SlideSeq")
rcc13 <- CreateSeuratObject(counts = raw_counts2, project = "SlideSeq")


mydata_s <- merge(rcc12, y = c(rcc13), add.cell.ids = c("rep12", "rep13"), project = "RCC_post")
mydata_s[["percent.mt"]] <- PercentageFeatureSet(mydata_s, pattern = "^MT-")
mydata_s <- SCTransform(mydata_s, verbose = FALSE)
mydata_s <- RunPCA(mydata_s)
mydata_s <- RunUMAP(mydata_s, dims = 1:30)
mydata_s <- FindNeighbors(mydata_s, dims = 1:30)
mydata_s <- FindClusters(mydata_s, resolution = 0.3, verbose = FALSE)
plot5 <- DimPlot(mydata_s, reduction = "umap", label = TRUE)
plot5

mydata_s.markers <- FindAllMarkers(mydata_s, min.pct = 0.25, logfc.threshold = 0.25)
mydata_s.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top20 <- mydata_s.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Seurat::DoHeatmap(mydata_s, features = top20$gene) + NoLegend()

FeaturePlot(mydata_s, features = c("TRBC2",'TRAC','CD8A','CD3D','CD4'))

SpatialDimPlot(mydata_s, cells.highlight = CellsByIdentities(object = mydata_s, idents = c(0,
                                                                                           1,2,3,4,5,6,7,8,9,10,11,12, 13,14,15)), facet.highlight = TRUE,
               cols.highlight=c("#DE2D26", "#EBECF0"),stroke=0)

plot7 <- SpatialFeaturePlot(mydata_s, features = c('TRBC2','TRAC')) + theme(legend.position = "right")

cluster_assignments_to_save <- data.frame(rownames(mydata_s@meta.data))
cluster_assignments_to_save['cluster_id'] <- mydata_s@meta.data$seurat_clusters
colnames(cluster_assignments_to_save) <- c('barcode','cluster')
write.csv(cluster_assignments_to_save,'/Volumes/broad_thechenlab/Fei/TCR/Wu_200909/Aggregated_clusterassignments_pretreatment.csv',row.names=FALSE)

index <- c(paste('rep8_',colnames(raw_counts),sep=''),
           paste('rep9_',colnames(raw_counts2),sep=''),
           paste('rep10_',colnames(raw_counts3),sep=''))

index <- c(paste('rep12_',colnames(raw_counts),sep=''),
           paste('rep13_',colnames(raw_counts2),sep=''))

replicate_status <- data.frame(c(replicate(length(colnames(raw_counts)),'rep8') ,
                                 replicate(length(colnames(raw_counts2)),'rep9') ,
                                 replicate(length(colnames(raw_counts2)),'rep10') ))

replicate_status <- data.frame(c(replicate(length(colnames(raw_counts)),'rep12') ,
                                 replicate(length(colnames(raw_counts2)),'rep13') ))
replicate_status

rownames(replicate_status) <- index
mydata_s$replicate <- replicate_status
# Next, switch the identity class of all cells to reflect replicate ID
Idents(mydata_s) <- "replicate"
DimPlot(mydata_s, reduction = "umap")

Idents(mydata_s) <- mydata_s@meta.data$seurat_clusters

# plot each replicate seprately
umap_tx = mydata_s@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = mydata_s$replicate)

ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point(fill=NA,alpha=0.1,size=1) + 
  scale_color_manual(values=c("rep9" = "gray", "rep10"="gray","rep8" = "darkblue"))

#########
#######



setwd('/Volumes/broad_macosko/data/libraries/2020-10-26_Puck_201006_08/hg19.exonic+intronic/alignment')

### copied from spleen code that works
raw_counts = data.table::fread(
  input = '/Volumes/broad_macosko/data/libraries/2020-10-26_Puck_201006_08/hg19.exonic+intronic/alignment/Puck_201006_08.digital_expression.txt',
  sep = "\t",
  data.table = FALSE)
rownames(x = raw_counts) <- toupper(x = raw_counts[, 1])
rownames(raw_counts) <- make.names(toupper(x=raw_counts[,1]), unique = TRUE)
colnames(x = raw_counts) <- toupper(x = colnames(x = raw_counts))
raw_counts <- raw_counts[, -1]

mydata_s <- CreateSeuratObject(counts = raw_counts, project = "SlideSeq",assay = 'Spatial')

positions <- read.csv( '../barcode_matching/Puck_201006_08_matched_bead_locations.txt', sep = '\t', header= FALSE)
barcodes <- read.csv('../barcode_matching/Puck_201006_08_matched_bead_barcodes.txt',sep = '\t', header= FALSE)
rownames(positions) <- barcodes$V1
positions <- data.frame(positions[, 2:3])
colnames(positions) = c('xcoord','ycoord')
mydata_s[['image']] <- new(Class = 'SlideSeq',assay = 'Spatial', coordinates = positions)

mydata_s[["percent.mt"]] <- PercentageFeatureSet(mydata_s, pattern = "^mt-")
mydata_s <- SCTransform(mydata_s, assay = "Spatial", ncells = 3000, verbose = FALSE)
mydata_s <- RunPCA(mydata_s)
mydata_s <- RunUMAP(mydata_s, dims = 1:30)
mydata_s <- FindNeighbors(mydata_s, dims = 1:30)
mydata_s <- FindClusters(mydata_s, resolution = 0.3, verbose = FALSE)
plot5 <- DimPlot(mydata_s, reduction = "umap", label = TRUE)
plot6 <- SpatialDimPlot(mydata_s, stroke = 0,pt.size.factor=1.2)
plot5 + plot6

mydata_s.markers <- FindAllMarkers(mydata_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mydata_s.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top20 <- mydata_s.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Seurat::DoHeatmap(mydata_s, features = top20$gene) + NoLegend()

FeaturePlot(mydata_s, features = c("TRBC2",'TRAC','CD8A','CD3D','CD4'))

SpatialDimPlot(mydata_s, cells.highlight = CellsByIdentities(object = mydata_s, idents = c(0,
                                                                                           1,2,3,4,5,6,7,8,9,10,11,12, 13,14,15)), facet.highlight = TRUE,
               cols.highlight=c("#DE2D26", "#EBECF0"),stroke=0)

plot7 <- SpatialFeaturePlot(mydata_s, features = c('TRBC2','TRAC')) + theme(legend.position = "right")

cluster_assignments_to_save <- data.frame(rownames(mydata_s@meta.data))
cluster_assignments_to_save['cluster_id'] <- mydata_s@meta.data$seurat_clusters
colnames(cluster_assignments_to_save) <- c('barcode','cluster')
write.csv(cluster_assignments_to_save,'2020-10-26_Puck_201006_08_clusterassignments.csv',row.names=FALSE)

save(mydata_s,file='/Volumes/broad_thechenlab/Fei/TCR/Wu_200909/rctd_all_pucks_merged_seurat_object_pre_treatment.Rdata')
mydata_s <- get(load('/Volumes/broad_thechenlab/Fei/TCR/Wu_200909/vectorized_figures/rctd_all_pucks_merged_seurat_object.Rdata'))

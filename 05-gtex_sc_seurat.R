################################################################################
################################################################################
################################################################################
################################################################################
###########################  GTEx SC ANALYSIS (SEURAT) #########################

#################################### SETUP #####################################

library(tidyverse)
library(scater)
library(edgeR)
library(Seurat)
library(SeuratObject)
library(rtracklayer)

################################## LOAD DATA ###################################

load("r_outputs/03-gtex_seurat.norm.Rdata")

############################## VARIABLE FEATURES ###############################
######## Identification of highly variable features (feature selection) ########

####################### GTEX_1MCC2_5013_SM_HPJ3D (Breast) ######################

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindVariableFeatures(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
all.genes <- rownames(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  ScaleData(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, features = all.genes)

# Perform linear dimensional reduction
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <-
  RunPCA(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
         features = VariableFeatures(object = GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm))

print(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:2, reduction = "pca")

DimHeatmap(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:5, cells = 1000, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  JackStraw(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, num.replicate = 100)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  ScoreJackStraw(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:20)

JackStrawPlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:15)

# Choosing 11 PCs

# Cluster the cells
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindNeighbors(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:6)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindClusters(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, resolution = 0.5)

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  RunUMAP(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:6)

DimPlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GTEX_1MCC2_5013_SM_HPJ3D.markers <- 
  FindAllMarkers(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GTEX_1MCC2_5013_SM_HPJ3D.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

GTEX_1MCC2_5013_SM_HPJ3D.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top10
DoHeatmap(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, features = top10$gene) + NoLegend()

####################### GTEX_1CAMS_5015_SM_HPJ3C (Breast) ######################

GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  FindVariableFeatures(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  ScaleData(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, features = all.genes)

GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <-
  RunPCA(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, 
       features = VariableFeatures(object = GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm))

print(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:2, reduction = "pca")

DimHeatmap(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:5, cells = 1000, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  JackStraw(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, num.replicate = 100)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  ScoreJackStraw(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:20)

JackStrawPlot(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:15)

# Choosing 5 PCs

# Cluster the cells
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  FindNeighbors(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:4)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  FindClusters(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, resolution = 0.5)

GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- 
  RunUMAP(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, dims = 1:5)

DimPlot(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GTEX_1CAMS_5015_SM_HPJ3C.markers <- 
  FindAllMarkers(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GTEX_1CAMS_5015_SM_HPJ3C.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

GTEX_1CAMS_5015_SM_HPJ3C.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, features = top10$gene) + NoLegend()

####################### GTEX_1R9PN_5002_SM_HD2MC (Breast) ######################

GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  FindVariableFeatures(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
all.genes <- rownames(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm)
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  ScaleData(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, features = all.genes)

# Perform linear dimensional reduction
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <-
  RunPCA(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, 
         features = VariableFeatures(object = GTEX_1R9PN_5002_SM_HD2MC.seurat.norm))

print(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:2, reduction = "pca")

DimHeatmap(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:5, cells = 1000, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  JackStraw(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, num.replicate = 100)
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  ScoreJackStraw(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:20)

JackStrawPlot(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:15)

# Choosing 10 PCs

# Cluster the cells
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  FindNeighbors(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:10)
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  FindClusters(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, resolution = 0.5)

GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- 
  RunUMAP(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, dims = 1:5)

DimPlot(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GTEX_1R9PN_5002_SM_HD2MC.markers <- 
  FindAllMarkers(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GTEX_1R9PN_5002_SM_HD2MC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

GTEX_1R9PN_5002_SM_HD2MC.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, features = top10$gene) + NoLegend()


####################### GTEX_1MCC2_5013_SM_HPJ3D (Breast) ######################

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindVariableFeatures(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
all.genes <- rownames(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  ScaleData(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, features = all.genes)

# Perform linear dimensional reduction
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <-
  RunPCA(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
         features = VariableFeatures(object = GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm))

print(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:2, reduction = "pca")

DimHeatmap(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:5, cells = 1000, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  JackStraw(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, num.replicate = 100)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  ScoreJackStraw(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:20)

JackStrawPlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:15)

# Choosing 11 PCs

# Cluster the cells
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindNeighbors(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:6)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  FindClusters(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, resolution = 0.5)

GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- 
  RunUMAP(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, dims = 1:6)

DimPlot(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GTEX_1MCC2_5013_SM_HPJ3D.markers <- 
  FindAllMarkers(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GTEX_1MCC2_5013_SM_HPJ3D.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

GTEX_1MCC2_5013_SM_HPJ3D.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top10
DoHeatmap(GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm, features = top10$gene) + NoLegend()

######################## GTEX_1HSMQ_5005_SM_GKSJF (Lung) #######################

GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  FindVariableFeatures(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
pdf("plots/GTEX_1HSMQ_5005_SM_GKSJF_variable_features.pdf", height=7, width=7)
plot2
dev.off()

all.genes <- rownames(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  ScaleData(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, features = all.genes)

GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <-
  RunPCA(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, 
         features = VariableFeatures(object = GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm))

print(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm[["pca"]], dims = 1:5, nfeatures = 5)

pdf("plots/GTEX_1HSMQ_5005_SM_GKSJF_dim_loadings.pdf", height=6, width=8)
VizDimLoadings(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:2, reduction = "pca")
dev.off()

DimHeatmap(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:5, cells = 1000, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  JackStraw(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, num.replicate = 100)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  ScoreJackStraw(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:20)

JackStrawPlot(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:15)

# Choosing 8 PCs

# Cluster the cells
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  FindNeighbors(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:8)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  FindClusters(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, resolution = 0.5)

GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- 
  RunUMAP(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, dims = 1:5)

pdf("plots/GTEX_1HSMQ_5005_SM_GKSJF_umap.pdf", height=7, width=7)
DimPlot(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, reduction = "umap")
dev.off()

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(GTEX_1HSMQ_5005_SM_GKSJF.markers, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
GTEX_1HSMQ_5005_SM_GKSJF.markers <- 
  FindAllMarkers(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GTEX_1HSMQ_5005_SM_GKSJF.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

GTEX_1HSMQ_5005_SM_GKSJF.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top50

pdf("plots/GTEX_1HSMQ_5005_SM_GKSJF_heatmap.pdf", height=10, width=10)
DoHeatmap(GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, features = top50$gene) + NoLegend()
dev.off()

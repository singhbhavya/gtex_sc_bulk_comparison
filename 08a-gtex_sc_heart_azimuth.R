################################################################################
################################################################################
################################################################################
################################################################################
#############################  GTEx SC HEART (SEURAT) ##########################

#################################### SETUP #####################################

library(tidyverse)
library(Matrix)
library(scater)
library(Seurat)
library(SeuratData)
library(edgeR)
library(SeuratObject)
library(rtracklayer)
library(Azimuth)
library(patchwork)
library(cowplot)
library(data.table)

################################## LOAD DATA ###################################

load("r_outputs/03-gtex_seurat.norm.Rdata")

################################## MERGE DATA ##################################

# Merge all heart samples

heart.norm.merged <- merge(GTEX_13N11_5002_SM_H5JDV.seurat.norm, 
                         y = c(GTEX_15RIE_5015_SM_H8L6X.seurat.norm, 
                               GTEX_1ICG6_5003_SM_GHS9A.seurat.norm), 
                         add.cell.ids = c("GTEX_13N11_5002_SM_H5JDV", 
                                          "GTEX_15RIE_5015_SM_H8L6X",
                                          "GTEX_1ICG6_5003_SM_GHS9A"), 
                         project = "Heart",
                         merge.data = TRUE)

heart.norm.merged[['RNA']] <- 
  Seurat::AddMetaData(heart.norm.merged[['RNA']], 
                      GTEX_13N11_5002_SM_H5JDV.seurat.norm[['RNA']]@meta.features)

############################# DETERMINE DIMENSIONS ##############################

determine_dim <- function(sobj_norm){
  sobj_norm <- FindVariableFeatures(object = sobj_norm, assay = "RNA",
                                    selection.method = "vst", nfeatures = 2000)
  sobj_norm <- ScaleData(object = sobj_norm)
  sobj_norm <- RunPCA(object = sobj_norm)
  sobj_norm <- JackStraw(sobj_norm, num.replicate = 100)
  sobj_norm <- ScoreJackStraw(sobj_norm, dims = 1:20)
  return(JackStrawPlot(sobj_norm))
}

GTEX_13N11_5002_SM_H5JDV.seurat.p <- 
  determine_dim(GTEX_13N11_5002_SM_H5JDV.seurat.norm)+
  ggtitle("GTEX_13N11_5002_SM_H5JDV")
# All 5 PCs are significant
GTEX_15RIE_5015_SM_H8L6X.seurat.p <- 
  determine_dim(GTEX_15RIE_5015_SM_H8L6X.seurat.norm)+
  ggtitle("GTEX_15RIE_5015_SM_H8L6X")
# All 5 PCs are significant
GTEX_1ICG6_5003_SM_GHS9A.seurat.p <- 
  determine_dim(GTEX_1ICG6_5003_SM_GHS9A.seurat.norm)+
  ggtitle("GTEX_1ICG6_5003_SM_GHS9A")
# All 5 PCs are significant 


############################### STANDARD SEURAT ################################

standard_seurat <- function(sobj_norm){
  sobj_norm <- FindVariableFeatures(object = sobj_norm, 
                                    selection.method = "vst", 
                                    nfeatures = 5000)
  sobj_norm <- ScaleData(object = sobj_norm)
  sobj_norm <- RunPCA(object = sobj_norm)
  sobj_norm <- FindNeighbors(object = sobj_norm)
  sobj_norm <- FindClusters(object = sobj_norm)
  sobj_norm <- RunTSNE(object = sobj_norm)
  sobj_norm <- RunUMAP(object = sobj_norm, dims=1:15)
  return(DimPlot(object = sobj_norm, reduction = 'umap'))
}

GTEX_13N11_5002_SM_H5JDV.seurat.dim <- 
  standard_seurat(GTEX_13N11_5002_SM_H5JDV.seurat.norm)+
  ggtitle("GTEX_13N11_5002_SM_H5JDV")
GTEX_15RIE_5015_SM_H8L6X.seurat.dim <- 
  standard_seurat(GTEX_15RIE_5015_SM_H8L6X.seurat.norm)+
  ggtitle("GTEX_15RIE_5015_SM_H8L6X")
GTEX_1ICG6_5003_SM_GHS9A.dim <- standard_seurat(GTEX_1ICG6_5003_SM_GHS9A.seurat.norm)+
  ggtitle("GTEX_1ICG6_5003_SM_GHS9A")

pdf("plots/08-individual_heart_umaps.pdf", height=5, width=15)
plot_grid(GTEX_13N11_5002_SM_H5JDV.seurat.dim, GTEX_15RIE_5015_SM_H8L6X.seurat.dim, 
          GTEX_1ICG6_5003_SM_GHS9A.dim, nrow = 1)
dev.off()

heart.norm.merged.dim <- standard_seurat(heart.norm.merged)+
  ggtitle("All Heart")

pdf("plots/08-combined_heart_norm_merged_umap.pdf", height=5, width=5)
heart.norm.merged.dim
dev.off()


############################### AZIMUTH MAPPING ################################

mapped <- Azimuth::RunAzimuth(
  heart.norm.merged,
  reference = "heartref",
  do.adt = TRUE,
)

save(mapped, file = "r_outputs/08-azimuth_mapping_heart.RData")
load("r_outputs/08-azimuth_mapping_heart.RData")

results <- list()
if('impADT' %in% Assays(mapped)) {
  results$impADT <- mapped[['impADT']]
}
if('ref.umap' %in% Reductions(mapped)) {
  results$umap <- mapped[['ref.umap']]
}

results$pred.df <- mapped@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  dplyr::select(
    cell,
    dplyr::matches('predicted.celltype.l\\d$'),
    dplyr::matches('predicted.celltype.l\\d.score$'),
    mapping.score
  ) %>% as.data.frame

save(results, file = "r_outputs/08-azimuth_result_heart.RData")
saveRDS(results, file = "r_outputs/08-azimuth_result_heart.Rds")
load("r_outputs/08-azimuth_result_heart.RData")

########################### ADD AZIMUTH TO METADATA ############################

heart.norm.merged <- Seurat::AddAzimuthResults(heart.norm.merged, 
                                             file="r_outputs/08-azimuth_result_heart.Rds")

# Examine metadata
heart.norm.merged@meta.data %>% head

# the predicted celltype, score, and mapping.score columns are added to the
# metadata but are all NA. Add it manually here:
heart.norm.merged@meta.data$predicted.celltype.l1 <- results$pred.df$predicted.celltype.l1
heart.norm.merged@meta.data$predicted.celltype.l2 <- results$pred.df$predicted.celltype.l2
heart.norm.merged@meta.data$predicted.celltype.l1.score <- results$pred.df$predicted.celltype.l1.score
heart.norm.merged@meta.data$predicted.celltype.l2.score <- results$pred.df$predicted.celltype.l2.score
heart.norm.merged@meta.data$mapping.score <- results$pred.df$mapping.score

# See that values are added
heart.norm.merged@meta.data %>% head

saveRDS(heart.norm.merged, file = 'r_outputs/08-GCB.heart.norm.merged.mapped.Rds')

############################### FIND L1 MARKERS ################################

#

Idents(object=heart.norm.merged) <- "predicted.celltype.l1"

l1.markers <- FindAllMarkers(
  heart.norm.merged,
  test.use = 'wilcox',    # default = 'wilcox'
  only.pos = FALSE,        # default = FALSE
  min.pct = 0.1,         # default = 0.1
  logfc.threshold = 0.25, # default = 0.25
  return.thresh = 0.01    # default = 0.01
)

# Add feature metadata
fmeta <- heart.norm.merged[['RNA']]@meta.features %>%
  tibble::rownames_to_column() %>%
  select(rowname, symbol, feattype, te_class, te_family)

stopifnot(all(l1.markers$gene %in% fmeta$rowname))
orig.rownames <- rownames(l1.markers)
l1.markers <- dplyr::left_join(l1.markers, fmeta, by=c('gene' = 'rowname'))
rownames(l1.markers) <- orig.rownames
l1.markers[is.na(l1.markers)] <- ''
rm(orig.rownames, fmeta)

saveRDS(l1.markers, file = 'r_outputs/08-heart.norm.merged.l1.markers.rds')

############################### FIND L2 MARKERS ################################

Idents(object=heart.norm.merged) <- "predicted.celltype.l2"

l2.markers <- FindAllMarkers(
  heart.norm.merged,
  test.use = 'wilcox',    # default = 'wilcox'
  only.pos = FALSE,        # default = FALSE
  min.pct = 0.1,         # default = 0.1
  logfc.threshold = 0.25, # default = 0.25
  return.thresh = 0.01    # default = 0.01
)

# Add feature metadata
fmeta <- heart.norm.merged[['RNA']]@meta.features %>%
  tibble::rownames_to_column() %>%
  select(rowname, symbol, feattype, te_class, te_family)

stopifnot(all(l2.markers$gene %in% fmeta$rowname))
orig.rownames <- rownames(l2.markers)
l2.markers <- dplyr::left_join(l2.markers, fmeta, by=c('gene' = 'rowname'))
rownames(l2.markers) <- orig.rownames
l2.markers[is.na(l2.markers)] <- ''
rm(orig.rownames, fmeta)

saveRDS(l2.markers, file = 'r_outputs/08-heart.norm.merged.l2.markers.rds')

################################# PLOT UMAPs ###################################

umap1_lim <- range(heart.norm.merged[['umap.proj']]@cell.embeddings[,1])
umap2_lim <- range(heart.norm.merged[['umap.proj']]@cell.embeddings[,2])

pdf('plots/08-heart.norm.merged.mapped.l1.pdf', width=10, height=10)
nclust <- length(unique(heart.norm.merged@meta.data$predicted.celltype.l1))
p <- DimPlot(heart.norm.merged,
             repel=TRUE,
             group.by = "predicted.celltype.l1",
             label = TRUE,
             label.size = 9,
             cols=Seurat::DiscretePalette(nclust, 'polychrome'),
             raster = FALSE,
             raster.dpi = c(1200, 1200)
)
p + xlim(umap1_lim) + ylim(umap2_lim) + NoLegend() + theme(plot.title=element_blank())
dev.off()

nclust <- length(unique(heart.norm.merged@meta.data$predicted.celltype.l2))
p_l2 <- DimPlot(heart.norm.merged,
                repel = TRUE,
                label.size = 4,
                pt.size = 1,
                group.by = "predicted.celltype.l2",
                label = TRUE,
                label.box = FALSE,
                cols=c(Seurat::DiscretePalette(32, c('glasbey')),
                       Seurat::DiscretePalette(36, c('polychrome'))),
                raster = FALSE) 

pdf('plots/08-heart.norm.merged.mapped.l2_nolab.pdf', width=12, height=8)
p_l2 + theme_cowplot() +
  xlim(umap1_lim) + ylim(umap2_lim) + NoLegend() +  theme(plot.title=element_blank()) 
dev.off()

pdf('plots/08-heart.norm.merged.mapped.l2_lab.pdf', width=15, height=8) 
p_l2 + xlim(umap1_lim) + ylim(umap2_lim) + theme_cowplot() + 
  theme(plot.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2) ) 
dev.off()

############################## L1 FEATURE PLOTS ################################
# condition_pct.1: percentage of cells where the gene is detected in the cluster 
# for condition
# condition_pct.2: percentage of cells where the gene is detected on average in 
# the other clusters for condition


herv_markers_l1 <- l1.markers[l1.markers$te_class.x=='LTR',]

nrow(herv_markers_l1) # There are 236 HERV markers.
length(unique(herv_markers_l1$gene)) # There are 181 unique HERVs.

herv_markers_l2 <- l2.markers[l2.markers$te_class=='LTR',]

nrow(herv_markers_l2) # There are 286 HERV markers.
length(unique(herv_markers_l2$gene)) # There are 201 unique HERVs.

############################## L2 FEATURE PLOTS ################################

fpcols <- c('#eeeeeeFF', viridis::viridis(6))

## Plasmablast 
Smooth_muscle_all <- herv_markers_l2[herv_markers_l2$cluster == 'Smooth Muscle', 'gene']
lapply(Smooth_muscle_all, function(x) herv_markers_l2[herv_markers_l2$gene == x,])
b_x <- c('HERVH-1p36.', 'LTR46-1p36.13', 'MER4-1p36.13', 'PRIMA4-12p11.21b', 
         'LTR46-Xq11.1', 'MER101-Xq23a', 'HERVH-Xq13.2c', 'HERV3-8q11.23', 
         'HARLEQUIN-Xq23b')
pdf('plots/08-heart_combined_l2_herv_markers_smooth_muscle.pdf', width=9, height=9)
p <- FeaturePlot(heart.norm.merged, b_x, cols=fpcols, ncol=2, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & theme(plot.title=element_text(size=10)) 
dev.off()

############################ ALL HEART UNSUPERVISED ############################


heart.norm.merged <- 
  FindVariableFeatures(heart.norm.merged, 
                       selection.method = "vst",
                       nfeatures = 2000)

# Identify the 10 most highly variable genes
top50_heart_merged <- head(VariableFeatures(heart.norm.merged), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(heart.norm.merged)
plot2 <- LabelPoints(plot = plot1, points = top50_heart_merged, repel = TRUE)
pdf("plots/heart_merged_variable_features.pdf", height=7, width=7)
plot2
dev.off()

all.genes <- rownames(heart.norm.merged)
heart.norm.merged <- 
  ScaleData(heart.norm.merged, features = all.genes)

heart.norm.merged <-
  RunPCA(heart.norm.merged, 
         features = VariableFeatures(object = heart.norm.merged))

print(heart.norm.merged[["pca"]], dims = 1:5, nfeatures = 5)

pdf("plots/heart_merged_dim_loadings.pdf", height=6, width=8)
VizDimLoadings(heart.norm.merged, dims = 1:2, reduction = "pca")
dev.off()

pdf("plots/heart_merged_dim_heatmap.pdf", height=6, width=9)
DimHeatmap(heart.norm.merged, dims = 1:5, cells = 1000, balanced = TRUE)
dev.off()

# Determine the ‘dimensionality’ of the dataset

heart.norm.merged <- 
  JackStraw(heart.norm.merged, num.replicate = 100)
heart.norm.merged <- 
  ScoreJackStraw(heart.norm.merged, dims = 1:20)

JackStrawPlot(heart.norm.merged, dims = 1:15)

# Choosing 15 PCs

# Cluster the cells
heart.norm.merged <- 
  FindNeighbors(heart.norm.merged, dims = 1:15)
heart.norm.merged <- 
  FindClusters(heart.norm.merged, resolution = 0.5)

heart.norm.merged <- 
  RunUMAP(heart.norm.merged, dims = 1:5)

pdf("plots/heart.norm.merged_umap.pdf", height=5, width=5)
DimPlot(heart.norm.merged, reduction = "umap")
dev.off()

# Finding differentially expressed features (cluster biomarkers)

cluster2.markers <- FindMarkers(heart.norm.merged, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
heart.norm.merged.markers <- 
  FindAllMarkers(heart.norm.merged, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart.norm.merged.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

heart.norm.merged.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top50

pdf("plots/heart.norm.merged.heatmap.pdf", height=16, width=22)
DoHeatmap(heart.norm.merged, features = top50$gene) 
dev.off()

pdf('plots/08-heart_combined_unsupervised_LTR46-Xq11.1.pdf', height=5, width=5)
FeaturePlot(heart.norm.merged, c("LTR46-Xq11.1"), cols=fpcols, ncol=2, raster=FALSE) + NoLegend()
dev.off()
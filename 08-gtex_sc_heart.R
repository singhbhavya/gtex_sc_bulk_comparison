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

saveRDS(GCB.norm.merged, file = 'r_outputs/04-GCB.norm.merged.mapped.Rds')

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

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
load("r_outputs/02a-GTEx_8_snRNAseq_atlas.seurat.RData")

# Remove unwanted Seurat objects
remove(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm, GTEX_13N11_5030_SM_H5JDW.seurat.norm,
       GTEX_144GM_5010_SM_HD2M8.seurat.norm, GTEX_145ME_5005_SM_H8L6T.seurat.norm,
       GTEX_145ME_5018_SM_G8XQB.seurat.norm, GTEX_15RIE_5021_SM_H8L6Y.seurat.norm,
       GTEX_15EOM_5003_SM_G64IH.seurat.norm, GTEX_1I1GU_5006_SM_G8XQC.seurat.norm,
       GTEX_16BQI_5013_SM_H8SUW.seurat.norm, GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm,
       GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm,
       GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm,
       GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm, GTEX_1ICG6_5014_SM_GHS9D.seurat.norm,
       GTEX_15CHR_5014_SM_H5JDU.seurat.norm, GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm,
       GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, GTEX_15CHR_5005_SM_H5JDT.seurat.norm,
       GTEX_15SB6_5008_SM_H8L72.seurat.norm, GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm)

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

################################# ADD METADATA #################################

all_gtex_metadata <- GTEx_8_snRNAseq_atlas.seurat@meta.data
gtex_metadata <- GTEx_8_snRNAseq_atlas.seurat@meta.data
gtex_metadata <- gtex_metadata[, c("tissue", "Sample.ID", "barcode", "Broad.cell.type",
                                   "Granular.cell.type", "prep")]
gtex_metadata <- gtex_metadata[gtex_metadata$tissue == "heart" & 
                                 gtex_metadata$prep == "CST",]
gtex_metadata$Sample.ID <- gsub("-", "_", gtex_metadata$Sample.ID)
gtex_metadata$Cell <- paste0(gtex_metadata$Sample.ID, "_", gtex_metadata$barcode)
gtex_metadata <- gtex_metadata[!duplicated(gtex_metadata), ]
row.names(gtex_metadata) <- gtex_metadata$Cell

remove(GTEx_8_snRNAseq_atlas.seurat)

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

pdf("plots/08b-individual_heart_umaps.pdf", height=5, width=15)
plot_grid(GTEX_13N11_5002_SM_H5JDV.seurat.dim, GTEX_15RIE_5015_SM_H8L6X.seurat.dim, 
          GTEX_1ICG6_5003_SM_GHS9A.dim, nrow = 1)
dev.off()

############################ ALL HEART UNSUPERVISED ############################

heart.norm.merged <- 
  FindVariableFeatures(heart.norm.merged, 
                       selection.method = "vst",
                       nfeatures = 5000)

# Identify the 10 most highly variable genes
top50_heart_merged <- head(VariableFeatures(heart.norm.merged), 50)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(heart.norm.merged, 
                             cols = c("black", "red"),
                             pt.size = 1,
                             log = NULL,
                             selection.method = NULL,
                             assay = NULL,
                             raster = NULL,
                             raster.dpi = c(512, 512))

plot2 <- LabelPoints(plot = plot1, 
                     points = top50_heart_merged, 
                     repel = TRUE)
pdf("plots/08b-heart_merged_variable_features.pdf", height=7, width=7)
plot2
dev.off()

all.genes <- rownames(heart.norm.merged)
heart.norm.merged <- 
  ScaleData(heart.norm.merged, features = all.genes)

heart.norm.merged <-
  RunPCA(heart.norm.merged, 
         features = VariableFeatures(object = heart.norm.merged))

print(heart.norm.merged[["pca"]], dims = 1:5, nfeatures = 5)

pdf("plots/08b-heart_merged_dim_loadings.pdf", height=6, width=8)
VizDimLoadings(heart.norm.merged, dims = 1:2, reduction = "pca")
dev.off()

pdf("plots/08b-heart_merged_dim_heatmap.pdf", height=6, width=9)
DimHeatmap(heart.norm.merged, dims = 1:5, cells = 1000, balanced = TRUE)
dev.off()

# Determine the ‘dimensionality’ of the dataset
heart.norm.merged <- 
  JackStraw(heart.norm.merged, num.replicate = 100)
heart.norm.merged <- 
  ScoreJackStraw(heart.norm.merged, dims = 1:20)
JackStrawPlot(heart.norm.merged, dims = 1:15)
# keep 11 PCs

################################# CLUSTER CELLS ################################

# Cluster the cells
heart.norm.merged <- 
  FindNeighbors(heart.norm.merged, dims = 1:11)
heart.norm.merged <- 
  FindClusters(heart.norm.merged, resolution = 0.5)

heart.norm.merged <- 
  RunUMAP(heart.norm.merged, dims = 1:11)

pdf("plots/08b-heart.norm.merged_umap.pdf", height=5, width=5)
DimPlot(heart.norm.merged, reduction = "umap", 
        cols=Seurat::DiscretePalette(13, 'glasbey')[1:13], pt.size = 0.3)
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
  top_n(n = 13, wt = avg_log2FC) -> top50

pdf("plots/08b-heart.norm.merged.heatmap.pdf", height=20, width=22)
DoHeatmap(heart.norm.merged, features = top50$gene) 
dev.off()

################################# FEATURE PLOTS ################################

fpcols <- c('#eeeeeeFF', viridis::viridis(6))

# Clusters 0, 2, 3, 4, 7, 9: Cardiomyocytes
pdf("plots/08b-heart_cluster0_3_4.pdf", height=15, width=18)
FeaturePlot(heart.norm.merged, c("TTN",  "MYBPC3", "TNNT2", "RYR2", "PLN", "SLC8A1", 
                                 "MYH7","MYL2", "L1FLI-10q25.1b", "L1FLnI-1q43a", 
                                 "HERV3-8q11.23", "L1FLnI-5q12.3v"), 
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

# Clusters 1 & 6: Myofibroblast
pdf("plots/08b-heart_cluster1.pdf", height=15, width=18)
FeaturePlot(heart.norm.merged, c("COL3A1", "COL1A2", "COL1A1", "CCDC80", "BGN",
                                 "SERPINE2", "MGP", "POSTN", 
                                 "L1FLnI-Xq21.1yc", "L1FLnI-13q31.1u", "L1FLnI-15q22.2e", 
                                 "L1FLnI-17q24.2c"), 
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 2: Cytoplasmic cardiomyocytes?
pdf("plots/08b-heart_cluster2.pdf", height=6, width=9)
FeaturePlot(heart.norm.merged, c("MYH7", "CKM", "NDUFA4", "TNNC1","MYL3", "S100A1"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 3: LncRNA and TE-rich cardiomyocytes
pdf("plots/08b-heart_cluster3.pdf", height=15, width=18)
FeaturePlot(heart.norm.merged, c("DONSON", "SORBS2", "TTN", "PLIN5", "SNHG14",
                                 "XIST", "L1FLnI-5q33.2i", "L1FLnI-12p13.33d",
                                 "L1FLI-3q13.13", "L1FLnI-14q12cb", "L1FLnI-4q12l",
                                 "HERVL-3p21.31a"), 
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 5: Immune (DC/macrophage)

# Cluster 8: vascular endothelial cells

# Cluster 10: Pericyte / SMC

# Cluster 11: Adipocyte

# Cluster 12: lymphatic endothelial cells


################################## ADD IDENTS ##################################

heart.norm.merged.gtex <- heart.norm.merged
heart.norm.merged.gtex <- AddMetaData(object = heart.norm.merged.gtex, 
                                         metadata = gtex_metadata)
heart.norm.merged.gtex[["old.ident"]] <- Idents(object = heart.norm.merged.gtex)

# Set identity classes
Idents(object = heart.norm.merged.gtex) <- 
  heart.norm.merged.gtex@meta.data$Broad.cell.type

pdf("plots/09-prostate.norm.merged_umap_original_idents.pdf", height=5, width=8)
DimPlot(heart.norm.merged.gtex, reduction = "umap", 
        cols=Seurat::DiscretePalette(12, 'glasbey'))
dev.off()

############################# ADD GRANULAR IDENTS ##############################

heart.norm.merged.azimuth <- heart.norm.merged
heart.norm.merged.azimuth <- AddMetaData(object = heart.norm.merged.gtex, 
                                      metadata = gtex_metadata)
heart.norm.merged.gtex[["old.ident"]] <- Idents(object = heart.norm.merged.gtex)

# Set identity classes
Idents(object = heart.norm.merged.gtex) <- 
  heart.norm.merged.gtex@meta.data$Broad.cell.type

pdf("plots/08b-heart.norm.merged_umap_original_idents.pdf", height=5, width=8)
DimPlot(heart.norm.merged.gtex, reduction = "umap", pt.size = 0.3,
        cols=Seurat::DiscretePalette(12, 'glasbey'))
dev.off()

# Set identity classes
Idents(object = heart.norm.merged.gtex) <- 
  heart.norm.merged.gtex@meta.data$Granular.cell.type

pdf("plots/08b-heart.norm.merged_umap_original_idents_granular.pdf", height=5, width=12)
DimPlot(heart.norm.merged.gtex, reduction = "umap", pt.size = 0.3,
        cols=Seurat::DiscretePalette(20, 'glasbey'))
dev.off()


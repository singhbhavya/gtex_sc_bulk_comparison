################################################################################
################################################################################
################################################################################
################################################################################
###########################  GTEx SC PROSTATE (SEURAT) #########################

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
load("r_outputs/04d-prostate_unique_shared_te_cpm_0.5.RData")
load("r_outputs/01-counts.Rdata")
load("r_outputs/02a-GTEx_8_snRNAseq_atlas.seurat.RData")
load("r_outputs/03-scgtex_seurat_counts.Rdata")

# Remove unwanted Seurat objects
remove(GTEX_13N11_5002_SM_H5JDV.seurat.norm, GTEX_13N11_5030_SM_H5JDW.seurat.norm,
       GTEX_144GM_5010_SM_HD2M8.seurat.norm, GTEX_145ME_5005_SM_H8L6T.seurat.norm,
       GTEX_145ME_5018_SM_G8XQB.seurat.norm, GTEX_15RIE_5021_SM_H8L6Y.seurat.norm,
       GTEX_15EOM_5003_SM_G64IH.seurat.norm, GTEX_15RIE_5015_SM_H8L6X.seurat.norm, 
       GTEX_16BQI_5013_SM_H8SUW.seurat.norm, GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm,
       GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm,
       GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm,
       GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm, GTEX_1ICG6_5014_SM_GHS9D.seurat.norm,
       GTEX_1ICG6_5003_SM_GHS9A.seurat.norm, GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm,
       GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, GTEX_15CHR_5005_SM_H5JDT.seurat.norm,
       GTEX_15SB6_5008_SM_H8L72.seurat.norm)

################################## MERGE DATA ##################################

# Merge all prostate samples

prostate.norm.merged <- merge(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm, 
                           y = c(GTEX_15CHR_5014_SM_H5JDU.seurat.norm, 
                                 GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm,
                                 GTEX_1I1GU_5006_SM_G8XQC.seurat.norm), 
                           add.cell.ids = c("GTEX_12BJ1_5007_SM_H8L6U", 
                                            "GTEX_15CHR_5014_SM_H5JDU",
                                            "GTEX_1HSMQ_5014_SM_GKSJI",
                                            "GTEX_1I1GU_5006_SM_G8XQC"), 
                           project = "Prostate",
                           merge.data = TRUE)

prostate.norm.merged[['RNA']] <- 
  Seurat::AddMetaData(prostate.norm.merged[['RNA']], 
                      GTEX_12BJ1_5007_SM_H8L6U.seurat.norm[['RNA']]@meta.features)

################################# ADD METADATA #################################

all_gtex_metadata <- GTEx_8_snRNAseq_atlas.seurat@meta.data
gtex_metadata <- GTEx_8_snRNAseq_atlas.seurat@meta.data
gtex_metadata <- gtex_metadata[, c("tissue", "Sample.ID", "barcode", "Broad.cell.type",
                                   "Granular.cell.type", "prep")]
gtex_metadata <- gtex_metadata[gtex_metadata$tissue == "prostate" & 
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

GTEX_12BJ1_5007_SM_H8L6U.seurat.p <- 
  determine_dim(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm)+
  ggtitle("GTEX_12BJ1_5007_SM_H8L6U")
# All 5 PCs are significant
GTEX_15CHR_5014_SM_H5JDU.seurat.p <- 
  determine_dim(GTEX_15CHR_5014_SM_H5JDU.seurat.norm)+
  ggtitle("GTEX_15CHR_5014_SM_H5JDU")
# First 4 PCs are significant 
GTEX_1HSMQ_5014_SM_GKSJI.seurat.p <- 
  determine_dim(GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm)+
  ggtitle("GTEX_1HSMQ_5014_SM_GKSJI")
# All 5 PCs are significant 
GTEX_1I1GU_5006_SM_G8XQC.seurat.p <- 
  determine_dim(GTEX_1I1GU_5006_SM_G8XQC.seurat.norm)+
  ggtitle("GTEX_1I1GU_5006_SM_G8XQC")
# All 5 PCs are significant 

prostate.norm.merged.p <-
  determine_dim(prostate.norm.merged)+
  ggtitle("All Prostate")
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

GTEX_12BJ1_5007_SM_H8L6U.seurat.dim <- 
  standard_seurat(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm)+
  ggtitle("GTEX_12BJ1_5007_SM_H8L6U")
GTEX_15CHR_5014_SM_H5JDU.seurat.dim <- 
  standard_seurat(GTEX_15CHR_5014_SM_H5JDU.seurat.norm)+
  ggtitle("GTEX_15CHR_5014_SM_H5JDU")
GTEX_1HSMQ_5014_SM_GKSJI.seurat.dim <- 
  standard_seurat(GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm)+
  ggtitle("GTEX_1HSMQ_5014_SM_GKSJI")
GTEX_1I1GU_5006_SM_G8XQC.seurat.dim <- 
  standard_seurat(GTEX_1I1GU_5006_SM_G8XQC.seurat.norm)+
  ggtitle("GTEX_1I1GU_5006_SM_G8XQC")

pdf("plots/09-individual_prostate_umaps.pdf", height=10, width=10)
plot_grid(GTEX_12BJ1_5007_SM_H8L6U.seurat.dim, GTEX_15CHR_5014_SM_H5JDU.seurat.dim, 
          GTEX_1HSMQ_5014_SM_GKSJI.seurat.dim, GTEX_1I1GU_5006_SM_G8XQC.seurat.dim,
          nrow = 2)
dev.off()

########################### ALL PROSTATE UNSUPERVISED ##########################

prostate.norm.merged <- 
  FindVariableFeatures(prostate.norm.merged, 
                       selection.method = "vst",
                       nfeatures = 5000)

# Identify the 10 most highly variable genes
top50_prostate_merged <- head(VariableFeatures(prostate.norm.merged), 50)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(prostate.norm.merged, 
                             cols = c("black", "red"),
                             pt.size = 1,
                             log = NULL,
                             selection.method = NULL,
                             assay = NULL,
                             raster = NULL,
                             raster.dpi = c(512, 512))

plot2 <- LabelPoints(plot = plot1, 
                     points = top50_prostate_merged, 
                     repel = TRUE)
pdf("plots/09-prostate_merged_variable_features.pdf", height=7, width=7)
plot2
dev.off()

all.genes <- rownames(prostate.norm.merged)
prostate.norm.merged <- 
  ScaleData(prostate.norm.merged, features = all.genes)

prostate.norm.merged <-
  RunPCA(prostate.norm.merged, 
         features = VariableFeatures(object = prostate.norm.merged))

print(prostate.norm.merged[["pca"]], dims = 1:5, nfeatures = 5)

pdf("plots/09-prostate_merged_dim_loadings.pdf", height=6, width=8)
VizDimLoadings(prostate.norm.merged, dims = 1:2, reduction = "pca")
dev.off()

pdf("plots/09-prostate_merged_dim_heatmap.pdf", height=6, width=9)
DimHeatmap(prostate.norm.merged, dims = 1:5, cells = 1000, balanced = TRUE)
dev.off()

# Determine the ‘dimensionality’ of the dataset
prostate.norm.merged <- 
  JackStraw(prostate.norm.merged, num.replicate = 100)
prostate.norm.merged <- 
  ScoreJackStraw(prostate.norm.merged, dims = 1:20)
JackStrawPlot(prostate.norm.merged, dims = 1:15)
# keep 11 PCs

################################# CLUSTER CELLS ################################

# Cluster the cells
prostate.norm.merged <- 
  FindNeighbors(prostate.norm.merged, dims = 1:11)
prostate.norm.merged <- 
  FindClusters(prostate.norm.merged, resolution = 0.5)

prostate.norm.merged <- 
  RunUMAP(prostate.norm.merged, dims = 1:11)

pdf("plots/09-prostate.norm.merged_umap.pdf", height=5, width=5)
DimPlot(prostate.norm.merged, reduction = "umap", 
        cols=Seurat::DiscretePalette(10, 'polychrome')[1:10])
dev.off()

# Finding differentially expressed features (cluster biomarkers)
cluster2.markers <- FindMarkers(prostate.norm.merged, 
                                ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# Finding differentially expressed features (cluster biomarkers)
cluster6.markers <- FindMarkers(prostate.norm.merged, 
                                ident.1 = "Luminal Neuroendocrine epithelia", 
                                min.pct = 0.25)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
prostate.norm.merged.markers <- 
  FindAllMarkers(prostate.norm.merged, 
                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
prostate.norm.merged.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

prostate.norm.merged.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top50

pdf("plots/09-prostate.norm.merged.heatmap.pdf", height=16, width=22)
DoHeatmap(prostate.norm.merged, features = top50$gene) 
dev.off()

################################# FEATURE PLOTS ################################

fpcols <- c('#eeeeeeFF', viridis::viridis(6))


# Cluster 5: Myocyte (smooth muscle)
pdf("plots/09-prostate_cluster5.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("TAGLN", "MYL9", "TPM2", "ACTA2", "CNN1",
                                    "ACTG2"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 

dev.off()

# Cluster 9: Immune (DC/Macrophage)
pdf("plots/09-prostate_cluster9.pdf", height=15, width=15)
FeaturePlot(prostate.norm.merged, c("RGS1", "C1QA", "C1QB", "C1QC", "CD163", 
                                    "HAVCR2", "CXCR4",
                                    "ERV316A3-2q22.2b", "L1FLnI-12p13.31c", ), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 7: Endothelia (vascular/lymphatic)
pdf("plots/09-prostate_cluster7.pdf", height=15, width=15)
FeaturePlot(prostate.norm.merged, c("IFI27", "ACKR1", "CLDN5", "PECAM1", "AQP1", "ENG",
                                    "PLVAP", "VWF", "ADGRL4"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 0: Luminal epithelia
pdf("plots/09-prostate_cluster0.pdf", height=15, width=15)
FeaturePlot(prostate.norm.merged, c("PLA2G2A", "CD24", "MSMB", "KLK3", "AZGP1",
                                    "P4HB","BASP1","MER4-2q21.1b", "MER4-14q11.2c"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 6: Luminal Neuroendocrine epithelia
pdf("plots/09-prostate_cluster6.pdf", height=15, width=19)
FeaturePlot(prostate.norm.merged, c("FASN", "CPE", "CEBPD", "HP",
                                    "SLC29A4", "ASRGL1", "PCSK1N", "SEC11C",
                                    "SERPINA3", "HML4-8q24.3", "HERVH-15q26.1c"), 
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 3: Club epithelia + Hillock epithelia
pdf("plots/09-prostate_cluster3.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("KRT13","KRT19","SCGB3A1", 
                                    "OLFM4", "PIGR", "SCUBE2"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 4: Basal epithelia
pdf("plots/09-prostate_cluster4.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("CHCHD2", "CEBPD", "BTF3", "HINT1", "JUNB",
                                    "HERVH-Xq13.2c"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# Cluster 8: Fibroblast
pdf("plots/09-prostate_cluster8.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("SCN7A", "LUM", "DCN", "CCDC80", "C7", "MMP2"), 
            cols=fpcols, ncol=3, raster=TRUE,pt.size = 2) 
dev.off()

# Cluster 1: HERV-rich?
pdf("plots/09-prostate_cluster1.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("HERVH-Xq13.2c", "LTR46-Xq11.1", 
                                    "L1FLnI-2p13.2a", "HERV3-19q13.42b",
                                    "MER101-16p12.2a", "HERVIP10FH-Yq11.222b"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2)
dev.off()

# Cluster 2: LINE-RICH luminal epithelium
pdf("plots/09-prostate_cluster2.pdf", height=15, width=15)
FeaturePlot(prostate.norm.merged, c("KLK2", "KLK3", "POTEM",
                                    "L1ORF2-18q11.2a", "L1FLnI-5q32ba", 
                                    "L1FLnI-4q13.1va", "L1FLI-Xp21.1a",
                                    "L1FLnI-6q24.3h", "L1FLnI-4q28.1ba"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

# 22q11.23
pdf("plots/09-prostate_22q11.23.pdf", height=10, width=15)
FeaturePlot(prostate.norm.merged, c("HML2-22q11.23", "ERVLE-22q11.23b", 
                                    "HERVL40-22q11.23", "PRIMA41-22q11.23", 
                                    "HML2-1q21.3", "HML2-3q12.3"), 
            cols=fpcols, ncol=3, raster=TRUE, pt.size = 2) 
dev.off()

pdf("plots/09-prostate_HML2-22q11.23.pdf", height=5, width=10)
FeaturePlot(prostate.norm.merged, c("HML2-22q11.23", "ASRGL1"), 
            cols=fpcols, ncol=2, raster=TRUE, pt.size = 2) 

FeaturePlot(prostate.norm.merged, c("HML2-11q12.3", "ASRGL1"), 
            cols=fpcols, ncol=2, raster=TRUE, pt.size = 2) 
dev.off()

# Name clusters
new.cluster.ids <- c("Luminal epithelia", "HERV-rich luminal epithelia", "LINE-rich luminal epithelia", 
                     "Club/Hillock/basal epithelia", "Basal epithelia", 
                     "Myocyte (smooth muscle)", "Luminal epithelia sub-type unknown", "Endothelia",
                     "Fibroblast", "Immune (DC/Macrophage)")
names(new.cluster.ids) <- levels(prostate.norm.merged)
prostate.norm.merged <- RenameIdents(prostate.norm.merged, new.cluster.ids)

# Re-run UMAP and dim plot with new cluster numbers
pdf("plots/09-prostate.norm.merged_umap_labelled.pdf", height=5, width=8)
DimPlot(prostate.norm.merged, reduction = "umap", 
        cols=Seurat::DiscretePalette(10, 'polychrome')[1:10])
dev.off()

pdf("plots/09-prostate.norm.merged.heatmap_labelled.pdf", height=20, width=22)
DoHeatmap(prostate.norm.merged, features = top50$gene) 
dev.off()


# Epithelial cell markers
pdf("plots/09-prostate_epithelial_markers.pdf", height=15, width=18)
FeaturePlot(prostate.norm.merged, c("MUC1", "KRT14", "CLDN1", "OCLN", "EPCAM", 
                                    "KLK3", "KLK2", "KLK1", "CDH1", "TXNDC2", 
                                    "MTRNR2L1"), 
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 3) 
dev.off()

################################## ADD IDENTS ##################################

prostate.norm.merged.gtex <- prostate.norm.merged
prostate.norm.merged.gtex <- AddMetaData(object = prostate.norm.merged.gtex, 
                                    metadata = gtex_metadata)
prostate.norm.merged.gtex[["old.ident"]] <- Idents(object = prostate.norm.merged.gtex)

# Set identity classes
Idents(object = prostate.norm.merged.gtex) <- 
  prostate.norm.merged.gtex@meta.data$Broad.cell.type

pdf("plots/09-prostate.norm.merged_umap_original_idents.pdf", height=5, width=8)
DimPlot(prostate.norm.merged.gtex, reduction = "umap", 
        cols=c("#66013d", "#B00068", "#FE00FA", "#f288f0", "#855284", "#5A5156",
               "#1CFFCE", "#90AD1C", "#3283FE", "#0c8999"))
dev.off()

############################### HERVS BULK ONLY ################################


prostate.cpm.herv <- counts.cpm.herv[, c("GTEX-12BJ1-1226-SM-5LUAE", "GTEX-15CHR-1226-SM-79OON",
                                    "GTEX-1HSMQ-1826-SM-A9SMJ", "GTEX-1I1GU-1126-SM-A96S9")]

prostate.cpm.herv.sc <- pseudobulk.herv.cpm.raw[, c("GTEX_12BJ1_5007_SM_H8L6U", "GTEX_15CHR_5014_SM_H5JDU",
                                              "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1I1GU_5006_SM_G8XQC")]

row.names(Prostate_bulk_only) <- Prostate_bulk_only$TE
row.names(Prostate_both) <- Prostate_both$TE
row.names(Prostate_sc_only) <- Prostate_sc_only$TE
row.names(Prostate_sc_only)<- gsub("_", "-", row.names(Prostate_sc_only))

prostate.cpm.herv.bulk.only <-
  as.data.frame(prostate.cpm.herv[rownames(prostate.cpm.herv) %in% 
                                    row.names(Prostate_bulk_only), ])

prostate.cpm.herv.both <-
  as.data.frame(prostate.cpm.herv[rownames(prostate.cpm.herv) %in% 
                                    row.names(Prostate_both), ])

prostate.cpm.herv.sc.only <-
  as.data.frame(prostate.cpm.herv.sc[rownames(prostate.cpm.herv.sc) %in% 
                                    row.names(Prostate_sc_only), ])


pdf("plots/09-prostate_bulk_only_herv_cpm.pdf", height=11, width=10)
ggplot(melt(prostate.cpm.herv.bulk.only), 
       aes(factor(variable), value)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free") +
  theme_pubclean() +
  xlab("Prostate Bulk Sample") +
  ylab("CPM of HERVs Unique to Bulk") +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 5)) +
  theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.y = element_text(size=10))
dev.off()

pdf("plots/09-prostate_shared_herv_cpm.pdf", height=11, width=10)
ggplot(melt(prostate.cpm.herv.both), 
       aes(factor(variable), value)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free") +
  theme_pubclean() +
  xlab("Prostate Bulk Sample") +
  ylab("CPM of HERVs Shared in Bulk and Single Cell")  +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 5)) +
  theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.y = element_text(size=10))
dev.off()

pdf("plots/09-prostate_sc_only_herv_cpm.pdf", height=11, width=10)
ggplot(melt(prostate.cpm.herv.sc.only), 
       aes(factor(variable), value)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free") +
  theme_pubclean() +
  xlab("Prostate Single Cell Sample") +
  ylab("CPM of HERVs Unique to Single Cell") +fpcols
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 5)) +
  theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.y = element_text(size=10))
dev.off()

############################ TOP BULK HERVs in SC ##############################


pdf("plots/09-top_bulk_in_sc_feature_plots_GTEX−12BJ1.pdf", height=15, width=18)
FeaturePlot(prostate.norm.merged, c("LTR46-Xq11.1", "HML6-19q13.43b", "HML4-16p13.3", 
                                 "HERVL-22q13.31", "HERVK11-15q15.1", "MER4-22q12.3",
                                 "ERVLB4-2q13c", "MER4-3q29f", "ERV316A3-6q24.1a",
                                 "HML5-12q23.1"),
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 

dev.off()

pdf("plots/09-top_bulk_in_sc_violin_plots_GTEX−12BJ1.pdf", height=15, width=18)
VlnPlot(prostate.norm.merged, features = c("LTR46-Xq11.1", "HML6-19q13.43b", "HML4-16p13.3", 
                                            "HERVL-22q13.31", "HERVK11-15q15.1", "MER4-22q12.3",
                                            "ERVLB4-2q13c", "MER4-3q29f", "ERV316A3-6q24.1a",
                                            "HML5-12q23.1"))
dev.off()

pdf("plots/09-top_bulk_in_sc_feature_plots_GTEX−15CHR.pdf", height=15, width=18)
FeaturePlot(prostate.norm.merged, c("LTR46-Xq11.1", "HML6-19q13.43b", "HML2-22q11.23", 
                                    "HERVL-22q13.31", "MER4-3q29f", "HARLEQUIN-17q21.31",
                                    "HERVK11-15q15.1", "HML4-8q24.3", "HML2-17p13.1",
                                    "ERV316A3-6q24.1a"),
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

pdf("plots/09-top_bulk_in_sc_violin_plots_GTEX−15CHR.pdf", height=15, width=18)
VlnPlot(prostate.norm.merged, features = c("LTR46-Xq11.1", "HML6-19q13.43b", "HML2-22q11.23", 
                                           "HERVL-22q13.31", "MER4-3q29f", "HARLEQUIN-17q21.31",
                                           "HERVK11-15q15.1", "HML4-8q24.3", "HML2-17p13.1",
                                           "ERV316A3-6q24.1a"))
dev.off()


pdf("plots/09-top_bulk_in_sc_feature_plots_GTEX−1HSMQ.pdf", height=15, width=18)
FeaturePlot(prostate.norm.merged, c("MER4-3q29f", "HERVL-22q13.31", "HML6-19q13.43b", 
                                    "HERVK11-15q15.1", "LTR46-Xq11.1", "ERVLE-15q25.3b",
                                    "HARLEQUIN-17q21.31", "HML4-16p13.3", "HML2-3q12.3",
                                    "HML4-8q24.3"),
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()

pdf("plots/09-top_bulk_in_sc_violin_plots_GTEX−1HSMQ.pdf", height=15, width=18)
VlnPlot(prostate.norm.merged, features = c("MER4-3q29f", "HERVL-22q13.31", "HML6-19q13.43b", 
                                           "HERVK11-15q15.1", "LTR46-Xq11.1", "ERVLE-15q25.3b",
                                           "HARLEQUIN-17q21.31", "HML4-16p13.3", "HML2-3q12.3",
                                           "HML4-8q24.3"))
dev.off()

pdf("plots/09-top_bulk_in_sc_feature_plots_GTEX−1I1GU.pdf", height=15, width=18)
FeaturePlot(prostate.norm.merged, c("HERVE-12q13.13", "HERVH-Xq25b", "HML6-19q13.43b", 
                                    "LTR46-Xq11.1", "HERVE-20p11.21b", "HERVE-17q11.2",
                                    "HML2-17p13.1", "HERVS71-19p12a", "ERV316A3-6q24.1a",
                                    "HML4-16p13.3"),
            cols=fpcols, ncol=4, raster=TRUE, pt.size = 2) 
dev.off()


pdf("plots/09-top_bulk_in_sc_violin_plots_GTEX−1I1GU.pdf", height=15, width=18)
VlnPlot(prostate.norm.merged, features = c("HERVE-12q13.13", "HERVH-Xq25b", "HML6-19q13.43b", 
                                           "LTR46-Xq11.1", "HERVE-20p11.21b", "HERVE-17q11.2",
                                           "HML2-17p13.1", "HERVS71-19p12a", "ERV316A3-6q24.1a",
                                           "HML4-16p13.3"))
dev.off()

FeaturePlot(prostate.norm.merged, c("TUBA3D", "HP", "ZDHHC8P1", "SLC29A4", "GAL", "CD177",
                                    "GAS1", "SCNN1B", "MT1M", "SLC7A2", "PDK4"))

############################## SAVE R DATA FILES ###############################

save(prostate.norm.merged, prostate.norm.merged.gtex, prostate.norm.merged.markers,
     all_gtex_metadata, gtex_metadata, file = "r_outputs/09-prostate_merged.RData")

load("r_outputs/09-prostate_merged.RData")

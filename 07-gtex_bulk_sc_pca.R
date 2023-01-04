################################################################################
################################################################################
################################################################################
################################################################################
##############################  GTEx: Bulk / SC PCA ############################

#################################### SETUP #####################################

library(tidyverse)
library(tibble)
library(gridExtra)
library(cowplot)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(PCAtools)
library(DESeq2)

################################## LOAD DATA ###################################

load("r_outputs/01-mean_cpm_by_tissue_type.RData")
load("r_outputs/01-counts.Rdata")
load("r_outputs/03-mean_raw_scTE_cpm_by_tissue_type.RData")
load("r_outputs/03-scgtex_seurat_counts.Rdata")

################################ METADATA SETUP ################################

row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

bulk_samples <- samples
sc_samples <-samples

sc_samples <- tibble::rownames_to_column(sc_samples, "remove")
rownames(sc_samples) <- sc_samples$sn_RNAseq_names
sc_samples <- subset(sc_samples, select = -c(remove))
sc_samples$type <- "sc"

bulk_samples <- tibble::rownames_to_column(bulk_samples, "remove")
rownames(bulk_samples) <- sc_samples$bulk_RNAseq
bulk_samples <- subset(bulk_samples, select = -c(remove))
bulk_samples$type <- "bulk"

metadata <- rbind(bulk_samples, sc_samples)

################################## DATA SETUP ##################################

row.names(pseudobulk.rtx.counts.raw) <- gsub("-", "_", row.names(pseudobulk.rtx.counts.raw))

combined_counts <- merge(counts.rtx, pseudobulk.rtx.counts.raw, 
                    by.x="row.names", by.y="row.names",
                    all.y=TRUE)

rownames(combined_counts) <- combined_counts[,1]
combined_counts[,1] <- NULL

stopifnot(all(colnames(combined_counts) == rownames(metadata)))

################################# COLOR SETUP ##################################

mypal = pal_npg("nrc", alpha = 0.7)(8)
mycols <- c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2",
            "#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2")

tissues <- c("Prostate", "Heart", "Lung", "E_Muscularis",
             "E_Mucosa" , "Sk_muscle", "Skin", "Breast")

tissue_color<- data.frame(tissues, mycols)

samples$color <- 
  tissue_color$mycols[match(samples$tissue, tissue_color$tissues)]

#################################### DESEQ #####################################

dds <- DESeq2::DESeqDataSetFromMatrix(countData = combined_counts,
                                      colData = metadata,
                                      design = ~ 1)

dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

normalized_counts <- counts(dds, normalized=TRUE)

# dds accounting for sequencing type
dds_acc <- DESeq2::DESeqDataSetFromMatrix(countData = combined_counts,
                                      colData = metadata,
                                      design = ~ type + tissue)

dds_acc <- DESeq2::DESeq(dds_acc, parallel=T)
tform_acc <- DESeq2::varianceStabilizingTransformation(dds_acc, blind=FALSE)

normalized_counts_acc <- counts(dds_acc, normalized=TRUE)

save(dds, tform, normalized_counts,
     dds_acc, tform_acc, normalized_counts_acc, 
     metadata,
     file="r_outputs/07-deseq_bulk_pseudobulk.Rdata")

##################################### PCA ######################################

removeVar <- 0.1
pca.obj <- 
  PCAtools::pca(assay(tform), metadata=metadata, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
            removeVar*100, length(pca.obj$xvars)))

varline <- 50
varline.x <- min(which(cumsum(pca.obj$variance) >= varline))

horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
elbow <- PCAtools::findElbowPoint(pca.obj$variance)

PCAtools::screeplot(pca.obj,
                    axisLabSize = 6,
                    components = getComponents(pca.obj, 1:30),
                    title="Retrotranscriptome SCREE NCI + TCGA",
                    hline=varline, vline=c(varline.x, horn$n, elbow)
) +
  geom_label(aes(x=varline.x+1, y=50, 
                 label = paste0(varline, '% var'), vjust = -1)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1))


cat(sprintf('%d PCs for Elbow method\n', elbow)) # 5
cat(sprintf('%d PCs for Horn method\n', horn$n)) # 9 
cat(sprintf('%d PCs needed to explain %d percent of variation\n', varline.x, varline)) # 5


################################### PCA ACC ####################################

removeVar <- 0.1
pca.obj.acc <- 
  PCAtools::pca(assay(tform_acc), metadata=metadata, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
            removeVar*100, length(pca.obj.acc$xvars)))

varline <- 50
varline.x <- min(which(cumsum(pca.obj.acc$variance) >= varline))

horn <- PCAtools::parallelPCA(assay(tform_acc), removeVar = removeVar)
elbow <- PCAtools::findElbowPoint(pca.obj.acc$variance)

PCAtools::screeplot(pca.obj.acc,
                    axisLabSize = 6,
                    components = getComponents(pca.obj.acc, 1:30),
                    title="Retrotranscriptome SCREE NCI + TCGA",
                    hline=varline, vline=c(varline.x, horn$n, elbow)
) +
  geom_label(aes(x=varline.x+1, y=50, 
                 label = paste0(varline, '% var'), vjust = -1)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1))


cat(sprintf('%d PCs for Elbow method\n', elbow)) # 4
cat(sprintf('%d PCs for Horn method\n', horn$n)) # 9 
cat(sprintf('%d PCs needed to explain %d percent of variation\n', varline.x, varline)) # 3

################################## PAIRS PLOT ###################################

pairsplot(pca.obj, getComponents(pca.obj, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "tissue",
          colkey = c("Prostate" = "#E64B35B2", "Heart" = "#4DBBD5B2",
                     "Lung" = "#00A087B2", "E_Muscularis" = "3C5488B2",
                     "E_Mucosa" = "#F39B7FB2", "Sk_muscle" = "#8491B4B2",
                     "Skin" = "#91D1C2B2", "Breast" = "#DC0000B2"),
          trianglelabSize = 12,
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_tissue.pdf", height=8, width=9)

pairsplot(pca.obj, getComponents(pca.obj, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "type",
          colkey = pal_tron("legacy", alpha = 0.7)(2),
          trianglelabSize = 12,
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_seqtype.pdf", height=8, width=9)

pairsplot(pca.obj, getComponents(pca.obj, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "participant_id",
          colkey = pal_igv(alpha=1)(16),
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_participant.pdf", height=8, width=9)

################################ PAIRS PLOT ACC ################################

pairsplot(pca.obj.acc, getComponents(pca.obj.acc, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "tissue",
          colkey = c("Prostate" = "#E64B35B2", "Heart" = "#4DBBD5B2",
                     "Lung" = "#00A087B2", "E_Muscularis" = "3C5488B2",
                     "E_Mucosa" = "#F39B7FB2", "Sk_muscle" = "#8491B4B2",
                     "Skin" = "#91D1C2B2", "Breast" = "#DC0000B2"),
          trianglelabSize = 12,
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_tissue_acc.pdf", height=8, width=9)

pairsplot(pca.obj.acc, getComponents(pca.obj.acc, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "type",
          colkey = pal_tron("legacy", alpha = 0.7)(2),
          trianglelabSize = 12,
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_seqtype_acc.pdf", height=8, width=9)

pairsplot(pca.obj.acc, getComponents(pca.obj.acc, seq_len(5)),
          labSize = 2, 
          pointSize = 0.9,
          colby = "participant_id",
          colkey = pal_igv(alpha=1)(16),
          gridlines.major = FALSE,
          gridlines.minor = FALSE)

ggsave("plots/pairsplot_participant_acc.pdf", height=8, width=9)


################################### BI PLOTS ###################################

biplot(pca.obj, 
       x="PC3",
       y="PC4",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "tissue",
       colkey = c("Prostate" = "#E64B35B2", "Heart" = "#4DBBD5B2",
                  "Lung" = "#00A087B2", "E_Muscularis" = "3C5488B2",
                  "E_Mucosa" = "#F39B7FB2", "Sk_muscle" = "#8491B4B2",
                  "Skin" = "#91D1C2B2", "Breast" = "#DC0000B2"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_tissue_pc3_pc4.pdf", height=5, width=7)


biplot(pca.obj, 
       x="PC3",
       y="PC4",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "type",
       colkey = pal_tron("legacy", alpha = 0.7)(2),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_seqtype_pc3_pc4.pdf", height=5, width=7)

biplot(pca.obj, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "type",
       ellipse = TRUE,
       ellipseConf = 0.95,
       xlim = c(-150, 200),
       ylim= c(-150,220),
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       colkey = pal_tron("legacy", alpha = 0.7)(2),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_seqtype_pc1_pc2.pdf", height=5, width=7)

biplot(pca.obj, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       xlim = c(-150, 200),
       ylim= c(-150,220),
       colby = "tissue",
       colkey = c("Prostate" = "#E64B35B2", "Heart" = "#4DBBD5B2",
                  "Lung" = "#00A087B2", "E_Muscularis" = "3C5488B2",
                  "E_Mucosa" = "#F39B7FB2", "Sk_muscle" = "#8491B4B2",
                  "Skin" = "#91D1C2B2", "Breast" = "#DC0000B2"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_tissue_pc1_pc2.pdf", height=5, width=7)

################################# BI PLOTS ACC #################################

biplot(pca.obj.acc, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "type",
       ellipse = TRUE,
       ellipseConf = 0.95,
       xlim = c(-100, 150),
       ylim= c(-75,120),
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       colkey = pal_tron("legacy", alpha = 0.7)(2),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_seqtype_pc1_pc2_acc.pdf", height=5, width=7)

biplot(pca.obj.acc, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       
       colby = "tissue",
       xlim = c(-100, 150),
       ylim= c(-75,120),
       colkey = c("Prostate" = "#E64B35B2", "Heart" = "#4DBBD5B2",
                  "Lung" = "#00A087B2", "E_Muscularis" = "3C5488B2",
                  "E_Mucosa" = "#F39B7FB2", "Sk_muscle" = "#8491B4B2",
                  "Skin" = "#91D1C2B2", "Breast" = "#DC0000B2"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right")

ggsave("plots/biplot_tissue_pc1_pc2_acc.pdf", height=5, width=7)


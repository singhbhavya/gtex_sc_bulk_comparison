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

bulk_samples <- tibble::rownames_to_column(bulk_samples, "remove")
rownames(bulk_samples) <- sc_samples$bulk_RNAseq
bulk_samples <- subset(bulk_samples, select = -c(remove))

metadata <- rbind(bulk_samples, sc_samples)

################################## DATA SETUP ##################################

combined_counts <- cbind(counts.rtx, pseudobulk.rtx.counts.raw)

stopifnot(all(colnames(combined_cpm) == rownames(metadata)))

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
                                      design = ~ tissues)

dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

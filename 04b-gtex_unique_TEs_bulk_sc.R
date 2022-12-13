################################################################################
################################################################################
################################################################################
################################################################################
#########################  GTEx: Common TEs in Bulk and SC #####################

#################################### SETUP #####################################

library(tidyverse)
library(tibble)
library(gridExtra)
library(cowplot)
library(wesanderson)
library(ggsci)
library(ggplot2)
library(ggtext)
library(ggvenn)

################################## LOAD DATA ###################################

load("r_outputs/01-counts.Rdata")
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

# Merge counts
row.names(pseudobulk.rtx.counts.raw) <- gsub("-", "_", row.names(pseudobulk.rtx.counts.raw))

combined_counts <- merge(counts.rtx, pseudobulk.rtx.counts.raw, 
                         by.x="row.names", by.y="row.names",
                         all.y=TRUE)

rownames(combined_counts) <- combined_counts[,1]
combined_counts[,1] <- NULL

stopifnot(all(colnames(combined_counts) == rownames(metadata)))

# Merge cpm
row.names(pseudobulk.rtx.cpm.raw) <- gsub("-", "_", row.names(pseudobulk.rtx.cpm.raw))

combined_cpm <- merge(counts.cpm, pseudobulk.rtx.cpm.raw, 
                         by.x="row.names", by.y="row.names",
                         all.y=TRUE)

rownames(combined_cpm) <- combined_cpm[,1]
combined_cpm[,1] <- NULL

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

################################ UNIQUE COUNTS #################################

long_counts <- combined_counts
long_counts$te <- rownames(combined_counts)

long_counts <-
  long_counts %>% pivot_longer(
  cols = 1:50,
  names_to = "sample",
  values_to  = "value"
  )

long_counts$tissue <- metadata$tissue[match(long_counts$sample, row.names(metadata))]
long_counts$type <- metadata$type[match(long_counts$sample, row.names(metadata))]

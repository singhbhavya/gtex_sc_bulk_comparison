################################################################################
################################################################################
################################################################################
################################################################################
##########################  GTEx: Bulk / SC Correlation ########################

#################################### SETUP #####################################

library(tidyverse)
library(tibble)
library(gridExtra)
library(cowplot)
library(wesanderson)
library(ggsci)
library(ggplot2)
library(ggpubr)

################################## LOAD DATA ###################################

load("r_outputs/01-mean_cpm_by_tissue_type.RData")
load("r_outputs/01-counts.Rdata")
load("r_outputs/03-mean_raw_scTE_cpm_by_tissue_type.RData")
load("r_outputs/03-scgtex_seurat_counts.Rdata")

################################## DATA SETUP ##################################

row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

# Keep only genes 

################################# COLOR SETUP ##################################

mypal = pal_npg("nrc", alpha = 0.7)(8)
mycols <- c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2",
            "#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2")

tissues <- c("Prostate", "Heart", "Lung", "E_Muscularis",
             "E_Mucosa" , "Sk_muscle", "Skin", "Breast")

tissue_color<- data.frame(tissues, mycols)

samples$color <- 
  tissue_color$mycols[match(samples$tissue, tissue_color$tissues)]

############################# SC / BULK CORRELATION ############################
################################ SAMPLE LEVEL TE ###############################

all_bulk_sc <- data.frame()

bulk_sc_cor <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_tes <- as.data.frame(counts.cpm.rtx[,bulk_samp])
  colnames(bulk_tes) <- c(bulk_samp)
  
  sc_tes <- as.data.frame(pseudobulk.rtx.cpm.raw[,sc_samp, drop=FALSE])
  colnames(sc_tes) <- c(sc_samp)
  
  samp_tes <- cbind(bulk_tes, sc_tes)
  
  samp_te_df <- samp_tes
  samp_te_df$bulk <- bulk_samp
  samp_te_df$sc <- sc_samp
  samp_te_df$tissue <- samples$tissue[i]
  samp_te_df$participant <- samples$participant_id[i]
  colnames(samp_te_df) <- c("bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                            "participant")
  row.names(samp_te_df) <- NULL
  
  all_bulk_sc <<- rbind(all_bulk_sc, samp_te_df)
  
  cor.test(samp_tes[,1], samp_tes[,2], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_tes, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            color=samples$color[samples$bulk_RNAseq == bulk_samp],
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
  
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_corr.pdf", height=20, width=20)

########################## SC / BULK CORRELATION ALL TE #########################

ggscatter(all_bulk_sc, x = "bulk_counts", y = "sc_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          exact=FALSE,
          xlab = "Bulk CPM", 
          ylab = "Single Cell CPM") +
  font("xlab", size=10) +
  font("ylab", size=10)

ggsave("plots/all_bulk_sc_corr.pdf", height=4, width=4)

############################# SC / BULK CORRELATION ############################
############################### SAMPLE LEVEL GENE ##############################

all_bulk_sc_gene <- data.frame()

bulk_sc_cor_gene <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_tes <- as.data.frame(counts.cpm.rtx[,bulk_samp])
  colnames(bulk_tes) <- c(bulk_samp)
  
  sc_tes <- as.data.frame(pseudobulk.rtx.cpm.raw[,sc_samp, drop=FALSE])
  colnames(sc_tes) <- c(sc_samp)
  
  samp_tes <- cbind(bulk_tes, sc_tes)
  
  samp_te_df <- samp_tes
  samp_te_df$bulk <- bulk_samp
  samp_te_df$sc <- sc_samp
  samp_te_df$tissue <- samples$tissue[i]
  samp_te_df$participant <- samples$participant_id[i]
  colnames(samp_te_df) <- c("bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                            "participant")
  row.names(samp_te_df) <- NULL
  
  all_bulk_sc <<- rbind(all_bulk_sc, samp_te_df)
  
  cor.test(samp_tes[,1], samp_tes[,2], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_tes, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            color=samples$color[samples$bulk_RNAseq == bulk_samp],
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
  
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_corr.pdf", height=20, width=20)
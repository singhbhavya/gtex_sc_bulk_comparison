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

################################## DATA SETUP ##################################

row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

############################# SC / BULK CORRELATION ############################

bulk_sc_cor <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_tes <- as.data.frame(counts.cpm.rtx[,bulk_samp])
  colnames(bulk_tes) <- c(bulk_samp)
  
  sc_tes <- as.data.frame(pseudobulk.rtx.cpm.raw[,sc_samp, drop=FALSE])
  colnames(sc_tes) <- c(sc_samp)
  
  samp_tes <- cbind(bulk_tes, sc_tes)
  
  cor.test(samp_tes[,1], samp_tes[,2], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_tes, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_corr.pdf", height=20, width=20)
################################################################################
################################################################################
################################################################################
################################################################################
##########################  GTEx: Top TEs in Bulk and SC #######################

#################################### SETUP #####################################

library(tidyverse)
library(tibble)
library(gridExtra)
library(cowplot)
library(wesanderson)
library(ggsci)

################################## LOAD DATA ###################################

load("r_outputs/01-mean_cpm_by_tissue_type.RData")
load("r_outputs/01-counts.Rdata")
load("r_outputs/03-mean_raw_scTE_cpm_by_tissue_type.RData")

################################## DATA SETUP ##################################

row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

################################# COLOR SETUP ##################################

mypal = pal_npg("nrc", alpha = 0.7)(8)
mycols <- c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2",
                   "#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2")

tissues <- c("Prostate", "Heart", "Lung", "E_Muscularis",
             "E_Mucosa" , "Sk_muscle", "Skin", "Breast")

tissue_color<- data.frame(tissues, mycols)

samples$color <- 
  tissue_color$mycols[match(samples$tissue, tissue_color$tissues)]

######################## TOP TEs PER TISSUE TYPE (BULK) ########################
  
top_tes_per_tissue <- function(i) {
  
  tissue_type_df <- tibble::rownames_to_column(bulk_list[[i]], "TEs")
  
  ggplot(tissue_type_df[1:10,], 
         aes(x= reorder(TEs, -mean_counts), 
             y=mean_counts)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("Average CPM") +
    ggtitle(names(bulk_list)[i]) + 
    geom_bar(stat='identity', 
             fill = tissue_color$mycols[tissues == names(bulk_list)[i]]) 
  
  }
  
bulk_list <- list(Breast = Breast.cpm, E_Mucosa = E_Mucosa.cpm, 
                    E_Muscularis = E_Muscularis.cpm, Heart = Heart.cpm,
                    Lung = Lung.cpm, Prostate = Prostate.cpm, 
                    Sk_muscle = Sk_muscle.cpm, Skin = Skin.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(bulk_list), top_tes_per_tissue), 
                   ncol=4)
  
ggsave("plots/top_tes_per_tissue_bulk.pdf", height=10, width=17)
remove(bulk_list)


######################### TOP TEs PER TISSUE TYPE (SC) #########################

top_tes_per_tissue_sc <- function(i) {
  
  tissue_type_df <- tibble::rownames_to_column(sc_list[[i]], "TEs")
  
  ggplot(tissue_type_df[1:10,], 
         aes(x= reorder(TEs, -mean_counts), 
             y=mean_counts)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("Average CPM") +
    ggtitle(names(sc_list)[i]) + 
    geom_bar(stat='identity', 
             fill = tissue_color$mycols[tissues == names(sc_list)[i]]) 
  
}

sc_list <- list(Breast = Breast.sc.cpm, E_Mucosa = E_Mucosa.sc.cpm, 
                  E_Muscularis = E_Muscularis.sc.cpm, Heart = Heart.sc.cpm,
                  Lung = Lung.sc.cpm, Prostate = Prostate.sc.cpm, 
                  Sk_muscle = Sk_muscle.sc.cpm, Skin = Skin.sc.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(sc_list), top_tes_per_tissue_sc), 
                   ncol=4)

ggsave("plots/top_tes_per_tissue_sc.pdf", height=10, width=17)
remove(sc_list)

########################### TOP TEs PER SAMPLE (BULK) ##########################

top_tes_per_sample <- function(i) {
  
  sample_counts <- as.data.frame(counts.cpm.rtx[,i])
  colnames(sample_counts) <- "counts"
  sample_counts <- tibble::rownames_to_column(sample_counts, "TEs")
  sample_counts <- sample_counts[order(sample_counts$counts, decreasing = TRUE),] 
  
  ggplot(sample_counts[1:10,], 
         aes(x= reorder(TEs, -counts), 
             y=counts)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(i) + 
    geom_bar(stat='identity',
             fill=samples$color[samples$bulk_RNAseq == i]) 
}

cowplot::plot_grid(plotlist = lapply(samples$bulk_RNAseq, top_tes_per_sample), 
                   ncol=5)

ggsave("plots/top_tes_per_sample_bulk.pdf", height=20, width=20)

############################ TOP TEs PER SAMPLE (SC) ###########################

top_tes_per_sample_sc <- function(i) {
  
  sample_counts <- as.data.frame(pseudobulk.rtx.cpm.raw[,i, 
                                                        drop=FALSE])
  colnames(sample_counts) <- "counts"
  sample_counts <- tibble::rownames_to_column(sample_counts, "TEs")
  sample_counts <- sample_counts[order(sample_counts$counts, decreasing = TRUE),] 
  
  ggplot(sample_counts[1:10,], 
         aes(x= reorder(TEs, -counts), 
             y=counts)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(i) + 
    geom_bar(stat='identity',
             fill=samples$color[samples$sn_RNAseq_names == i]) 
}

cowplot::plot_grid(plotlist = lapply(samples$sn_RNAseq_names, top_tes_per_sample_sc), 
                   ncol=5)

ggsave("plots/top_tes_per_sample_sc.pdf", height=20, width=20)


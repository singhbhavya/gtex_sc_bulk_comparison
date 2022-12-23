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
library(ggplot2)
library(ggtext)

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

####################### TOP HERVs PER TISSUE TYPE (BULK) #######################

bulk_list <- list(Breast = Breast.herv.cpm, E_Mucosa = E_Mucosa.herv.cpm, 
                  E_Muscularis = E_Muscularis.herv.cpm, Heart = Heart.herv.cpm,
                  Lung = Lung.herv.cpm, Prostate = Prostate.herv.cpm, 
                  Sk_muscle = Sk_muscle.herv.cpm, Skin = Skin.herv.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(bulk_list), top_tes_per_tissue), 
                   ncol=4)

ggsave("plots/top_HERVs_per_tissue_bulk.pdf", height=10, width=17)
remove(bulk_list)

####################### TOP LINEs PER TISSUE TYPE (BULK) #######################

bulk_list <- list(Breast = Breast.l1.cpm, E_Mucosa = E_Mucosa.l1.cpm, 
                  E_Muscularis = E_Muscularis.l1.cpm, Heart = Heart.l1.cpm,
                  Lung = Lung.l1.cpm, Prostate = Prostate.l1.cpm, 
                  Sk_muscle = Sk_muscle.l1.cpm, Skin = Skin.l1.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(bulk_list), top_tes_per_tissue), 
                   ncol=4)

ggsave("plots/top_LINEs_per_tissue_bulk.pdf", height=10, width=17)
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

######################## TOP HERVs PER TISSUE TYPE (SC) ########################

sc_list <- list(Breast = Breast.sc.herv.cpm, E_Mucosa = E_Mucosa.sc.herv.cpm, 
                E_Muscularis = E_Muscularis.sc.herv.cpm, Heart = Heart.sc.herv.cpm,
                Lung = Lung.sc.herv.cpm, Prostate = Prostate.sc.herv.cpm, 
                Sk_muscle = Sk_muscle.sc.herv.cpm, Skin = Skin.sc.herv.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(sc_list), top_tes_per_tissue_sc), 
                   ncol=4)

ggsave("plots/top_HERVs_per_tissue_sc.pdf", height=10, width=17)
remove(sc_list)

######################## TOP LINEs PER TISSUE TYPE (SC) ########################

sc_list <- list(Breast = Breast.sc.l1.cpm, E_Mucosa = E_Mucosa.sc.l1.cpm, 
                E_Muscularis = E_Muscularis.sc.l1.cpm, Heart = Heart.sc.l1.cpm,
                Lung = Lung.sc.l1.cpm, Prostate = Prostate.sc.l1.cpm, 
                Sk_muscle = Sk_muscle.sc.l1.cpm, Skin = Skin.sc.l1.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(sc_list), top_tes_per_tissue_sc), 
                   ncol=4)

ggsave("plots/top_LINEs_per_tissue_sc.pdf", height=10, width=17)
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

########################## TOP HERVs PER SAMPLE (BULK) #########################

top_hervs_per_sample <- function(i) {
  
  sample_counts <- as.data.frame(counts.cpm.herv[,i])
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

cowplot::plot_grid(plotlist = lapply(samples$bulk_RNAseq, top_hervs_per_sample), 
                   ncol=5)

ggsave("plots/top_HERVs_per_sample_bulk.pdf", height=20, width=20)

########################## TOP LINEs PER SAMPLE (BULK) #########################

top_lines_per_sample <- function(i) {
  
  sample_counts <- as.data.frame(counts.cpm.l1[,i])
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

cowplot::plot_grid(plotlist = lapply(samples$bulk_RNAseq, top_lines_per_sample), 
                   ncol=5)

ggsave("plots/top_LINEs_per_sample_bulk.pdf", height=20, width=20)

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

########################### TOP HERVs PER SAMPLE (SC) ##########################

top_hervs_per_sample_sc <- function(i) {
  
  sample_counts <- as.data.frame(pseudobulk.herv.cpm.raw[,i, 
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

cowplot::plot_grid(plotlist = lapply(samples$sn_RNAseq_names, top_hervs_per_sample_sc), 
                   ncol=5)

ggsave("plots/top_HERVs_per_sample_sc.pdf", height=20, width=20)

########################### TOP LINEs PER SAMPLE (SC) ##########################

top_lines_per_sample_sc <- function(i) {
  
  sample_counts <- as.data.frame(pseudobulk.l1.cpm.raw[,i, 
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

cowplot::plot_grid(plotlist = lapply(samples$sn_RNAseq_names, top_lines_per_sample_sc), 
                   ncol=5)

ggsave("plots/top_LINEs_per_sample_sc.pdf", height=20, width=20)

######################### TOP TEs PER SAMPLE (BULK & SC) #######################

top_bulk_tes_per_sample <- function(i) {
  
  sample_counts_top_bulk <- as.data.frame(counts.cpm.rtx[,samples$bulk_RNAseq[i]])
  sample_counts_sc <- as.data.frame(pseudobulk.rtx.cpm.raw[,samples$sn_RNAseq_names[i], 
                                                         drop = FALSE])
  
  row.names(sample_counts_sc) <- gsub("-", "_", row.names(sample_counts_sc))
  colnames(sample_counts_top_bulk) <- "top_bulk"
  colnames(sample_counts_sc) <- "sc_for_top_bulk"
  sample_counts_top_bulk <- tibble::rownames_to_column(sample_counts_top_bulk, "TEs")
  sample_counts_sc <- tibble::rownames_to_column(sample_counts_sc, "TEs")
  sample_counts <- merge(sample_counts_top_bulk,sample_counts_sc, by = "TEs")
  sample_counts <- sample_counts[order(sample_counts$top_bulk, 
                                     decreasing = TRUE),] 
  sample_counts <- sample_counts[1:10,]
  sample_counts <- 
    sample_counts %>% pivot_longer(
      cols = c("top_bulk", "sc_for_top_bulk"),
      names_to = c("count_type"),
      values_to = "counts") 
  
  sample_counts$TEs <- factor(sample_counts$TEs, levels = unique(sample_counts$TEs))
  sample_counts$count_type <- factor(sample_counts$count_type, levels =
                                     c("top_bulk", "sc_for_top_bulk"))
  
  ggplot(sample_counts, 
         aes(x= TEs, y=counts, fill=count_type))+
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(paste0(samples$participant_id[i], " ", samples$tissue[i])) +
    theme(plot.title = element_textbox_simple()) +
    geom_bar(stat='identity', position="dodge") +
    theme(legend.position="none")

}

bulk_plot <-
  cowplot::plot_grid(plotlist = lapply(seq_along(1:25), top_bulk_tes_per_sample))

legend_bulk <- get_legend(
  top_bulk_tes_per_sample(1) + 
    theme(legend.position="right") + 
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(bulk_plot, legend_bulk, rel_widths = c(3, .4))

ggsave("plots/top_bulk_TEs_per_sample_plus_sc.pdf", height=20, width=20)

######################### TOP TEs PER SAMPLE (SC & BULK) #######################

top_sc_tes_per_sample <- function(i) {
  
  sample_counts_top_sc <- as.data.frame(pseudobulk.rtx.cpm.raw[,samples$sn_RNAseq_names[i], 
                                                               drop = FALSE])
  sample_counts_bulk <- as.data.frame(counts.cpm.rtx[,samples$bulk_RNAseq[i]])
  
  row.names(sample_counts_top_sc) <- gsub("-", "_", row.names(sample_counts_top_sc))
  colnames(sample_counts_bulk) <- "bulk_for_top_sc"
  colnames(sample_counts_top_sc) <- "top_sc"
  
  sample_counts_top_sc <- tibble::rownames_to_column(sample_counts_top_sc, "TEs")
  sample_counts_bulk <- tibble::rownames_to_column(sample_counts_bulk, "TEs")
  sample_counts <- merge(sample_counts_top_sc,sample_counts_bulk, by = "TEs")
  sample_counts <- sample_counts[order(sample_counts$top_sc, 
                                       decreasing = TRUE),] 
  
  sample_counts <- sample_counts[1:10,]
  sample_counts <- 
    sample_counts %>% pivot_longer(
      cols = c("top_sc", "bulk_for_top_sc"),
      names_to = c("count_type"),
      values_to = "counts") 
  
  sample_counts$TEs <- factor(sample_counts$TEs, levels = unique(sample_counts$TEs))
  sample_counts$count_type <- factor(sample_counts$count_type, levels =
                                       c("bulk_for_top_sc", "top_sc"))
  
  ggplot(sample_counts, 
         aes(x= TEs, y=counts, fill=count_type))+
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(paste0(samples$participant_id[i], " ", samples$tissue[i])) +
    theme(plot.title = element_textbox_simple()) +
    geom_bar(stat='identity', position="dodge") +
    theme(legend.position="none")
}

sc_plot <-
  cowplot::plot_grid(plotlist = lapply(seq_along(1:25), top_sc_tes_per_sample))

legend_sc <- get_legend(
  top_sc_tes_per_sample(1) + 
    theme(legend.position="right") + 
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(sc_plot, legend_sc, rel_widths = c(3, .4))

ggsave("plots/top_sc_TEs_per_sample_plus_bulk.pdf", height=20, width=20)

######################## TOP HERVs PER SAMPLE (BULK & SC) ######################

top_bulk_hervs_per_sample <- function(i) {
  
  sample_counts_top_bulk <- as.data.frame(counts.cpm.herv[,samples$bulk_RNAseq[i]])
  sample_counts_sc <- as.data.frame(pseudobulk.herv.cpm.raw[,samples$sn_RNAseq_names[i], 
                                                           drop = FALSE])
  
  row.names(sample_counts_sc) <- gsub("-", "_", row.names(sample_counts_sc))
  colnames(sample_counts_top_bulk) <- "top_bulk"
  colnames(sample_counts_sc) <- "sc_for_top_bulk"
  sample_counts_top_bulk <- tibble::rownames_to_column(sample_counts_top_bulk, "TEs")
  sample_counts_sc <- tibble::rownames_to_column(sample_counts_sc, "TEs")
  sample_counts <- merge(sample_counts_top_bulk,sample_counts_sc, by = "TEs")
  sample_counts <- sample_counts[order(sample_counts$top_bulk, 
                                       decreasing = TRUE),] 
  sample_counts <- sample_counts[1:10,]
  sample_counts <- 
    sample_counts %>% pivot_longer(
      cols = c("top_bulk", "sc_for_top_bulk"),
      names_to = c("count_type"),
      values_to = "counts") 
  
  sample_counts$TEs <- factor(sample_counts$TEs, levels = unique(sample_counts$TEs))
  sample_counts$count_type <- factor(sample_counts$count_type, levels =
                                       c("top_bulk", "sc_for_top_bulk"))
  
  ggplot(sample_counts, 
         aes(x= TEs, y=counts, fill=count_type))+
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(paste0(samples$participant_id[i], " ", samples$tissue[i])) +
    theme(plot.title = element_textbox_simple()) +
    geom_bar(stat='identity', position="dodge") +
    theme(legend.position="none")
  
}

bulk_herv_plot <-
  cowplot::plot_grid(plotlist = lapply(seq_along(1:25), top_bulk_hervs_per_sample))

legend_bulk_herv <- get_legend(
  top_bulk_hervs_per_sample(1) + 
    theme(legend.position="right") + 
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(bulk_herv_plot, legend_bulk_herv, rel_widths = c(3, .4))

ggsave("plots/top_bulk_hervs_per_sample_plus_sc.pdf", height=20, width=20)

######################### TOP HERVs PER SAMPLE (SC & BULK) #######################

top_sc_hervs_per_sample <- function(i) {
  
  sample_counts_top_sc <- as.data.frame(pseudobulk.herv.cpm.raw[,samples$sn_RNAseq_names[i], 
                                                               drop = FALSE])
  sample_counts_bulk <- as.data.frame(counts.cpm.herv[,samples$bulk_RNAseq[i]])
  
  row.names(sample_counts_top_sc) <- gsub("-", "_", row.names(sample_counts_top_sc))
  colnames(sample_counts_bulk) <- "bulk_for_top_sc"
  colnames(sample_counts_top_sc) <- "top_sc"
  
  sample_counts_top_sc <- tibble::rownames_to_column(sample_counts_top_sc, "TEs")
  sample_counts_bulk <- tibble::rownames_to_column(sample_counts_bulk, "TEs")
  sample_counts <- merge(sample_counts_top_sc,sample_counts_bulk, by = "TEs")
  sample_counts <- sample_counts[order(sample_counts$top_sc, 
                                       decreasing = TRUE),] 
  
  sample_counts <- sample_counts[1:10,]
  sample_counts <- 
    sample_counts %>% pivot_longer(
      cols = c("top_sc", "bulk_for_top_sc"),
      names_to = c("count_type"),
      values_to = "counts") 
  
  sample_counts$TEs <- factor(sample_counts$TEs, levels = unique(sample_counts$TEs))
  sample_counts$count_type <- factor(sample_counts$count_type, levels =
                                       c("bulk_for_top_sc", "top_sc"))
  
  ggplot(sample_counts, 
         aes(x= TEs, y=counts, fill=count_type))+
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(
            angle = 30, vjust = 1, hjust=1, size = 13)) +
    ylab("CPM") +
    ggtitle(paste0(samples$participant_id[i], " ", samples$tissue[i])) +
    theme(plot.title = element_textbox_simple()) +
    geom_bar(stat='identity', position="dodge") +
    theme(legend.position="none")
}

sc_herv_plot <-
  cowplot::plot_grid(plotlist = lapply(seq_along(1:25), top_sc_hervs_per_sample))

legend_sc_hervs <- get_legend(
  top_sc_hervs_per_sample(1) + 
    theme(legend.position="right") + 
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(sc_herv_plot, legend_sc_hervs, rel_widths = c(3, .4))

ggsave("plots/top_sc_hervs_per_sample_plus_bulk.pdf", height=20, width=20)



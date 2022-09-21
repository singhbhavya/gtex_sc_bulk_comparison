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
    geom_bar(stat='identity') 
  
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
    geom_bar(stat='identity') 
  
}

sc_list <- list(Breast = Breast.sc.cpm, E_Mucosa = E_Mucosa.sc.cpm, 
                  E_Muscularis = E_Muscularis.sc.cpm, Heart = Heart.sc.cpm,
                  Lung = Lung.sc.cpm, Prostate = Prostate.sc.cpm, 
                  Sk_muscle = Sk_muscle.sc.cpm, Skin = Skin.sc.cpm)

cowplot::plot_grid(plotlist = lapply(seq_along(sc_list), top_tes_per_tissue_sc), 
                   ncol=4)

ggsave("plots/top_tes_per_tissue_sc.pdf", height=10, width=17)
remove(sc_list)

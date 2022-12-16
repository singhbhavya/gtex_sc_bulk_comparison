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
library(ggVennDiagram)
library(UpSetR)
library(ComplexUpset)
library(ggsci)

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

################################ ANNOTATIONS ###################################
## load annotation
retro.hg38.v1 <-
  readr::read_tsv(
    "https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/genes.tsv.gz",
    na=c('.'))
retro.hg38.v1 <- retro.hg38.v1 %>%
  tidyr::separate(locus, c("family"), sep='_', remove=F, extra='drop') %>%
  dplyr::mutate(
    te_class = factor(ifelse(is.na(l1base_id), 'LTR', 'LINE'), levels=c('LTR','LINE')),
  )

# Remove the confounding LINE element (L1FLnI_Xq21.1db) that has a poly A tail
# in the middle of it:

retro.hg38.v1<-
  retro.hg38.v1[!(retro.hg38.v1$locus=="L1FLnI_Xq21.1db"),]

retro.annot <- retro.hg38.v1
row.names(retro.annot) <- retro.annot$locus
row.names(retro.annot) <- gsub("_", "-", row.names(retro.annot))


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

bulk <- sort(unique(long_counts$te[long_counts$type == "bulk"]))
sc <- sort(unique(long_counts$te[long_counts$type == "sc"]))

# Function to get unique TE lists

unique_tes_per_tissue <- function(tiss) {
  
  tiss_bulk <-sort(unique(long_counts$te[long_counts$type == "bulk" &
                               long_counts$tissue == tiss]))
  tiss_sc <- sort(unique(long_counts$te[long_counts$type == "sc" &
                                          long_counts$tissue == tiss]))

  assign(paste0(tiss, "_bulk"), tiss_bulk, envir = .GlobalEnv)
  assign(paste0(tiss, "_sc"), tiss_sc, envir = .GlobalEnv)
  
}

lapply(unique(metadata$tissue), unique_tes_per_tissue)

# Set of bulk TEs
bulk.sets <- list(Heart_bulk, Sk_muscle_bulk, Skin_bulk, E_Mucosa_bulk, 
                  E_Muscularis_bulk, Prostate_bulk, Breast_bulk, Lung_bulk)
names(bulk.sets) <- c("Heart","Skeletal muscle", "Skin", "E Mucosa",
                      "E Muscularis", "Prostate", "Breast", "Lung")

# Set of single cell TEs
sc.sets <- list(Heart_sc, Sk_muscle_sc, Skin_sc, E_Mucosa_sc, 
                E_Muscularis_sc, Prostate_sc, Breast_sc, Lung_sc)
names(sc.sets) <- c("Heart","Skeletal muscle", "Skin", "E Mucosa",
                    "E Muscularis", "Prostate", "Breast", "Lung")


################################ VENN DIAGRAMS #################################

# Venn diagram of bulk versus single cell
pdf("plots/04b-bulk_sc_venn.pdf", height=2, width=4)
ggVennDiagram(x = list(bulk,sc),  
              category.names = c("Bulk","Single cell"),
              fill_color = c("#f8766d", "#00bfc4"),
              set_color = c("#f8766d", "#00bfc4"),
              label_alpha = 0,
              label_color = "white") + 
  scale_fill_brewer(values = c("#f8766d", "#00bfc4"))
dev.off()

a <- list(`Bulk` = sort(unique(long_counts$te[long_counts$type == "bulk"])),
          `Single Cell` = sort(unique(long_counts$te[long_counts$type == "sc"])))

pdf("plots/04b-bulk_sc_venn.pdf", height=4, width=5)
ggvenn(a, c("Bulk", "Single Cell"),
       fill_color = c("#f8766d", "#00bfc4"))
dev.off()

################################ COMPLEX UPSET #################################

# Setting up dataframes
tissue_counts_bulk <- long_counts[!(long_counts$type=="sc"),]
tissue_counts_sc <- long_counts[!(long_counts$type=="bulk"),]

tissue_counts_bulk <- tissue_counts_bulk[,c("te", "tissue", "value")]
tissue_counts_bulk$value <- TRUE
tissue_counts_bulk <- unique(tissue_counts_bulk)
tissue_counts_bulk <- tissue_counts_bulk %>%
  pivot_wider(names_from = "tissue",
              values_from = "value", 
              values_fill = FALSE)
tissue_counts_bulk <- tissue_counts_bulk[,c(colnames(tissue_counts_bulk[2:9]))]

tissue_counts_sc <- tissue_counts_sc[,c("te", "tissue", "value")]
tissue_counts_sc$value <- TRUE
tissue_counts_sc <- unique(tissue_counts_sc)
tissue_counts_sc <- tissue_counts_sc %>%
  pivot_wider(names_from = "tissue",
              values_from = "value", 
              values_fill = FALSE)
tissue_counts_sc <- tissue_counts_sc[,c(colnames(tissue_counts_sc[2:9]))]

tissues = colnames(tissue_counts_bulk)

# Upset plot for bulk
pdf("plots/04b-bulk_tissue_upset.pdf", height=4, width=7)
upset(tissue_counts_bulk, tissues, name='Shared and Unique TEs Per Tissue Type in Bulk',
      intersections = list( 
        c("Heart","Sk_muscle", "Skin", "E_Mucosa",
             "E_Muscularis", "Prostate", "Breast", "Lung"), 
        c("Heart"), 
        c("Sk_muscle"), 
        c("Skin"),
        c("E_Mucosa"), 
        c("E_Muscularis"),
        c("Prostate"),
        c("Breast"),
        c("Lung")), 
      queries = list(
        upset_query(set=c("Prostate"), color = "#E64B35B2", fill = "#E64B35B2"),
        upset_query(set=c("Heart"), color = "#4DBBD5B2", fill = "#4DBBD5B2"),
        upset_query(set=c("E_Muscularis"), color = "#3C5488B2", fill = "#3C5488B2"),
        upset_query(set=c("E_Mucosa"), color = "#F39B7FB2", fill = "#F39B7FB2"),
        upset_query(set=c("Sk_muscle"), color = "#8491B4B2", fill = "#8491B4B2"),
        upset_query(set=c("Skin"), color = "#91D1C2B2", fill = "#91D1C2B2"),
        upset_query(set=c("Lung"), color = "#00A087B2", fill = "#00A087B2"),
        upset_query(set=c("Breast"), color = "#DC0000B2", fill = "#DC0000B2")
        ),
      set_sizes=(
        upset_set_size()
        + ylab('Total TEs Per Tissue') ))

dev.off()

pdf("plots/04b-sc_tissue_upset.pdf", height=4, width=7)
# Upset for single cell
upset(tissue_counts_sc, tissues, name='Shared and Unique TEs Per Tissue Type in Single Cell',
      intersections = list( 
        c("Heart","Sk_muscle", "Skin", "E_Mucosa",
          "E_Muscularis", "Prostate", "Breast", "Lung"), 
        c("Heart"), 
        c("Sk_muscle"), 
        c("Skin"),
        c("E_Mucosa"), 
        c("E_Muscularis"),
        c("Prostate"),
        c("Breast"),
        c("Lung")), 
      queries = list(
        upset_query(set=c("Prostate"), color = "#E64B35B2", fill = "#E64B35B2"),
        upset_query(set=c("Heart"), color = "#4DBBD5B2", fill = "#4DBBD5B2"),
        upset_query(set=c("E_Muscularis"), color = "#3C5488B2", fill = "#3C5488B2"),
        upset_query(set=c("E_Mucosa"), color = "#F39B7FB2", fill = "#F39B7FB2"),
        upset_query(set=c("Sk_muscle"), color = "#8491B4B2", fill = "#8491B4B2"),
        upset_query(set=c("Skin"), color = "#91D1C2B2", fill = "#91D1C2B2"),
        upset_query(set=c("Lung"), color = "#00A087B2", fill = "#00A087B2"),
        upset_query(set=c("Breast"), color = "#DC0000B2", fill = "#DC0000B2")
      ),
      set_sizes=(
        upset_set_size()
        + ylab('Total TEs Per Tissue') ))
dev.off()

###################### FAMILY-LEVEL SHARED AND UNIQUE TEs ######################

# Setting up dataframes

both <- as.data.frame(intersect(bulk, sc))
colnames(both) <- c("TE")
both$family <- retro.annot$family[match(both$TE, retro.annot$locus)]
both_only_families <-
  both %>% dplyr::count(family, sort = TRUE) 
colnames(both_only_families) <- c("family", "both")

bulk_only <- as.data.frame(setdiff(bulk, sc))
colnames(bulk_only) <- c("TE")
bulk_only$family <- retro.annot$family[match(bulk_only$TE, retro.annot$locus)]
bulk_only_families <-
  bulk_only %>% dplyr::count(family, sort = TRUE) 
colnames(bulk_only_families) <- c("family", "bulk")

sc_only <- as.data.frame(setdiff(sc, bulk))
colnames(sc_only) <- c("TE")
sc_only$family <- retro.annot$family[match(sc_only$TE, retro.annot$locus)]
sc_only_families <-
  sc_only %>% dplyr::count(family, sort = TRUE) 
colnames(sc_only_families) <- c("family", "sc")

merged_family_breakdown <- merge(bulk_only_families, 
                                 sc_only_families,
                                 by = "family", all = TRUE) %>%
  merge(both_only_families, by = "family", all = TRUE)

merged_family_breakdown[is.na(merged_family_breakdown)] <- 0

merged_family_breakdown <-
  merged_family_breakdown %>%
  pivot_longer(
    cols = c("bulk", "sc", "both"),
    names_to = "Identified_in",
    values_to = "TE_count"
  )

total_merged_family_breakdown <- 
  merged_family_breakdown %>%
  group_by(Identified_in) %>%
  summarize(total = sum(TE_count))

pdf("plots/04b-te_families_bulk_sc.pdf", height=8, width=5)
ggplot(merged_family_breakdown, aes(fill=reorder(family, -TE_count), y=Identified_in, x=TE_count)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7))) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("TEs identified") +
  xlab("Proportion of TEs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(ncol = 2))

dev.off()

############################### TE DISTRIBUTION ################################

pdf("plots/04b-violin_bulk_sc.pdf", height=3, width=4)
ggplot(long_counts, aes(type, value, fill=type)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=23, size=2) + 
  scale_y_continuous(trans = "log10") +
  theme_cowplot() +
  xlab("GTEx Sequencing Type") +
  ylab("TE Expression") + 
  theme(legend.position = c("None"))
dev.off()
  
pdf("plots/04b-violin_bulk_sc_tissue.pdf", height=3, width=5)
ggplot(long_counts, aes(tissue, value, fill=type)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=2, size=2) + 
  scale_y_continuous(trans = "log10") +
  theme_cowplot() +
  xlab("GTEx Sequencing Type") +
  ylab("TE Expression") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf("plots/04b-violin_tissue_bulk_sc.pdf", height=3, width=5)
ggplot(long_counts, aes(type, value, fill=tissue)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=23, size=2) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values=c("#DC0000B2", "#F39B7FB2", "#3C5488B2", "#4DBBD5B2",
                             "#00A087B2", "#E64B35B2", "#8491B4B2", "#91D1C2B2")) +
  xlab("GTEx Sequencing Tissue and Type") +
  ylab("TE Expression") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme_cowplot()

dev.off()

####################### TISSUE-LEVEL SHARED AND UNIQUE TEs #####################

# heart
heart_venn <-
  ggvenn(list(`Bulk` = Heart_bulk, `Single Cell` = Heart_sc),
       c("Bulk", "Single Cell"),
       fill_color = c("#f8766d", "#00bfc4"),
       set_name_size = 2, text_size = 2) + 
  ggtitle("Heart") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# lung
lung_venn <-
  ggvenn(list(`Bulk` = Lung_bulk, `Single Cell` = Lung_sc),
       c("Bulk", "Single Cell"),
       fill_color = c("#f8766d", "#00bfc4"),
       set_name_size = 2, text_size = 2) + 
  ggtitle("Lung") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# prostate 
prostate_venn <- 
  ggvenn(list(`Bulk` = Prostate_bulk, `Single Cell` = Prostate_sc),
       c("Bulk", "Single Cell"),
       fill_color = c("#f8766d", "#00bfc4"),
       set_name_size = 2, text_size = 2) + 
  ggtitle("Prostate") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# Sk_muscle
sk_muscle_venn <-
  ggvenn(list(`Bulk` = Sk_muscle_bulk, `Single Cell` = Sk_muscle_sc),
         c("Bulk", "Single Cell"),
         fill_color = c("#f8766d", "#00bfc4"),
         set_name_size = 2, text_size = 2) + 
  ggtitle("Sk muscle") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# E_muscularis
E_muscularis_venn <-
  ggvenn(list(`Bulk` = E_Muscularis_bulk, `Single Cell` = E_Muscularis_sc),
         c("Bulk", "Single Cell"),
         fill_color = c("#f8766d", "#00bfc4"),
         set_name_size = 2, text_size = 2) + 
  ggtitle("E. muscularis") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# E_mucosa
E_mucosa_venn <-
  ggvenn(list(`Bulk` = E_Mucosa_bulk, `Single Cell` = E_Mucosa_sc),
         c("Bulk", "Single Cell"),
         fill_color = c("#f8766d", "#00bfc4"),
         set_name_size = 2, text_size = 2) + 
  ggtitle("E. mucosa") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# Breast
Breast_venn <-
  ggvenn(list(`Bulk` = Breast_bulk, `Single Cell` = Breast_sc),
         c("Bulk", "Single Cell"), 
         fill_color = c("#f8766d", "#00bfc4"),
         set_name_size = 2, text_size = 2) + 
  ggtitle("Breast") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# Skin
Skin_venn <-
  ggvenn(list(`Bulk` = Skin_bulk, `Single Cell` = Skin_sc),
         c("Bulk", "Single Cell"),
         fill_color = c("#f8766d", "#00bfc4"),
         set_name_size = 2, text_size = 2) + 
  ggtitle("Skin") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

pdf("plots/04b-tissue_venns.pdf", height=4, width=7.6)
cowplot::plot_grid(heart_venn, lung_venn, prostate_venn, sk_muscle_venn,
                   E_muscularis_venn, E_mucosa_venn, Breast_venn, Skin_venn,
                   ncol=4)

dev.off()

####################### TISSUE-LEVEL SHARED AND UNIQUE TEs #####################

shared_unique_tissues <-
  function(bulk_tissue, sc_tissue, tissue_name) {
    
    both <- as.data.frame(intersect(bulk_tissue, sc_tissue))
    colnames(both) <- c("TE")
    both$family <- retro.annot$family[match(both$TE, retro.annot$locus)]
    both_only_families <-
      both %>% dplyr::count(family, sort = TRUE) 
    colnames(both_only_families) <- c("family", "both")
    
    bulk_only <- as.data.frame(setdiff(bulk_tissue, sc_tissue))
    colnames(bulk_only) <- c("TE")
    bulk_only$family <- retro.annot$family[match(bulk_only$TE, retro.annot$locus)]
    bulk_only_families <-
      bulk_only %>% dplyr::count(family, sort = TRUE) 
    colnames(bulk_only_families) <- c("family", "bulk")
    
    sc_only <- as.data.frame(setdiff(sc_tissue, bulk_tissue))
    colnames(sc_only) <- c("TE")
    sc_only$family <- retro.annot$family[match(sc_only$TE, retro.annot$locus)]
    sc_only_families <-
      sc_only %>% dplyr::count(family, sort = TRUE) 
    colnames(sc_only_families) <- c("family", "sc")
    
    merged_family_breakdown <- merge(bulk_only_families, 
                                     sc_only_families,
                                     by = "family", all = TRUE) %>%
      merge(both_only_families, by = "family", all = TRUE)
    
    merged_family_breakdown[is.na(merged_family_breakdown)] <- 0
    
    merged_family_breakdown <-
      merged_family_breakdown %>%
      pivot_longer(
        cols = c("bulk", "sc", "both"),
        names_to = "Identified_in",
        values_to = "TE_count"
      )
    
    total_merged_family_breakdown <- 
      merged_family_breakdown %>%
      group_by(Identified_in) %>%
      summarize(total = sum(TE_count))
    
    ggplot(merged_family_breakdown, aes(fill=reorder(family, -TE_count), y=Identified_in, x=TE_count)) + 
      geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
      scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                                   pal_npg("nrc", alpha = 0.7)(10),
                                   pal_jco("default", alpha=0.7)(10),
                                   pal_nejm("default", alpha=0.7)(8),
                                   pal_tron("legacy", alpha=0.7)(7),
                                   pal_lancet("lanonc", alpha=0.7)(9),
                                   pal_startrek("uniform", alpha=0.7)(7))) + 
      coord_flip() +
      theme_pubclean() + 
      guides(fill=guide_legend(title="TE Family")) +
      ylab("TEs identified") +
      xlab("Proportion of TEs") + 
      theme(legend.position = c("right")) + 
      guides(fill = guide_legend(ncol = 2)) +
      ggtitle(tissue_name) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.title=element_blank())

}

pdf("plots/04b-te_families_heart_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Heart_bulk, Heart_sc, "Heart")
dev.off()

pdf("plots/04b-te_families_lung_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Lung_bulk, Lung_sc, "Lung")
dev.off()

pdf("plots/04b-te_families_e_muscularis_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(E_Muscularis_bulk, E_Muscularis_sc, "E. muscolaris")
dev.off()

pdf("plots/04b-te_families_e_mucosa_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(E_Mucosa_bulk, E_Mucosa_sc, "E. mucosa")
dev.off()

pdf("plots/04b-te_families_prostate_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Prostate_bulk, Prostate_sc, "Prostate")
dev.off()

pdf("plots/04b-te_families_breast_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Breast_bulk, Breast_sc, "Breast")
dev.off()

pdf("plots/04b-te_families_skin_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Skin_bulk, Skin_sc, "Skin")
dev.off()

pdf("plots/04b-te_families_sk_muscle_bulk_sc.pdf", height=8, width=5)
shared_unique_tissues(Sk_muscle_bulk, Sk_muscle_sc, "Sk muscle")
dev.off()
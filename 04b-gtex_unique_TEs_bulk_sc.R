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
              label_alpha = 0,
              label_color = "white") + 
  scale_color_brewer(palette = "Paired") 
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
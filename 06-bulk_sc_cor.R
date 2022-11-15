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

# Fix sample row names
row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

# Fix dashes
row.names(pseudobulk.rtx.counts.raw) <- 
  gsub("-", "_", row.names(pseudobulk.rtx.counts.raw))
row.names(pseudobulk.rtx.cpm.raw) <- 
  gsub("-", "_", row.names(pseudobulk.rtx.cpm.raw))

################################# COLOR SETUP ##################################

mypal = pal_npg("nrc", alpha = 0.7)(8)
mycols <- c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2",
            "#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2")

tissues <- c("Prostate", "Heart", "Lung", "E_Muscularis",
             "E_Mucosa" , "Sk_muscle", "Skin", "Breast")

tissue_color<- data.frame(tissues, mycols)

samples$color <- 
  tissue_color$mycols[match(samples$tissue, tissue_color$tissues)]

############################### METADATA SETUP ##################################

# load annotation
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

# Make sure that both counts.rtx and retro.annot have the same loci
retro.annot <- retro.hg38.v1
remove(retro.hg38.v1)
row.names(retro.annot) <- retro.annot$locus

# Load gtf 

gtf <- rtracklayer::import("refs//gencode.v38.annotation.gtf")
gtf_df=as.data.frame(gtf)
gtf_df <-
  gtf_df[,
         c("gene_id", "seqnames", "start", "end", "strand", "width", "gene_name",
           "gene_type")]


colnames(gtf_df) <- c("gene_id", "chrom", "start", "end", "strand", "length",
                      "gene_name", "gene_type")

gene_table <- 
  gtf_df[!duplicated(gtf_df[,c(1,7)]), ] %>% 
  dplyr::select('gene_id', 'gene_name', 'gene_type')

rownames(gene_table) <- gene_table$gene_id

gene_table$display <- gene_table$gene_name
gene_table[duplicated(gene_table$gene_name), 'display'] <- 
  paste(gene_table[duplicated(gene_table$gene_name), 'display'], 
        gene_table[duplicated(gene_table$gene_name), 'gene_id'], sep='|')


############################## GENE DATA SETUP #################################

# Change gene IDs to gene names
rownames(counts.cpm) <- gene_table[rownames(counts.cpm), 'display']

# Keep only the genes, remove TEs
bulk_gene_counts.cpm <- counts.cpm[gene_table$display,]

# Find all common genes between sc and bulk 
common <- intersect(rownames(bulk_gene_counts.cpm), rownames(pseudobulk.cpm.raw))  
bulk_genes_cpm <- counts.cpm[common,] 
pseudobulk_genes_cpm <- pseudobulk.cpm.raw[common,]

############################## HERV DATA SETUP #################################

rownames(pseudobulk.herv.cpm.raw) <- gsub("-", "_", 
                                          rownames(pseudobulk.herv.cpm.raw))

# Find all common genes between sc and bulk 
common_hervs <- intersect(rownames(counts.cpm.herv), 
                          rownames(pseudobulk.herv.cpm.raw))  
bulk_hervs_cpm <- counts.cpm.herv[common_hervs,] 
pseudobulk_hervs_cpm <- pseudobulk.herv.cpm.raw[common_hervs,]

############################## LINE DATA SETUP #################################

rownames(pseudobulk.l1.cpm.raw) <- gsub("-", "_", 
                                          rownames(pseudobulk.l1.cpm.raw))

# Find all common genes between sc and bulk 
common_lines <- intersect(rownames(counts.cpm.l1), 
                          rownames(pseudobulk.l1.cpm.raw))  
bulk_l1_cpm <- counts.cpm.l1[common_lines,] 
pseudobulk_l1_cpm <- pseudobulk.l1.cpm.raw[common_lines,]

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
  
  samp_tes <- merge(bulk_tes, sc_tes, 
                    by.x="row.names", by.y="row.names",
                    all.y=TRUE)
  
  samp_te_df <- samp_tes
  samp_te_df$bulk <- bulk_samp
  samp_te_df$sc <- sc_samp
  samp_te_df$tissue <- samples$tissue[i]
  samp_te_df$participant <- samples$participant_id[i]
  colnames(samp_te_df) <- c("TE", "bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                            "participant")
  row.names(samp_te_df) <- NULL
  
  all_bulk_sc <<- rbind(all_bulk_sc, samp_te_df)
  
  cor.test(samp_tes[,2], samp_tes[,3], 
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

ggsave("plots/sample_level_bulk_sc_te_corr.pdf", height=20, width=20)

########################## SC / BULK CORRELATION ALL TE #########################

ggscatter(all_bulk_sc, x = "bulk_counts", y = "sc_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          exact=FALSE,
          xlab = "Bulk CPM", 
          ylab = "Single Cell CPM") +
  font("xlab", size=10) +
  font("ylab", size=10)

ggsave("plots/all_bulk_sc_te_corr.pdf", height=4, width=4)

############################# SC / BULK CORRELATION ############################
############################### SAMPLE LEVEL GENE ##############################

all_bulk_sc_gene <- data.frame()

bulk_sc_cor_gene <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_genes <- as.data.frame(bulk_genes_cpm[,bulk_samp])
  colnames(bulk_genes) <- c(bulk_samp)
  
  sc_genes <- as.data.frame(pseudobulk_genes_cpm[,sc_samp, drop=FALSE])
  colnames(sc_genes) <- c(sc_samp)
  
  samp_genes <- merge(bulk_genes, sc_genes, 
                      by.x="row.names", by.y="row.names")
  
  samp_genes_df <- samp_genes
  samp_genes_df$bulk <- bulk_samp
  samp_genes_df$sc <- sc_samp
  samp_genes_df$tissue <- samples$tissue[i]
  samp_genes_df$participant <- samples$participant_id[i]
  colnames(samp_genes_df) <- c("gene", "bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                            "participant")
  row.names(samp_genes_df) <- NULL
  
  all_bulk_sc_gene <<- rbind(all_bulk_sc_gene, samp_genes_df)
  
  cor.test(samp_genes[,2], samp_genes[,3], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_genes, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            color=samples$color[samples$bulk_RNAseq == bulk_samp],
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
  
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor_gene), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_gene_corr.pdf", height=20, width=20)

########################## SC / BULK CORRELATION ALL GENES #########################

ggscatter(all_bulk_sc_gene, x = "bulk_counts", y = "sc_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          exact=FALSE,
          xlab = "Bulk Gene CPM", 
          ylab = "Single Cell Gene CPM") +
  font("xlab", size=10) +
  font("ylab", size=10)

ggsave("plots/all_bulk_sc_gene_corr.pdf", height=4, width=4)

############################# SC / BULK CORRELATION ############################
############################### SAMPLE LEVEL HERVs #############################

all_bulk_sc_hervs <- data.frame()

bulk_sc_cor_hervs <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_hervs <- as.data.frame(bulk_hervs_cpm[,bulk_samp])
  colnames(bulk_hervs) <- c(bulk_samp)
  
  sc_hervs <- as.data.frame(pseudobulk_hervs_cpm[,sc_samp, drop=FALSE])
  colnames(sc_hervs) <- c(sc_samp)
  
  samp_hervs <- merge(bulk_hervs, sc_hervs, 
                      by.x="row.names", by.y="row.names")
  
  samp_hervs_df <- samp_hervs
  samp_hervs_df$bulk <- bulk_samp
  samp_hervs_df$sc <- sc_samp
  samp_hervs_df$tissue <- samples$tissue[i]
  samp_hervs_df$participant <- samples$participant_id[i]
  colnames(samp_hervs_df) <- c("herv", "bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                               "participant")
  row.names(samp_hervs_df) <- NULL
  
  all_bulk_sc_hervs <<- rbind(all_bulk_sc_hervs, samp_hervs_df)
  
  cor.test(samp_hervs[,2], samp_hervs[,3], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_hervs, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            color=samples$color[samples$bulk_RNAseq == bulk_samp],
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
  
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor_hervs), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_hervs_corr.pdf", height=20, width=20)

########################## SC / BULK CORRELATION ALL HERVs #########################

ggscatter(all_bulk_sc_hervs, x = "bulk_counts", y = "sc_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          exact=FALSE,
          xlab = "Bulk HERV CPM", 
          ylab = "Single Cell HERV CPM") +
  font("xlab", size=10) +
  font("ylab", size=10)

ggsave("plots/all_bulk_sc_herv_corr.pdf", height=4, width=4)

############################# SC / BULK CORRELATION ############################
############################### SAMPLE LEVEL LINEs #############################

all_bulk_sc_l1 <- data.frame()

bulk_sc_cor_l1 <- function(i) {
  
  bulk_samp <- samples$bulk_RNAseq[i]
  sc_samp <- samples$sn_RNAseq_names[i]
  
  bulk_l1 <- as.data.frame(bulk_l1_cpm[,bulk_samp])
  colnames(bulk_l1) <- c(bulk_samp)
  
  sc_l1 <- as.data.frame(pseudobulk_l1_cpm[,sc_samp, drop=FALSE])
  colnames(sc_l1) <- c(sc_samp)
  
  samp_l1 <- merge(bulk_l1, sc_l1, 
                      by.x="row.names", by.y="row.names")
  
  samp_l1_df <- samp_l1
  samp_l1_df$bulk <- bulk_samp
  samp_l1_df$sc <- sc_samp
  samp_l1_df$tissue <- samples$tissue[i]
  samp_l1_df$participant <- samples$participant_id[i]
  colnames(samp_l1_df) <- c("herv", "bulk_counts", "sc_counts", "bulk", "sc", "tissue",
                               "participant")
  row.names(samp_l1_df) <- NULL
  
  all_bulk_sc_l1 <<- rbind(all_bulk_sc_l1, samp_l1_df)
  
  cor.test(samp_l1[,2], samp_l1[,3], 
           method="spearman", exact = FALSE)
  
  ggscatter(samp_l1, x = bulk_samp, y = sc_samp, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            exact=FALSE,
            color=samples$color[samples$bulk_RNAseq == bulk_samp],
            xlab = paste(bulk_samp, "Bulk", sep = " "), 
            ylab = paste(sc_samp, "SC", sep = " "))
  
}

cowplot::plot_grid(plotlist = lapply(seq_along(samples$bulk_RNAseq), bulk_sc_cor_l1), 
                   ncol=5)

ggsave("plots/sample_level_bulk_sc_l1_corr.pdf", height=20, width=20)

######################## SC / BULK CORRELATION ALL LINEs #######################

ggscatter(all_bulk_sc_l1, x = "bulk_counts", y = "sc_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          exact=FALSE,
          xlab = "Bulk L1 CPM", 
          ylab = "Single Cell L1 CPM") +
  font("xlab", size=10) +
  font("ylab", size=10)

ggsave("plots/all_bulk_sc_l1_corr.pdf", height=4, width=4)

################################### SAVE FILES #################################


save(all_bulk_sc, all_bulk_sc_hervs, all_bulk_sc_l1, all_bulk_sc_gene,
     file="r_outputs/06-bulk_sc_direct_count_comparison.RData")

################################################################################
################################################################################
################################################################################
################################################################################
#############################  GTEx SC Analysis Plan ###########################

# - Make Seurat objects for each sample
# - Do the normalization
# - Make UMAPs
# - Make feature plots of the top features
# - Use the atlas files to get the clusters
# - Start comparing the single cell to bulk
# 
#   Figure 4:
#     - Bar plots of the highest TEs
#     - Top 10 TEs per tissue type, compared to the proportion of the TEs within 
#       cell types in each tissue type (Generate figures for each sample AND for)
#       each tissue. 
#     - What are the top 10 TEs in the single cell sample, and see where they are 
#       in bulk.
#     - Get the percent TE from the bulk and the single cells:
#         - Compare sample by sample.
#         - Compare tissue by tissue.
#     - Is there a difference? 
#     - Does more TE expression mean in general more immune targets?
#     - Top 50 loci in bulk and single-cell pseudobulk, do a union / intersection
#       of both, see the trend. 
#
# - Display the concordance between bulk and single cell:
#     - What about statistically? 
#     - Spearman - rank the genes, raw counts.
#         - How will do they correlate?
#         - Get the residuals to see which ones are detected in one and not the
#           other. (Use the residuals R function, it works on any fit, can plot)
#           the residuals as well. 
# - PCA of 25 pseudobulks vs 25 actual bulks, connecting lines between samples
# - UMAP of pseudobulk vs bulk

#   Things to look up / to do:
#     - Base mean for each gene in the Seurat objects
#     - Get cpm for both pseudobulk and the bulk (same as counts, ranks don't
#       change, but we can now compare sample to sample).
#     - Don't use the Seurat or the DESeq2 normalization for pseudobulk/bulk.
#
#   Comparisons to make
#     - Sample-wise (true pseudobulk): One column for each sample, rows for each 
#       feature.
#     - Compare to single-cell matrix, One column for each cell (all the cells),
#       rows for each feature
#     


################################################################################
################################################################################
################################################################################
################################################################################
#################################### SETUP #####################################

library(tidyr)
library(scopetools)
library(scater)
library(Seurat)
library(edgeR)
library(Seurat)
library(SeuratObject)
library(rtracklayer)

################################### METADATA ####################################
# Read in sample metadata
samples <- read.csv("metadata/GTEx_samples.tsv", sep = "\t")
row.names(samples) <- samples$sn_RNAseq
samples$sn_RNAseq_names <- gsub("-", "_", samples$sn_RNAseq)

# sample names
sample_names <- samples$sn_RNAseq_names

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

gene_table <- 
  rbind(gene_table, data.frame(gene_id=retro.annot$locus, 
                               gene_name=retro.annot$locus, 
                               gene_type=retro.annot$te_class))
rownames(gene_table) <- gene_table$gene_id

gene_table$display <- gene_table$gene_name
gene_table[duplicated(gene_table$gene_name), 'display'] <- 
  paste(gene_table[duplicated(gene_table$gene_name), 'display'], 
        gene_table[duplicated(gene_table$gene_name), 'gene_id'], sep='|')


################################# LOAD SEURAT ##################################

# Function to load Seurat objects for all samples
load_all_seurat <- function(i) {
  sample_name <- samples$sn_RNAseq_names[samples$sn_RNAseq == i]
  
  seurat_object <- 
    scopetools::load_stellarscope_seurat(stellarscope_dir = 
                                         paste("results/stellarscope_pseudobulk/", i, "_rep1_U/", sep = ""),
                                       starsolo_dir = 
                                         paste("results/starsolo_alignment/", i, "/", i, ".Solo.out/Gene/filtered/", sep = ""),
                                       exp_tag = paste(i, "_pseudobulk_U", sep = ""),
                                       use.symbols = TRUE,
                                       TE_metadata = retro.hg38.v1)

  assign(paste(sample_name, "seurat", sep="."), seurat_object, envir=.GlobalEnv)
  }

# Apply function to all samples 
invisible(lapply(samples$sn_RNAseq, load_all_seurat))

# Add seurat obj names to sample sheet
samples$seurat_obj_og <- paste(samples$sn_RNAseq_names, ".seurat", sep="")

############################ REMOVE L1FLnI_Xq21.1db #############################

remove_artifacts <- function(seurat.object) {
  counts <- GetAssayData(seurat.object, assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c('L1FLnI-Xq21.1db'))),]
  seurat.object <<- subset(seurat.object, features = rownames(counts))
}

GTEX_12BJ1_5007_SM_H8L6U.seurat <- remove_artifacts(GTEX_12BJ1_5007_SM_H8L6U.seurat)
GTEX_13N11_5002_SM_H5JDV.seurat <- remove_artifacts(GTEX_13N11_5002_SM_H5JDV.seurat)
GTEX_13N11_5030_SM_H5JDW.seurat <- remove_artifacts(GTEX_13N11_5030_SM_H5JDW.seurat)
GTEX_144GM_5010_SM_HD2M8.seurat <- remove_artifacts(GTEX_144GM_5010_SM_HD2M8.seurat)
GTEX_145ME_5005_SM_H8L6T.seurat <- remove_artifacts(GTEX_145ME_5005_SM_H8L6T.seurat)
GTEX_145ME_5018_SM_G8XQB.seurat <- remove_artifacts(GTEX_145ME_5018_SM_G8XQB.seurat)
GTEX_15CHR_5005_SM_H5JDT.seurat <- remove_artifacts(GTEX_15CHR_5005_SM_H5JDT.seurat)
GTEX_15CHR_5014_SM_H5JDU.seurat <- remove_artifacts(GTEX_15CHR_5014_SM_H5JDU.seurat)
GTEX_15EOM_5003_SM_G64IH.seurat <- remove_artifacts(GTEX_15EOM_5003_SM_G64IH.seurat)
GTEX_15RIE_5015_SM_H8L6X.seurat <- remove_artifacts(GTEX_15RIE_5015_SM_H8L6X.seurat)
GTEX_15RIE_5021_SM_H8L6Y.seurat <- remove_artifacts(GTEX_15RIE_5021_SM_H8L6Y.seurat)
GTEX_15SB6_5008_SM_H8L72.seurat <- remove_artifacts(GTEX_15SB6_5008_SM_H8L72.seurat)
GTEX_16BQI_5013_SM_H8SUW.seurat <- remove_artifacts(GTEX_16BQI_5013_SM_H8SUW.seurat)
GTEX_1CAMR_5015_SM_HPJ3B.seurat <- remove_artifacts(GTEX_1CAMR_5015_SM_HPJ3B.seurat)
GTEX_1CAMS_5015_SM_HPJ3C.seurat <- remove_artifacts(GTEX_1CAMS_5015_SM_HPJ3C.seurat)
GTEX_1HSMQ_5021_SM_HD2MA.seurat <- remove_artifacts(GTEX_1HSMQ_5021_SM_HD2MA.seurat)
GTEX_1HSMQ_5005_SM_GKSJF.seurat <- remove_artifacts(GTEX_1HSMQ_5005_SM_GKSJF.seurat)
GTEX_1HSMQ_5011_SM_GKSJH.seurat <- remove_artifacts(GTEX_1HSMQ_5011_SM_GKSJH.seurat)
GTEX_1HSMQ_5014_SM_GKSJI.seurat <- remove_artifacts(GTEX_1HSMQ_5014_SM_GKSJI.seurat)
GTEX_1HSMQ_5007_SM_GKSJG.seurat <- remove_artifacts(GTEX_1HSMQ_5007_SM_GKSJG.seurat)
GTEX_1I1GU_5006_SM_G8XQC.seurat <- remove_artifacts(GTEX_1I1GU_5006_SM_G8XQC.seurat)
GTEX_1ICG6_5014_SM_GHS9D.seurat <- remove_artifacts(GTEX_1ICG6_5014_SM_GHS9D.seurat)
GTEX_1ICG6_5003_SM_GHS9A.seurat <- remove_artifacts(GTEX_1ICG6_5003_SM_GHS9A.seurat)
GTEX_1MCC2_5013_SM_HPJ3D.seurat <- remove_artifacts(GTEX_1MCC2_5013_SM_HPJ3D.seurat)
GTEX_1R9PN_5002_SM_HD2MC.seurat <- remove_artifacts(GTEX_1R9PN_5002_SM_HD2MC.seurat)

############################### STELLARSCOPE QC ################################

GTEX_12BJ1_5007_SM_H8L6U.seurat.qc <- stellarscope_cell_qc(GTEX_12BJ1_5007_SM_H8L6U.seurat)
GTEX_13N11_5002_SM_H5JDV.seurat.qc <- stellarscope_cell_qc(GTEX_13N11_5002_SM_H5JDV.seurat)
GTEX_13N11_5030_SM_H5JDW.seurat.qc <- stellarscope_cell_qc(GTEX_13N11_5030_SM_H5JDW.seurat)
GTEX_144GM_5010_SM_HD2M8.seurat.qc <- stellarscope_cell_qc(GTEX_144GM_5010_SM_HD2M8.seurat)
GTEX_145ME_5005_SM_H8L6T.seurat.qc <- stellarscope_cell_qc(GTEX_145ME_5005_SM_H8L6T.seurat)
GTEX_145ME_5018_SM_G8XQB.seurat.qc <- stellarscope_cell_qc(GTEX_145ME_5018_SM_G8XQB.seurat)
GTEX_15CHR_5005_SM_H5JDT.seurat.qc <- stellarscope_cell_qc(GTEX_15CHR_5005_SM_H5JDT.seurat)
GTEX_15CHR_5014_SM_H5JDU.seurat.qc <- stellarscope_cell_qc(GTEX_15CHR_5014_SM_H5JDU.seurat)
GTEX_15EOM_5003_SM_G64IH.seurat.qc <- stellarscope_cell_qc(GTEX_15EOM_5003_SM_G64IH.seurat)
GTEX_15RIE_5015_SM_H8L6X.seurat.qc <- stellarscope_cell_qc(GTEX_15RIE_5015_SM_H8L6X.seurat)
GTEX_15RIE_5021_SM_H8L6Y.seurat.qc <- stellarscope_cell_qc(GTEX_15RIE_5021_SM_H8L6Y.seurat)
GTEX_15SB6_5008_SM_H8L72.seurat.qc <- stellarscope_cell_qc(GTEX_15SB6_5008_SM_H8L72.seurat)
GTEX_16BQI_5013_SM_H8SUW.seurat.qc <- stellarscope_cell_qc(GTEX_16BQI_5013_SM_H8SUW.seurat)
GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc <- stellarscope_cell_qc(GTEX_1CAMR_5015_SM_HPJ3B.seurat)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc <- stellarscope_cell_qc(GTEX_1CAMS_5015_SM_HPJ3C.seurat)
GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc <- stellarscope_cell_qc(GTEX_1HSMQ_5021_SM_HD2MA.seurat)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc <- stellarscope_cell_qc(GTEX_1HSMQ_5005_SM_GKSJF.seurat)
GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc <- stellarscope_cell_qc(GTEX_1HSMQ_5011_SM_GKSJH.seurat)
GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc <- stellarscope_cell_qc(GTEX_1HSMQ_5014_SM_GKSJI.seurat)
GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc <- stellarscope_cell_qc(GTEX_1HSMQ_5007_SM_GKSJG.seurat)
GTEX_1I1GU_5006_SM_G8XQC.seurat.qc <- stellarscope_cell_qc(GTEX_1I1GU_5006_SM_G8XQC.seurat)
GTEX_1ICG6_5014_SM_GHS9D.seurat.qc <- stellarscope_cell_qc(GTEX_1ICG6_5014_SM_GHS9D.seurat)
GTEX_1ICG6_5003_SM_GHS9A.seurat.qc <- stellarscope_cell_qc(GTEX_1ICG6_5003_SM_GHS9A.seurat)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc <- stellarscope_cell_qc(GTEX_1MCC2_5013_SM_HPJ3D.seurat)
GTEX_1R9PN_5002_SM_HD2MC.seurat.qc <- stellarscope_cell_qc(GTEX_1R9PN_5002_SM_HD2MC.seurat)



############################# SAVE RDATA OBJECTS ###############################

# Post-QC Seurat Objects

save(GTEX_12BJ1_5007_SM_H8L6U.seurat.qc,GTEX_13N11_5002_SM_H5JDV.seurat.qc,
     GTEX_13N11_5030_SM_H5JDW.seurat.qc, GTEX_144GM_5010_SM_HD2M8.seurat.qc,
     GTEX_145ME_5005_SM_H8L6T.seurat.qc, GTEX_145ME_5018_SM_G8XQB.seurat.qc,
     GTEX_15CHR_5005_SM_H5JDT.seurat.qc, GTEX_15CHR_5014_SM_H5JDU.seurat.qc,
     GTEX_15EOM_5003_SM_G64IH.seurat.qc, GTEX_15RIE_5015_SM_H8L6X.seurat.qc,
     GTEX_15RIE_5021_SM_H8L6Y.seurat.qc, GTEX_15SB6_5008_SM_H8L72.seurat.qc,
     GTEX_16BQI_5013_SM_H8SUW.seurat.qc, GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc,
     GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc, GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc,
     GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc, GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc,
     GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc, GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc,
     GTEX_1I1GU_5006_SM_G8XQC.seurat.qc, GTEX_1ICG6_5014_SM_GHS9D.seurat.qc,
     GTEX_1ICG6_5003_SM_GHS9A.seurat.qc, GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc,
     GTEX_1R9PN_5002_SM_HD2MC.seurat.qc, samples,
     file = "r_outputs/03-gtex_seurat_qc.Rdata")

# Original Seurat Objects

save(GTEX_12BJ1_5007_SM_H8L6U.seurat,GTEX_13N11_5002_SM_H5JDV.seurat,
     GTEX_13N11_5030_SM_H5JDW.seurat, GTEX_144GM_5010_SM_HD2M8.seurat,
     GTEX_145ME_5005_SM_H8L6T.seurat, GTEX_145ME_5018_SM_G8XQB.seurat,
     GTEX_15CHR_5005_SM_H5JDT.seurat, GTEX_15CHR_5014_SM_H5JDU.seurat,
     GTEX_15EOM_5003_SM_G64IH.seurat, GTEX_15RIE_5015_SM_H8L6X.seurat,
     GTEX_15RIE_5021_SM_H8L6Y.seurat, GTEX_15SB6_5008_SM_H8L72.seurat,
     GTEX_16BQI_5013_SM_H8SUW.seurat, GTEX_1CAMR_5015_SM_HPJ3B.seurat,
     GTEX_1CAMS_5015_SM_HPJ3C.seurat, GTEX_1HSMQ_5021_SM_HD2MA.seurat,
     GTEX_1HSMQ_5005_SM_GKSJF.seurat, GTEX_1HSMQ_5011_SM_GKSJH.seurat,
     GTEX_1HSMQ_5014_SM_GKSJI.seurat, GTEX_1HSMQ_5007_SM_GKSJG.seurat,
     GTEX_1I1GU_5006_SM_G8XQC.seurat, GTEX_1ICG6_5014_SM_GHS9D.seurat,
     GTEX_1ICG6_5003_SM_GHS9A.seurat, GTEX_1MCC2_5013_SM_HPJ3D.seurat,
     GTEX_1R9PN_5002_SM_HD2MC.seurat, samples,
     file = "r_outputs/03-gtex_seurat.Rdata")


# Add Seurat object names to sample sheet

samples$seurat_obj_qc <- paste(samples$sn_RNAseq_names, ".seurat.qc", sep="")

############################## GET RAW TE COUNTS ################################

get_feature_counts_sum <- function(seurat.object) {
  
  # get matrix
  seurat.object@assays$RNA@counts %>% Matrix::rowSums() -> features.total.counts
  # keep only TE features
  # features.total.counts <- features.total.counts[!grepl("^ENSG", names(features.total.counts))]
  return(features.total.counts)
}

sc.counts.list <- 
  lapply(list(GTEX_12BJ1_5007_SM_H8L6U.seurat.qc,GTEX_13N11_5002_SM_H5JDV.seurat.qc,
              GTEX_13N11_5030_SM_H5JDW.seurat.qc, GTEX_144GM_5010_SM_HD2M8.seurat.qc,
              GTEX_145ME_5005_SM_H8L6T.seurat.qc, GTEX_145ME_5018_SM_G8XQB.seurat.qc,
              GTEX_15CHR_5005_SM_H5JDT.seurat.qc, GTEX_15CHR_5014_SM_H5JDU.seurat.qc,
              GTEX_15EOM_5003_SM_G64IH.seurat.qc, GTEX_15RIE_5015_SM_H8L6X.seurat.qc,
              GTEX_15RIE_5021_SM_H8L6Y.seurat.qc, GTEX_15SB6_5008_SM_H8L72.seurat.qc,
              GTEX_16BQI_5013_SM_H8SUW.seurat.qc, GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc,
              GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc, GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc,
              GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc, GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc,
              GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc, GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc,
              GTEX_1I1GU_5006_SM_G8XQC.seurat.qc, GTEX_1ICG6_5014_SM_GHS9D.seurat.qc,
              GTEX_1ICG6_5003_SM_GHS9A.seurat.qc, GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc,
              GTEX_1R9PN_5002_SM_HD2MC.seurat.qc),
         get_feature_counts_sum)

names(sc.counts.list) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
                            "GTEX_13N11_5030_SM_H5JDW", "GTEX_144GM_5010_SM_HD2M8",
                            "GTEX_145ME_5005_SM_H8L6T", "GTEX_145ME_5018_SM_G8XQB",
                            "GTEX_15CHR_5005_SM_H5JDT", "GTEX_15CHR_5014_SM_H5JDU",
                            "GTEX_15EOM_5003_SM_G64IH", "GTEX_15RIE_5015_SM_H8L6X",
                            "GTEX_15RIE_5021_SM_H8L6Y", "GTEX_15SB6_5008_SM_H8L72",
                            "GTEX_16BQI_5013_SM_H8SUW", "GTEX_1CAMR_5015_SM_HPJ3B",
                            "GTEX_1CAMS_5015_SM_HPJ3C", "GTEX_1HSMQ_5021_SM_HD2MA",
                            "GTEX_1HSMQ_5005_SM_GKSJF", "GTEX_1HSMQ_5011_SM_GKSJH",
                            "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1HSMQ_5007_SM_GKSJG",
                            "GTEX_1I1GU_5006_SM_G8XQC", "GTEX_1ICG6_5014_SM_GHS9D",
                            "GTEX_1ICG6_5003_SM_GHS9A", "GTEX_1MCC2_5013_SM_HPJ3D",
                            "GTEX_1R9PN_5002_SM_HD2MC")

pseudobulk.counts.raw <- data.frame(sc.counts.list) 
pseudobulk.cpm.raw <- edgeR::cpm(pseudobulk.counts.raw)

########################## GET NORMALIZED TE COUNTS ############################

get_norm_feature_counts_sum <- function(seurat.object) {
  
  # get matrix
  seurat.object@assays$RNA@data %>% Matrix::rowSums() -> features.total.counts
  # keep only TE features
  # features.total.counts <- features.total.counts[!grepl("^ENSG", names(features.total.counts))]
  return(features.total.counts)
}

GTEX_12BJ1_5007_SM_H8L6U.seurat.norm <- NormalizeData(GTEX_12BJ1_5007_SM_H8L6U.seurat.qc)
GTEX_13N11_5002_SM_H5JDV.seurat.norm <- NormalizeData(GTEX_13N11_5002_SM_H5JDV.seurat.qc)
GTEX_13N11_5030_SM_H5JDW.seurat.norm <- NormalizeData(GTEX_13N11_5030_SM_H5JDW.seurat.qc)
GTEX_144GM_5010_SM_HD2M8.seurat.norm <- NormalizeData(GTEX_144GM_5010_SM_HD2M8.seurat.qc)
GTEX_145ME_5005_SM_H8L6T.seurat.norm <- NormalizeData(GTEX_145ME_5005_SM_H8L6T.seurat.qc)
GTEX_145ME_5018_SM_G8XQB.seurat.norm <- NormalizeData(GTEX_145ME_5018_SM_G8XQB.seurat.qc)
GTEX_15CHR_5005_SM_H5JDT.seurat.norm <- NormalizeData(GTEX_15CHR_5005_SM_H5JDT.seurat.qc)
GTEX_15CHR_5014_SM_H5JDU.seurat.norm <- NormalizeData(GTEX_15CHR_5014_SM_H5JDU.seurat.qc)
GTEX_15EOM_5003_SM_G64IH.seurat.norm <- NormalizeData(GTEX_15EOM_5003_SM_G64IH.seurat.qc)
GTEX_15RIE_5015_SM_H8L6X.seurat.norm <- NormalizeData(GTEX_15RIE_5015_SM_H8L6X.seurat.qc)
GTEX_15RIE_5021_SM_H8L6Y.seurat.norm <- NormalizeData(GTEX_15RIE_5021_SM_H8L6Y.seurat.qc)
GTEX_15SB6_5008_SM_H8L72.seurat.norm <- NormalizeData(GTEX_15SB6_5008_SM_H8L72.seurat.qc)
GTEX_16BQI_5013_SM_H8SUW.seurat.norm <- NormalizeData(GTEX_16BQI_5013_SM_H8SUW.seurat.qc)
GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm <- NormalizeData(GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm <- NormalizeData(GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc)
GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm <- NormalizeData(GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm <- NormalizeData(GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc)
GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm <- NormalizeData(GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc)
GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm <- NormalizeData(GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc)
GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm <- NormalizeData(GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc)
GTEX_1I1GU_5006_SM_G8XQC.seurat.norm <- NormalizeData(GTEX_1I1GU_5006_SM_G8XQC.seurat.qc)
GTEX_1ICG6_5014_SM_GHS9D.seurat.norm <- NormalizeData(GTEX_1ICG6_5014_SM_GHS9D.seurat.qc)
GTEX_1ICG6_5003_SM_GHS9A.seurat.norm <- NormalizeData(GTEX_1ICG6_5003_SM_GHS9A.seurat.qc)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm <- NormalizeData(GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc)
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm <- NormalizeData(GTEX_1R9PN_5002_SM_HD2MC.seurat.qc)

sc.counts.norm.list <- 
  lapply(list(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm,GTEX_13N11_5002_SM_H5JDV.seurat.norm,
              GTEX_13N11_5030_SM_H5JDW.seurat.norm, GTEX_144GM_5010_SM_HD2M8.seurat.norm,
              GTEX_145ME_5005_SM_H8L6T.seurat.norm, GTEX_145ME_5018_SM_G8XQB.seurat.norm,
              GTEX_15CHR_5005_SM_H5JDT.seurat.norm, GTEX_15CHR_5014_SM_H5JDU.seurat.norm,
              GTEX_15EOM_5003_SM_G64IH.seurat.norm, GTEX_15RIE_5015_SM_H8L6X.seurat.norm,
              GTEX_15RIE_5021_SM_H8L6Y.seurat.norm, GTEX_15SB6_5008_SM_H8L72.seurat.norm,
              GTEX_16BQI_5013_SM_H8SUW.seurat.norm, GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm,
              GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm,
              GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm,
              GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm, GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm,
              GTEX_1I1GU_5006_SM_G8XQC.seurat.norm, GTEX_1ICG6_5014_SM_GHS9D.seurat.norm,
              GTEX_1ICG6_5003_SM_GHS9A.seurat.norm, GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm,
              GTEX_1R9PN_5002_SM_HD2MC.seurat.norm),
         get_norm_feature_counts_sum)

names(sc.counts.norm.list) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
                            "GTEX_13N11_5030_SM_H5JDW", "GTEX_144GM_5010_SM_HD2M8",
                            "GTEX_145ME_5005_SM_H8L6T", "GTEX_145ME_5018_SM_G8XQB",
                            "GTEX_15CHR_5005_SM_H5JDT", "GTEX_15CHR_5014_SM_H5JDU",
                            "GTEX_15EOM_5003_SM_G64IH", "GTEX_15RIE_5015_SM_H8L6X",
                            "GTEX_15RIE_5021_SM_H8L6Y", "GTEX_15SB6_5008_SM_H8L72",
                            "GTEX_16BQI_5013_SM_H8SUW", "GTEX_1CAMR_5015_SM_HPJ3B",
                            "GTEX_1CAMS_5015_SM_HPJ3C", "GTEX_1HSMQ_5021_SM_HD2MA",
                            "GTEX_1HSMQ_5005_SM_GKSJF", "GTEX_1HSMQ_5011_SM_GKSJH",
                            "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1HSMQ_5007_SM_GKSJG",
                            "GTEX_1I1GU_5006_SM_G8XQC", "GTEX_1ICG6_5014_SM_GHS9D",
                            "GTEX_1ICG6_5003_SM_GHS9A", "GTEX_1MCC2_5013_SM_HPJ3D",
                            "GTEX_1R9PN_5002_SM_HD2MC")

pseudobulk.counts.norm <- data.frame(sc.counts.norm.list) 

########################### GET CPM THROUGH SEURAT #############################

GTEX_12BJ1_5007_SM_H8L6U.seurat.norm.cpm <- NormalizeData(GTEX_12BJ1_5007_SM_H8L6U.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_13N11_5002_SM_H5JDV.seurat.norm.cpm <- NormalizeData(GTEX_13N11_5002_SM_H5JDV.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_13N11_5030_SM_H5JDW.seurat.norm.cpm <- NormalizeData(GTEX_13N11_5030_SM_H5JDW.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_144GM_5010_SM_HD2M8.seurat.norm.cpm <- NormalizeData(GTEX_144GM_5010_SM_HD2M8.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_145ME_5005_SM_H8L6T.seurat.norm.cpm <- NormalizeData(GTEX_145ME_5005_SM_H8L6T.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_145ME_5018_SM_G8XQB.seurat.norm.cpm <- NormalizeData(GTEX_145ME_5018_SM_G8XQB.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15CHR_5005_SM_H5JDT.seurat.norm.cpm <- NormalizeData(GTEX_15CHR_5005_SM_H5JDT.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15CHR_5014_SM_H5JDU.seurat.norm.cpm <- NormalizeData(GTEX_15CHR_5014_SM_H5JDU.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15EOM_5003_SM_G64IH.seurat.norm.cpm <- NormalizeData(GTEX_15EOM_5003_SM_G64IH.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15RIE_5015_SM_H8L6X.seurat.norm.cpm <- NormalizeData(GTEX_15RIE_5015_SM_H8L6X.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15RIE_5021_SM_H8L6Y.seurat.norm.cpm <- NormalizeData(GTEX_15RIE_5021_SM_H8L6Y.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_15SB6_5008_SM_H8L72.seurat.norm.cpm <- NormalizeData(GTEX_15SB6_5008_SM_H8L72.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_16BQI_5013_SM_H8SUW.seurat.norm.cpm <- NormalizeData(GTEX_16BQI_5013_SM_H8SUW.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm.cpm <- NormalizeData(GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm.cpm <- NormalizeData(GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm.cpm <- NormalizeData(GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm.cpm <- NormalizeData(GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm.cpm <- NormalizeData(GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm.cpm <- NormalizeData(GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm.cpm <- NormalizeData(GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1I1GU_5006_SM_G8XQC.seurat.norm.cpm <- NormalizeData(GTEX_1I1GU_5006_SM_G8XQC.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1ICG6_5014_SM_GHS9D.seurat.norm.cpm <- NormalizeData(GTEX_1ICG6_5014_SM_GHS9D.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1ICG6_5003_SM_GHS9A.seurat.norm.cpm <- NormalizeData(GTEX_1ICG6_5003_SM_GHS9A.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm.cpm <- NormalizeData(GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)
GTEX_1R9PN_5002_SM_HD2MC.seurat.norm.cpm <- NormalizeData(GTEX_1R9PN_5002_SM_HD2MC.seurat.qc,
                                                      normalization.method = "RC",
                                                      scale.factor = 1e6)

sc.seurat.cpm.list <- 
  lapply(list(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm.cpm,GTEX_13N11_5002_SM_H5JDV.seurat.norm.cpm,
              GTEX_13N11_5030_SM_H5JDW.seurat.norm.cpm, GTEX_144GM_5010_SM_HD2M8.seurat.norm.cpm,
              GTEX_145ME_5005_SM_H8L6T.seurat.norm.cpm, GTEX_145ME_5018_SM_G8XQB.seurat.norm.cpm,
              GTEX_15CHR_5005_SM_H5JDT.seurat.norm.cpm, GTEX_15CHR_5014_SM_H5JDU.seurat.norm.cpm,
              GTEX_15EOM_5003_SM_G64IH.seurat.norm.cpm, GTEX_15RIE_5015_SM_H8L6X.seurat.norm.cpm,
              GTEX_15RIE_5021_SM_H8L6Y.seurat.norm.cpm, GTEX_15SB6_5008_SM_H8L72.seurat.norm.cpm,
              GTEX_16BQI_5013_SM_H8SUW.seurat.norm.cpm, GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm.cpm,
              GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm.cpm, GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm.cpm,
              GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm.cpm, GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm.cpm,
              GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm.cpm, GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm.cpm,
              GTEX_1I1GU_5006_SM_G8XQC.seurat.norm.cpm, GTEX_1ICG6_5014_SM_GHS9D.seurat.norm.cpm,
              GTEX_1ICG6_5003_SM_GHS9A.seurat.norm.cpm, GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm.cpm,
              GTEX_1R9PN_5002_SM_HD2MC.seurat.norm.cpm),
         get_norm_feature_counts_sum)

names(sc.seurat.cpm.list) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
                                 "GTEX_13N11_5030_SM_H5JDW", "GTEX_144GM_5010_SM_HD2M8",
                                 "GTEX_145ME_5005_SM_H8L6T", "GTEX_145ME_5018_SM_G8XQB",
                                 "GTEX_15CHR_5005_SM_H5JDT", "GTEX_15CHR_5014_SM_H5JDU",
                                 "GTEX_15EOM_5003_SM_G64IH", "GTEX_15RIE_5015_SM_H8L6X",
                                 "GTEX_15RIE_5021_SM_H8L6Y", "GTEX_15SB6_5008_SM_H8L72",
                                 "GTEX_16BQI_5013_SM_H8SUW", "GTEX_1CAMR_5015_SM_HPJ3B",
                                 "GTEX_1CAMS_5015_SM_HPJ3C", "GTEX_1HSMQ_5021_SM_HD2MA",
                                 "GTEX_1HSMQ_5005_SM_GKSJF", "GTEX_1HSMQ_5011_SM_GKSJH",
                                 "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1HSMQ_5007_SM_GKSJG",
                                 "GTEX_1I1GU_5006_SM_G8XQC", "GTEX_1ICG6_5014_SM_GHS9D",
                                 "GTEX_1ICG6_5003_SM_GHS9A", "GTEX_1MCC2_5013_SM_HPJ3D",
                                 "GTEX_1R9PN_5002_SM_HD2MC")

pseudobulk.cpm.seurat <- data.frame(sc.seurat.cpm.list) 

################################# SUBSET TEs ###################################

# Remove all genes, and also remove the weird LINE element L1FLnI_Xq21.1db
feats_in_sc <- as.data.frame(
  GTEX_16BQI_5013_SM_H8SUW.seurat.qc$RNA@meta.features$id[
  !grepl("^ENSG|L1FLnI_Xq21.1db", 
         GTEX_16BQI_5013_SM_H8SUW.seurat.qc$RNA@meta.features$id)])

colnames(feats_in_sc) <- c("tes")
feats_in_sc$tes <- gsub("_", "-", feats_in_sc$tes)

# Make another df of only HERVs and LINEs
meta.features <- as.data.frame(
  GTEX_16BQI_5013_SM_H8SUW.seurat.qc$RNA@meta.features)

meta.features$id <- gsub("_", "-", meta.features$id)

# Subset RTX
pseudobulk.rtx.counts.raw <- pseudobulk.counts.raw[feats_in_sc$tes,]
pseudobulk.rtx.counts.norm <- pseudobulk.counts.norm[feats_in_sc$tes,]
pseudobulk.rtx.cpm.raw <- as.data.frame(pseudobulk.cpm.raw)[feats_in_sc$tes,]
pseudobulk.rtx.cpm.seurat <- as.data.frame(pseudobulk.cpm.seurat)[feats_in_sc$tes,]

# Subset HERVs only

pseudobulk.herv.counts.raw <- pseudobulk.counts.raw[meta.features$id[meta.features$te_class=="LTR"],]
pseudobulk.herv.counts.norm <- pseudobulk.counts.norm[meta.features$id[meta.features$te_class=="LTR"],]
pseudobulk.herv.cpm.raw <- as.data.frame(pseudobulk.cpm.raw)[meta.features$id[meta.features$te_class=="LTR"],]
pseudobulk.herv.cpm.seurat <- as.data.frame(pseudobulk.cpm.seurat)[meta.features$id[meta.features$te_class=="LTR"],]

# Subset LINEs only

pseudobulk.l1.counts.raw <- pseudobulk.counts.raw[meta.features$id[meta.features$te_class=="LINE"],]
pseudobulk.l1.counts.norm <- pseudobulk.counts.norm[meta.features$id[meta.features$te_class=="LINE"],]
pseudobulk.l1.cpm.raw <- as.data.frame(pseudobulk.cpm.raw)[meta.features$id[meta.features$te_class=="LINE"],]
pseudobulk.l1.cpm.seurat <- as.data.frame(pseudobulk.cpm.seurat)[meta.features$id[meta.features$te_class=="LINE"],]

###################### PSEUDOBULK MATRICES THROUGH SEURAT ######################

get_pseudobulk_matrices <- function(seurat.object) {

  # calculate pseudobulk expression
  AggregateExpression(seurat.object,
                      slot = "counts") -> pseudobulk_matrix
  return(pseudobulk_matrix)
}

sc.seurat.pseudobulk.list <- 
  lapply(list(GTEX_12BJ1_5007_SM_H8L6U.seurat.qc,GTEX_13N11_5002_SM_H5JDV.seurat.qc,
              GTEX_13N11_5030_SM_H5JDW.seurat.qc, GTEX_144GM_5010_SM_HD2M8.seurat.qc,
              GTEX_145ME_5005_SM_H8L6T.seurat.qc, GTEX_145ME_5018_SM_G8XQB.seurat.qc,
              GTEX_15CHR_5005_SM_H5JDT.seurat.qc, GTEX_15CHR_5014_SM_H5JDU.seurat.qc,
              GTEX_15EOM_5003_SM_G64IH.seurat.qc, GTEX_15RIE_5015_SM_H8L6X.seurat.qc,
              GTEX_15RIE_5021_SM_H8L6Y.seurat.qc, GTEX_15SB6_5008_SM_H8L72.seurat.qc,
              GTEX_16BQI_5013_SM_H8SUW.seurat.qc, GTEX_1CAMR_5015_SM_HPJ3B.seurat.qc,
              GTEX_1CAMS_5015_SM_HPJ3C.seurat.qc, GTEX_1HSMQ_5021_SM_HD2MA.seurat.qc,
              GTEX_1HSMQ_5005_SM_GKSJF.seurat.qc, GTEX_1HSMQ_5011_SM_GKSJH.seurat.qc,
              GTEX_1HSMQ_5014_SM_GKSJI.seurat.qc, GTEX_1HSMQ_5007_SM_GKSJG.seurat.qc,
              GTEX_1I1GU_5006_SM_G8XQC.seurat.qc, GTEX_1ICG6_5014_SM_GHS9D.seurat.qc,
              GTEX_1ICG6_5003_SM_GHS9A.seurat.qc, GTEX_1MCC2_5013_SM_HPJ3D.seurat.qc,
              GTEX_1R9PN_5002_SM_HD2MC.seurat.qc),
         get_pseudobulk_matrices)

pseudobulk.counts.aggregated.seurat <- data.frame(sc.seurat.pseudobulk.list) 
names(pseudobulk.counts.aggregated.seurat) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
                                      "GTEX_13N11_5030_SM_H5JDW", "GTEX_144GM_5010_SM_HD2M8",
                                      "GTEX_145ME_5005_SM_H8L6T", "GTEX_145ME_5018_SM_G8XQB",
                                      "GTEX_15CHR_5005_SM_H5JDT", "GTEX_15CHR_5014_SM_H5JDU",
                                      "GTEX_15EOM_5003_SM_G64IH", "GTEX_15RIE_5015_SM_H8L6X",
                                      "GTEX_15RIE_5021_SM_H8L6Y", "GTEX_15SB6_5008_SM_H8L72",
                                      "GTEX_16BQI_5013_SM_H8SUW", "GTEX_1CAMR_5015_SM_HPJ3B",
                                      "GTEX_1CAMS_5015_SM_HPJ3C", "GTEX_1HSMQ_5021_SM_HD2MA",
                                      "GTEX_1HSMQ_5005_SM_GKSJF", "GTEX_1HSMQ_5011_SM_GKSJH",
                                      "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1HSMQ_5007_SM_GKSJG",
                                      "GTEX_1I1GU_5006_SM_G8XQC", "GTEX_1ICG6_5014_SM_GHS9D",
                                      "GTEX_1ICG6_5003_SM_GHS9A", "GTEX_1MCC2_5013_SM_HPJ3D",
                                      "GTEX_1R9PN_5002_SM_HD2MC")

pseudobulk.rtx.counts.aggregated.seurat <- 
  pseudobulk.counts.aggregated.seurat[feats_in_sc$tes,]

pseudobulk.rtx.cpm.aggregated.seurat <- edgeR::cpm(pseudobulk.rtx.counts.aggregated.seurat)

# Error in cpm.default(pseudobulk.rtx.counts.aggregated.seurat) : 
#  library sizes should be finite and non-negative

###############################################################################
## Get average expression of the normalized data of every feature per tissue ##
###############################################################################

# Function to create df of counts and cpm for a given tissue type
get_tissue_samples <- function(i) {
  tissue_sample_names <- samples$sn_RNAseq_names[samples$tissue == i]
  
  df_raw <- as.data.frame(rowMeans(pseudobulk.rtx.counts.raw[, tissue_sample_names])) 
  colnames(df_raw) <- c("mean_counts")
  df_raw <- df_raw[order(-df_raw$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "raw", sep="."), df_raw, envir=.GlobalEnv)
  
  df_raw_herv <- as.data.frame(rowMeans(pseudobulk.herv.counts.raw[, tissue_sample_names])) 
  colnames(df_raw_herv) <- c("mean_counts")
  df_raw_herv <- df_raw_herv[order(-df_raw_herv$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "herv", "raw", sep="."), df_raw_herv, envir=.GlobalEnv)
  
  df_raw_l1 <- as.data.frame(rowMeans(pseudobulk.l1.counts.raw[, tissue_sample_names])) 
  colnames(df_raw_l1) <- c("mean_counts")
  df_raw_l1 <- df_raw_l1[order(-df_raw_l1$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "l1", "raw", sep="."), df_raw_l1, envir=.GlobalEnv)
  
  df_cpm <- as.data.frame(rowMeans(pseudobulk.rtx.cpm.raw[, tissue_sample_names])) 
  colnames(df_cpm) <- c("mean_counts")
  df_cpm <- df_cpm[order(-df_cpm$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "cpm", sep="."), df_cpm, envir=.GlobalEnv)
  
  df_cpm_herv <- as.data.frame(rowMeans(pseudobulk.herv.cpm.raw[, tissue_sample_names])) 
  colnames(df_cpm_herv) <- c("mean_counts")
  df_cpm_herv <- df_cpm_herv[order(-df_cpm_herv$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "herv", "cpm", sep="."), df_cpm_herv, envir=.GlobalEnv)
  
  df_cpm_l1 <- as.data.frame(rowMeans(pseudobulk.l1.cpm.raw[, tissue_sample_names])) 
  colnames(df_cpm_l1) <- c("mean_counts")
  df_cpm_l1 <- df_cpm_l1[order(-df_cpm_l1$mean_counts), , drop = FALSE]
  assign(paste(i, "sc", "l1", "cpm", sep="."), df_cpm_l1, envir=.GlobalEnv)
  
  remove(df_raw, df_cpm, tissue_sample_names)
}

# Apply function to all tissue types
invisible(lapply(unique(samples$tissue), get_tissue_samples))

################################# SAVE FILES ###################################

save(GTEX_12BJ1_5007_SM_H8L6U.seurat.norm,GTEX_13N11_5002_SM_H5JDV.seurat.norm,
     GTEX_13N11_5030_SM_H5JDW.seurat.norm, GTEX_144GM_5010_SM_HD2M8.seurat.norm,
     GTEX_145ME_5005_SM_H8L6T.seurat.norm, GTEX_145ME_5018_SM_G8XQB.seurat.norm,
     GTEX_15CHR_5005_SM_H5JDT.seurat.norm, GTEX_15CHR_5014_SM_H5JDU.seurat.norm,
     GTEX_15EOM_5003_SM_G64IH.seurat.norm, GTEX_15RIE_5015_SM_H8L6X.seurat.norm,
     GTEX_15RIE_5021_SM_H8L6Y.seurat.norm, GTEX_15SB6_5008_SM_H8L72.seurat.norm,
     GTEX_16BQI_5013_SM_H8SUW.seurat.norm, GTEX_1CAMR_5015_SM_HPJ3B.seurat.norm,
     GTEX_1CAMS_5015_SM_HPJ3C.seurat.norm, GTEX_1HSMQ_5021_SM_HD2MA.seurat.norm,
     GTEX_1HSMQ_5005_SM_GKSJF.seurat.norm, GTEX_1HSMQ_5011_SM_GKSJH.seurat.norm,
     GTEX_1HSMQ_5014_SM_GKSJI.seurat.norm, GTEX_1HSMQ_5007_SM_GKSJG.seurat.norm,
     GTEX_1I1GU_5006_SM_G8XQC.seurat.norm, GTEX_1ICG6_5014_SM_GHS9D.seurat.norm,
     GTEX_1ICG6_5003_SM_GHS9A.seurat.norm, GTEX_1MCC2_5013_SM_HPJ3D.seurat.norm,
     GTEX_1R9PN_5002_SM_HD2MC.seurat.norm, samples,
     file = "r_outputs/03-gtex_seurat.norm.Rdata")

# Save count files
save(pseudobulk.rtx.counts.raw, pseudobulk.rtx.counts.norm, 
     pseudobulk.rtx.cpm.raw, pseudobulk.counts.raw, pseudobulk.counts.norm,
     pseudobulk.cpm.raw, pseudobulk.herv.counts.raw, pseudobulk.herv.cpm.raw,
     pseudobulk.l1.counts.raw, pseudobulk.l1.cpm.raw,
     file = "r_outputs/03-scgtex_seurat_counts.Rdata")

# Save all ordered &  mean raw TE counts per tissue sample
save(pseudobulk.rtx.counts.raw, pseudobulk.herv.counts.raw, pseudobulk.l1.counts.raw,
     Breast.sc.raw, E_Mucosa.sc.raw, E_Muscularis.sc.raw, 
     Heart.sc.raw, Lung.sc.raw, Prostate.sc.raw, Sk_muscle.sc.raw, Skin.sc.raw,
     Breast.sc.herv.raw, E_Mucosa.sc.herv.raw, E_Muscularis.sc.herv.raw, 
     Heart.sc.herv.raw, Lung.sc.herv.raw, Prostate.sc.herv.raw, Sk_muscle.sc.herv.raw,
     Skin.sc.herv.raw,
     file="r_outputs/03-mean_raw_scTE_counts_by_tissue_type.RData")

# Save all ordered &  mean raw TE counts per tissue sample
save(pseudobulk.rtx.cpm.raw, pseudobulk.herv.cpm.raw, pseudobulk.l1.cpm.raw,
     Breast.sc.cpm, E_Mucosa.sc.cpm, E_Muscularis.sc.cpm, 
     Heart.sc.cpm, Lung.sc.cpm, Prostate.sc.cpm, Sk_muscle.sc.cpm, Skin.sc.cpm,
     Breast.sc.herv.cpm, E_Mucosa.sc.herv.cpm, E_Muscularis.sc.herv.cpm,
     Heart.sc.herv.cpm, Lung.sc.herv.cpm, Prostate.sc.herv.cpm, Sk_muscle.sc.herv.cpm,
     Skin.sc.herv.cpm, Breast.sc.l1.cpm, E_Mucosa.sc.l1.cpm, E_Muscularis.sc.l1.cpm,
     Heart.sc.l1.cpm, Lung.sc.l1.cpm, Prostate.sc.l1.cpm, Sk_muscle.sc.l1.cpm,
     Skin.sc.l1.cpm,
     file="r_outputs/03-mean_raw_scTE_cpm_by_tissue_type.RData")



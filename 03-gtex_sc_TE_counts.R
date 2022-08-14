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
#       cell types in each tissue type
#     - Look sample by sample
#     - What are the top 10 TEs in the single cell sample, and see where they are 
#       in bulk
#     - Get the percent TE from the bulk and the single cells:
#     - compare sample by sample
#     - Compare tissue by tissue
#     - Is there a difference? 
#     - Does more TE expression mean in general more immune targets?
#
# - Display the concordance between bulk and single cell:
#     - What about statistically? 
#     - Spearman - rank the genes, raw counts.
# - PCA of 25 pseudobulks vs 25 actual bulks, connecting lines between samples
# - UMAP of pseudobulk vs bulk

################################################################################
################################################################################
################################################################################
################################################################################
#################################### SETUP #####################################

library(tidyr)
library(scopetools)
library(scater)
library(Seurat)

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

################################# LOAD SEURAT ##################################

# Function to load Seurat objects for all samples
load_all_seurat <- function(i) {
  sample_name <- samples$sn_RNAseq_names[samples$sn_RNAseq == i]
  
  seurat_object <- 
    scopetools::load_stellarscope_seurat(stellarscope_dir = 
                                         paste("results/stellarscope_pseudobulk/", i, "/", sep = ""),
                                       TE_count_file = 
                                         paste("results/stellarscope_pseudobulk/", i, "/", i, "_pseudobulk-TE_counts_exclude.mtx", sep=""),
                                       starsolo_dir = 
                                         paste("results/starsolo_alignment/", i, "/", i, ".Solo.out/Gene/filtered/", sep = ""))
  
  assign(paste(sample_name, "seurat", sep="."), seurat_object, envir=.GlobalEnv)
  }

# Apply function to all samples 
invisible(lapply(samples$sn_RNAseq, load_all_seurat))

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

te.total.counts <- 
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

names(te.total.counts) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
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

pseudobulk_te.counts.raw <- data.frame(te.total.counts) 

########################## GET NORMALIZED TE COUNTS ############################

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

te.total.counts.norm <- 
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
         get_feature_counts_sum)

names(te.total.counts.norm) <- c("GTEX_12BJ1_5007_SM_H8L6U","GTEX_13N11_5002_SM_H5JDV",
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

pseudobulk_te.counts.norm <- data.frame(te.total.counts.norm) 

# Original Seurat Objects

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
save(pseudobulk_te.counts.norm, pseudobulk_te.counts.raw,
     file = "r_outputs/03-scgtex_seurat_counts.Rdata")


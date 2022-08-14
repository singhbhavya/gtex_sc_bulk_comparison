################################################################################
################################################################################
#################################### SETUP #####################################

library(reticulate)
library(anndata)
library(zellkonverter)
library(Seurat)

############################# CONVERT ATLAS FILES ##############################

GTEx_8_snRNAseq_atlas <-
  readH5AD(file = "literature/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")

GTEx_8_snRNAseq_immune_atlas <-
  readH5AD(file = "literature/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad")

GTEx_8_snRNAseq_atlas.seurat <- as.Seurat(GTEx_8_snRNAseq_atlas, data = NULL)

GTEx_8_snRNAseq_immune_atlas.seurat <- as.Seurat(GTEx_8_snRNAseq_immune_atlas, data = NULL)

save(GTEx_8_snRNAseq_atlas.seurat, file="r_outputs/02a-GTEx_8_snRNAseq_atlas.seurat.RData")
save(GTEx_8_snRNAseq_immune_atlas.seurat, file="r_outputs/02a-GTEx_8_snRNAseq_immune_atlas.seurat.RData")
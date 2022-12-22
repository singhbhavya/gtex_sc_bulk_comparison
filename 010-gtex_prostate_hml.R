################################################################################
################################################################################
################################################################################
################################################################################
###########################  GTEx SC PROSTATE (SEURAT) #########################

#################################### SETUP #####################################

library(tidyverse)
library(Matrix)
library(scater)
library(Seurat)
library(SeuratData)
library(edgeR)
library(SeuratObject)
library(rtracklayer)
library(Azimuth)
library(patchwork)
library(cowplot)
library(data.table)

################################## LOAD DATA ###################################

load("r_outputs/09-prostate_merged.RData")
load("r_outputs/01-annot.Rdata")

################################## AGGREGATED EXPRESSION ###################################

aggregated_expression <-
  AggregateExpression(
    prostate.norm.merged,
    group.by = "ident",
    slot = "counts"
  )


avg_expression <-
  AverageExpression(prostate.norm.merged,
                    group.by = "ident",
                    slot = "counts",
                    verbose = TRUE)


################################################################################
################################################################################
################################################################################
################################################################################
#############################  GTEx SC HEART (SEURAT) ##########################

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

load("r_outputs/03-gtex_seurat.norm.Rdata")


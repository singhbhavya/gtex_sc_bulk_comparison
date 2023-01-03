################################################################################
################################################################################
################################################################################
################################################################################
#############################  GTEx SC HEART (SEURAT) ##########################

#################################### SETUP #####################################
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
library(ggsci)
library(pheatmap)
library(grid) 

################################## LOAD DATA ###################################

load("r_outputs/08b-heart_merged.RData")
load("r_outputs/01-annot.Rdata")
load("r_outputs/01-counts.Rdata")
load("r_outputs/03-scgtex_seurat_counts.Rdata")

################################ METADATA SETUP ################################

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

metadata.heart <- metadata[c("GTEX-13N11-0426-SM-5KM3O", "GTEX-15RIE-1726-SM-7KUMU",
                                "GTEX-1ICG6-0626-SM-ACKWS", "GTEX_13N11_5002_SM_H5JDV",
                                "GTEX_15RIE_5015_SM_H8L6X", "GTEX_1ICG6_5003_SM_GHS9A"),]

metadata.heart <- metadata.heart[,c("participant_id", "type")]

################################ PHEATMAP SETUP ################################

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), 
                 vjust = 0.5,
                 hjust = 1, 
                 rot = 75, 
                 gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

############################ AGGREGATED EXPRESSION #############################

aggregated_expression <-
  AggregateExpression(
    heart.norm.merged,
    group.by = "ident",
    slot = "counts"
  )


avg_expression <-
  AverageExpression(heart.norm.merged,
                    group.by = "ident",
                    slot = "counts",
                    verbose = TRUE)


############################### HML-2 SC AVERAGE ###############################

hml2 <- retro.annot[retro.annot$family == "HML2",]
rownames(hml2) <- hml2$locus
rownames(hml2) <- gsub("_", "-", rownames(hml2) )

avg_expression_df <- as.data.frame(avg_expression$RNA)

avg_expression_hml <- avg_expression_df[rownames(hml2),]
avg_expression_hml <- na.omit(avg_expression_hml)
avg_expression_hml <- avg_expression_hml[rowSums
                                         (avg_expression_hml[])>0,]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

avg_exp_z_scores <- t(apply(avg_expression_hml, 1, cal_z_score))

pdf("plots/011_heart_hml2_clusters_avg.pdf", height=12, width=5)
pheatmap(avg_expression_hml,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

pdf("plots/011_heart_hml2_clusters_avg_z_scored.pdf", height=12, width=5)
pheatmap(avg_exp_z_scores,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

############################## HML-2 SC AGGREGATE ##############################

aggregated_expression_df <- as.data.frame(aggregated_expression$RNA)

aggregated_expression_df <- aggregated_expression_df[rownames(hml2),]
aggregated_expression_hml <- na.omit(aggregated_expression_df)
aggregated_expression_hml <- aggregated_expression_hml[rowSums
                                                       (aggregated_expression_hml[])>0,]

agg_exp_z_scores <- t(apply(aggregated_expression_hml, 1, cal_z_score))

pheatmap(agg_exp_z_scores)
pdf("plots/011_heart_hml2_clusters_agg.pdf", height=12, width=5)
pheatmap(agg_exp_z_scores,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

########################### HML-2 BULK & PSEUDO-BULK ###########################
###################################### CPM #####################################

# Merge cpm
row.names(pseudobulk.herv.cpm.raw) <- gsub("-", "_", row.names(pseudobulk.herv.cpm.raw))

combined_cpm <- merge(counts.cpm.herv, pseudobulk.herv.cpm.raw, 
                      by.x="row.names", by.y="row.names",
                      all.y=TRUE)

rownames(combined_cpm) <- combined_cpm[,1]
combined_cpm[,1] <- NULL

stopifnot(all(colnames(combined_cpm) == rownames(metadata)))

thresh <- combined_cpm > 0.5
# This produces a logical matrix with TRUEs and FALSE
head(thresh)
colSums(head(thresh))

keep <- rowSums(thresh) >= 1
summary(keep)

cpm.keep <- combined_cpm[keep,]
cpm.heart <- cpm.keep[, c("GTEX-13N11-0426-SM-5KM3O", "GTEX-15RIE-1726-SM-7KUMU",
                          "GTEX-1ICG6-0626-SM-ACKWS", "GTEX_13N11_5002_SM_H5JDV",
                          "GTEX_15RIE_5015_SM_H8L6X", "GTEX_1ICG6_5003_SM_GHS9A")]

hml2_underscores <- hml2
rownames(hml2_underscores) <- gsub("-", "_", rownames(hml2_underscores) )
cpm.heart <- cpm.heart[rownames(hml2_underscores),]
cpm.heart <- na.omit(cpm.heart)
cpm.heart <- cpm.heart[rowSums
                             (cpm.heart[])>0,]

cpm_prostate_z_scores <- t(apply(cpm.heart, 1, cal_z_score))

pheatmap(cpm_prostate_z_scores)
pdf("plots/011_heart_hml2_pseudobulk_cpm.pdf", height=13, width=5)
pheatmap(cpm.heart,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         annotation_col = metadata.heart,
         cluster_cols = FALSE,
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()


pdf("plots/011_heart_hml2_pseudobulk_cpm_z_scored.pdf", height=13, width=5)
pheatmap(cpm_prostate_z_scores,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         annotation_col = metadata.heart,
         cluster_cols = FALSE,
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

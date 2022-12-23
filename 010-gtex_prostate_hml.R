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
library(ggsci)
library(pheatmap)
library(grid) 

################################## LOAD DATA ###################################

load("r_outputs/09-prostate_merged.RData")
load("r_outputs/01-annot.Rdata")
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

metadata.prostate <- metadata[c("GTEX-12BJ1-1226-SM-5LUAE", "GTEX-15CHR-1226-SM-79OON",
                          "GTEX-1HSMQ-1826-SM-A9SMJ", "GTEX-1I1GU-1126-SM-A96S9",
                          "GTEX_12BJ1_5007_SM_H8L6U", "GTEX_15CHR_5014_SM_H5JDU",
                          "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1I1GU_5006_SM_G8XQC"),]

metadata.prostate <- metadata[,c("participant_id", "type")]

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
    prostate.norm.merged,
    group.by = "ident",
    slot = "counts"
  )


avg_expression <-
  AverageExpression(prostate.norm.merged,
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

pdf("plots/010_prostate_hml2_clusters_avg.pdf", height=12, width=5)
pheatmap(avg_expression_hml,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

pdf("plots/010_prostate_hml2_clusters_avg_z_scored.pdf", height=12, width=5)
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
pdf("plots/010_prostate_hml2_clusters_agg.pdf", height=12, width=5)
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
cpm.prostate <- cpm.keep[, c("GTEX-12BJ1-1226-SM-5LUAE", "GTEX-15CHR-1226-SM-79OON",
                             "GTEX-1HSMQ-1826-SM-A9SMJ", "GTEX-1I1GU-1126-SM-A96S9",
                             "GTEX_12BJ1_5007_SM_H8L6U", "GTEX_15CHR_5014_SM_H5JDU",
                             "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1I1GU_5006_SM_G8XQC")]

hml2_underscores <- hml2
rownames(hml2_underscores) <- gsub("-", "_", rownames(hml2_underscores) )
cpm.prostate <- cpm.prostate[rownames(hml2_underscores),]
cpm.prostate <- na.omit(cpm.prostate)
cpm.prostate <- cpm.prostate[rowSums
                             (cpm.prostate[])>0,]

cpm_prostate_z_scores <- t(apply(cpm.prostate, 1, cal_z_score))

pheatmap(cpm_prostate_z_scores)
pdf("plots/010_prostate_hml2_pseudobulk_cpm.pdf", height=13, width=5)
pheatmap(cpm.prostate,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         annotation_col = metadata.prostate,
         cluster_cols = FALSE,
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()


pdf("plots/010_prostate_hml2_pseudobulk_cpm_z_scored.pdf", height=13, width=5)
pheatmap(cpm_prostate_z_scores,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         annotation_col = metadata.prostate,
         cluster_cols = FALSE,
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()


#################################### COUNTS ####################################

row.names(pseudobulk.herv.counts.raw) <- gsub("-", "_", row.names(pseudobulk.herv.counts.raw))

combined_counts <- merge(counts.herv, pseudobulk.herv.counts.raw, 
                      by.x="row.names", by.y="row.names",
                      all.y=TRUE)

rownames(combined_counts) <- combined_counts[,1]
combined_counts[,1] <- NULL

stopifnot(all(colnames(combined_counts) == rownames(metadata)))

thresh <- combined_counts > 2
# This produces a logical matrix with TRUEs and FALSE
head(thresh)
colSums(head(thresh))

keep <- rowSums(thresh) >= 5
summary(keep)

counts.keep <- combined_counts[keep,]
counts.prostate <- counts.keep[, c("GTEX-12BJ1-1226-SM-5LUAE", "GTEX-15CHR-1226-SM-79OON",
                             "GTEX-1HSMQ-1826-SM-A9SMJ", "GTEX-1I1GU-1126-SM-A96S9",
                             "GTEX_12BJ1_5007_SM_H8L6U", "GTEX_15CHR_5014_SM_H5JDU",
                             "GTEX_1HSMQ_5014_SM_GKSJI", "GTEX_1I1GU_5006_SM_G8XQC")]

counts.prostate <- counts.prostate[rownames(hml2_underscores),]
counts.prostate <- na.omit(counts.prostate)
counts.prostate <- counts.prostate[rowSums
                             (counts.prostate[])>0,]

counts_prostate_z_scores <- t(apply(counts.prostate, 1, cal_z_score))

pheatmap(counts_prostate_z_scores)
pdf("plots/010_prostate_hml2_pseudobulk_counts.pdf", height=13, width=5)
pheatmap(counts_prostate_z_scores,
         color = cols,
         breaks=seq(-3,3,length.out=14), 
         annotation_col = metadata.prostate,
         cluster_cols = FALSE,
         cexRow=1,
         cexCol=1,
         trace="none")
dev.off()

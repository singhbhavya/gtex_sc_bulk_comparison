################################################################################
################################################################################
################################################################################
################################################################################
############################  GTEx Bulk Analysis Plan ##########################

# -   For each tissue site, we have varying number of subjects.
# -   Each subject has one single-cell sample, and one corresponding bulk.

# **To do:**
#   
# -   Use scope tools to pull in all the reports (count matrices for the TEs 
#     and the genes)
# -   Concatenate (rbind) the matrices
# -   We will have 27 to 25 columns, and roughly 80,0000 features in that matrix 
#     of counts
# -   Run DESEq2 on all of those together to normalize them (design will just be 
#     \~0 or a \~1, since this is an unsupervised matrix, and we are just trying 
#     to normalize the counts).
# -   Within each tissue, everything can be compared to everything.
# -   Subset the breast set, and do a rowSum / rowMeans to get the average 
#     expression of the normalized data of every feature. We'll have one 
#     normalized mean value for every feature within breast. And then, just rank 
#     those so that we know which TEs are expressed.
# -   Now we know what is expressed / not expressed. We'll set a threshold of 5 
#     for the normalized mean count?
# -   Of the ones that are expressed in bulk, do we see them in single cell, and 
#     what is the pattern?
#   
#   **Comparison:**
#   
#   -   Go down the list of ranked TEs that are expressed in the tissues. Are 
#       these TEs expressed in at least some of the cells in single cell?
#   -   We want to start looking at the highest expressed ones. Is it highly 
#       expressed also in single cell?
#   -   It will be interesting if there are things that are missed
#   -   if it is seen in bulk, why is it missed in single cell?
#   -   Is this actually calling everything at the single cell level?
#   
#   **Sample by sample comparison:**
#   
# -   In addition to looking at the mean of the breast and prostate, we can also 
#     do a sample to sample comparison.
# -   Do the sample process where everything is normalized together.
# -   We could do this on a count level (since we're just using ranks), on a 
#     sample-by-sample basis. Direct comparison of bulk samples to single cell 
#     samples.
# 
# **At the end, save an RData object of the following:**
# 
# -   Raw counts
# -   Normalized counts
# -   Averaged by tissue type
# 
# **Next analysis:**
# 
# -   Once I've started looking at the bulk, the next thing is to explore the h5 
#     files to find the published cell types

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#################################### SETUP #####################################

library(tidyr)
library(scopetools)
library(DESeq2)
library(edgeR)

################################### METADATA ####################################
# Read in sample metadata
samples <- read.csv("metadata/GTEx_samples.tsv", sep = "\t")
row.names(samples) <- samples$bulk_RNAseq

# sample names
sample_names <- samples$bulk_RNAseq

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

# Annotation directory for scopetools
ddir <- system.file("extdata", package="scopetools")

# Remove the confounding LINE element (L1FLnI_Xq21.1db) that has a poly A tail
# in the middle of it:

retro.hg38.v1<-
  retro.hg38.v1[!(retro.hg38.v1$locus=="L1FLnI_Xq21.1db"),]

################################# READ TELESCOPE #################################
# Telescope files to read in 
# Provide all_locs to include all annotated loci (preferred)
t_files <- file.path("results/telescope", 
                     paste0(sample_names, 
                            '/',
                            sample_names,
                            '_telescope.report.tsv'))
names(t_files) <- samples$bulk_RNAseq

# Load TE counts
counts.rtx <- load_telescope_reports(t_files, 
                                     all_locs=retro.hg38.v1$locus, 
                                     count_column = "count")

# remove unused variables
remove(t_files, ddir)

################################### READ STAR ##################################

# Create list of files to import, and import in bulk
star_files <- list.files( path = "results/star_alignment/bulk", 
                          pattern = "*ReadsPerGene.out.tab$", 
                          full.names = TRUE,
                          recursive = TRUE)
counts.files <- lapply( star_files, read.table, skip = 4 )

# Import second column (counts for unstranded data)
counts.tx <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )

# Clean up sample names
star_files <- gsub( "results/star_alignment/bulk/", "", star_files)
star_files <- gsub( "_GDC38.ReadsPerGene.out.tab", "", star_files)
star_files <-  gsub("(.*/\\s*(.*$))", "\\2", star_files)
colnames(counts.tx) <- star_files

row.names(counts.tx) <- counts.files[[1]]$V1

# remove unused variables
remove(counts.files, star_files, sample_names)

################################# Sanity Check #################################

# Make sure that both counts.tx and counts.rtx have the same samples
counts.tx <- counts.tx %>% dplyr::select(order(colnames(counts.tx)))
counts.rtx <- counts.rtx %>% dplyr::select(order(colnames(counts.rtx)))
stopifnot(all(names(counts.tx) == names(counts.rtx)))

# Make sure that both counts.rtx and retro.annot have the same loci
retro.annot <- retro.hg38.v1
remove(retro.hg38.v1)
row.names(retro.annot) <- retro.annot$locus
stopifnot(all(rownames(counts.rtx) == rownames(retro.annot)))

# reorder counts.tx by metadata rowname
reorder_idx_counts.tx <- match(rownames(samples), colnames(counts.tx))
counts.tx <- counts.tx[,reorder_idx_counts.tx]

# reorder counts.rts by metadata rowname
reorder_idx_counts.rtx <- match(rownames(samples), colnames(counts.rtx))
counts.rtx <- counts.rtx[ , reorder_idx_counts.rtx]

# sanity check
stopifnot(all(names(counts.tx) == names(counts.rtx)))

# Combine counts
counts <- rbind(counts.tx, counts.rtx)

# sanity check
stopifnot(all(names(counts) == row.names(samples)))

# Remove unused vairbles
remove(reorder_idx_counts.tx, reorder_idx_counts.rtx)

################################ REMOVE LINEs ##################################

herv.annot <- retro.annot[retro.annot$te_class == "LTR",]
counts.herv <-counts.rtx[herv.annot$locus,]

##################################### CPM ######################################

counts.cpm <- edgeR::cpm(counts)

#################################### DESeq2 ####################################

dds <- DESeq2::DESeqDataSetFromMatrix(
  counts, samples, ~1)

dds <- DESeq2::DESeq(dds, parallel=T)
dds <- estimateSizeFactors(dds)

counts.norm <- counts(dds, normalized=TRUE)

boxplot(counts.norm, main = "Normalized Counts", cex = 0.6)
boxplot(counts(dds), main = "Raw Counts", cex = 0.6)

###############################################################################
## Get average expression of the normalized data of every feature per tissue ##
###############################################################################

# Subset raw counts to only include HERVs
counts.herv <- counts[herv.annot$locus,]

# Subset raw counts to only include L1s
counts.l1 <- counts[retro.annot$locus[retro.annot$te_class=="LINE"],]

# Subset normalized counts to only include TEs
counts.norm.rtx <- counts.norm[row.names(retro.annot),]

# Subset normalized counts to only include HERVs
counts.norm.herv <- counts.norm[herv.annot$locus,]

# Subset normalized counts to only include L1s
counts.norm.l1 <- counts.norm[retro.annot$locus[retro.annot$te_class=="LINE"],]

# Subset cpm to only include TEs
counts.cpm.rtx <- counts.cpm[row.names(retro.annot),]

# Subset cpm to only include HERVs
counts.cpm.herv <- counts.cpm[herv.annot$locus,]

# Subset cpm to only include L1s
counts.cpm.l1 <- counts.cpm[retro.annot$locus[retro.annot$te_class=="LINE"],]

# Function to create df of normalized count for a given tissue type
get_tissue_samples <- function(i) {
  tissue_sample_names <- samples$bulk_RNAseq[samples$tissue == i]
  
  df <- as.data.frame(rowMeans(counts.norm.rtx[, tissue_sample_names])) 
  colnames(df) <- c("mean_counts")
  df <- df[order(-df$mean_counts), , drop = FALSE]
  assign(paste(i, "norm", sep="."), df, envir=.GlobalEnv)
  
  df_raw <- as.data.frame(rowMeans(counts.rtx[, tissue_sample_names])) 
  colnames(df_raw) <- c("mean_counts")
  df_raw <- df_raw[order(-df_raw$mean_counts), , drop = FALSE]
  assign(paste(i, "raw", sep="."), df_raw, envir=.GlobalEnv)
  
  df_cpm <- as.data.frame(rowMeans(counts.cpm.rtx[, tissue_sample_names])) 
  colnames(df_cpm) <- c("mean_counts")
  df_cpm <- df_cpm[order(-df_cpm$mean_counts), , drop = FALSE]
  assign(paste(i, "cpm", sep="."), df_cpm, envir=.GlobalEnv)
  
  df_cpm_l1 <- as.data.frame(rowMeans(counts.cpm.l1[, tissue_sample_names])) 
  colnames(df_cpm_l1) <- c("mean_counts")
  df_cpm_l1 <- df_cpm_l1[order(-df_cpm_l1$mean_counts), , drop = FALSE]
  assign(paste(i, "l1", "cpm", sep="."), df_cpm_l1, envir=.GlobalEnv)
  
  df_cpm_herv <- as.data.frame(rowMeans(counts.cpm.herv[, tissue_sample_names])) 
  colnames(df_cpm_herv) <- c("mean_counts")
  df_cpm_herv <- df_cpm_herv[order(-df_cpm_herv$mean_counts), , drop = FALSE]
  assign(paste(i, "herv", "cpm", sep="."), df_cpm_herv, envir=.GlobalEnv)
  
  remove(df, df_raw, df_cpm, tissue_sample_names)
}

# Apply function to all tissue types
invisible(lapply(unique(samples$tissue), get_tissue_samples))

# Save all ordered & mean normalized TE counts per tissue sample
save(counts.norm.rtx, Breast.norm, E_Mucosa.norm, E_Muscularis.norm, Heart.norm, 
     Lung.norm, Prostate.norm, Sk_muscle.norm, Skin.norm,
     file="r_outputs/01-mean_normalized_TE_counts_by_tissue_type.RData")

# Save all ordered &  mean raw TE counts per tissue sample
save(counts.rtx, Breast.raw, E_Mucosa.raw, E_Muscularis.raw, Heart.raw, 
     Lung.raw, Prostate.raw, Sk_muscle.raw, Skin.raw,
     file="r_outputs/01-mean_raw_TE_counts_by_tissue_type.RData")

# Save all ordered &  mean TE CPMs per tissue sample
save(counts.cpm.rtx, counts.cpm.herv, counts.cpm.l1,
     Breast.cpm, E_Mucosa.cpm, E_Muscularis.cpm, Heart.cpm, 
     Lung.cpm, Prostate.cpm, Sk_muscle.cpm, Skin.cpm,
     Breast.herv.cpm, E_Mucosa.herv.cpm, E_Muscularis.herv.cpm, Heart.herv.cpm,
     Lung.herv.cpm, Prostate.herv.cpm, Sk_muscle.herv.cpm, Skin.herv.cpm,
     Breast.l1.cpm, E_Mucosa.l1.cpm, E_Muscularis.l1.cpm, Heart.l1.cpm,
     Lung.l1.cpm, Prostate.l1.cpm, Sk_muscle.l1.cpm, Skin.l1.cpm,
     file="r_outputs/01-mean_cpm_by_tissue_type.RData")

save(counts, counts.rtx, counts.tx, counts.l1, counts.herv,
     counts.cpm, counts.cpm.rtx, counts.cpm.herv, counts.cpm.l1, samples, 
     file="r_outputs/01-counts.Rdata")

remove(Breast.norm, E_Mucosa.norm, E_Muscularis.norm, Heart.norm, 
       Lung.norm, Prostate.norm, Sk_muscle.norm, Skin.norm,
       Breast.raw, E_Mucosa.raw, E_Muscularis.raw, Heart.raw, 
       Lung.raw, Prostate.raw, Sk_muscle.raw, Skin.raw,
       Breast.cpm, E_Mucosa.cpm, E_Muscularis.cpm, Heart.cpm, 
       Lung.cpm, Prostate.cpm, Sk_muscle.cpm, Skin.cpm)

save(retro.annot, herv.annot, file="r_outputs/01-annot.Rdata")

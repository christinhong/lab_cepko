# Christin M. Hong
# Last updated: 2018-09
# Harvard Medical School, Connie Cepko Lab


#### INFRASTRUCTURE ####

# Addressing "max DLL" limit
library(R.utils)
R.utils::gcDLLs(quiet=F)
getLoadedDLLs()

# From https://stackoverflow.com/questions/49740401/increase-max-number-of-dlls-in-r-ubuntu-from-rstudio
# usethis::edit_r_environ() # Currently set R_MAX_NUM_DLLS=256 (default was 100)


#### Import libraries ####

# Data formatting/tidying
library(data.table)
#library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(knitr)
library(httr)
#library(forcats)

# Analysis
library(biomaRt)
library(DESeq2)
library(pheatmap)
library(dendsort)
#library(Rtsne)
#library(cytofkit)
library(sva)
library(limma)

# Visualization
library(Cairo) # For visualizing transparency in ggplots and PDFs
library(viridis)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(gtable)
#library(hexbin)


#### Options ####
options(digits = 5,
        stringsAsFactors = TRUE)
# Some people like setting stringsAsFactors = FALSE, but I find TRUE is more useful.  Still, has enough ramifications to make the option explicit.

# Color palette
colors.div <- brewer.pal(10, "PRGn") # For heatmap


# Seed for reproducibility
set.seed(200)


#### Set variables ####

# Working directory
wd <- "~/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_chickRFZ_rnaSeq/"
setwd(wd)


# Input
vFc <- "doc/featureCounts/featureCounts.txt"
vFcSum <- "doc/featureCounts/featureCounts.txt.summary"


# Output
out.name <- paste0("output/R_")


#### START ####

# # Import data
# fc <- read.table(vFc, header = TRUE)
#
# head(fc)
# tail(fc)
# names(fc)
#     # OK, so columns are Geneid, Chr, Start, End, Strand, Length, and then samples.  Rows are the genes/contigs.
#
#
# # Cleaning up sample names
# names(fc) <- stringr::str_replace(names(fc), "(.*).bam.(.*)", "\\2")
#
#
#
# #### Data annotation ####
#
# # Subset for geneid and counts
# fc2 <- fc[, -c(1:6)]
#
# head(fc2[,1:10])
#
#
# # Set Geneid as rownames
# rownames(fc2) <- fc[,1]
#
#
# # Transpose
# fc.t <- t(fc2)
#
# # Make rownames the first column
# fc.ann <- cbind(rownames(fc.t), data.frame(fc.t, row.names=NULL))
# colnames(fc.ann)[1] <- "Sample"
#
# head(fc.ann[,1:10])
#
#
# # Annotate by sample, batch, and lane
# samples <- fc.ann$Sample
#
# ## Split SeqFile name for annotations
# sample.ann <- stringr::str_split(samples, "_")
#
#
# ## Convert annotations (list of lists) into preferred format
# s.ann2 <- matrix(unlist(sample.ann), nrow=length(unlist(sample.ann[1])))
# head(s.ann2[,1:5])
#
# s.ann3 <- as.data.frame(t(s.ann2))
# head(s.ann3[1:5,])
#
#
# ## Name annotation columns
# colnames(s.ann3)[1] <- "Sample"
# colnames(s.ann3)[2] <- "SeqID"
# colnames(s.ann3)[3] <- "Batch"
# colnames(s.ann3)[4] <- "Lane"
#
#
# ## Repeat for sample annotations
# sample.ann2 <- stringr::str_split(s.ann3[,1], "[.]")
#
# s.ann4 <- matrix(unlist(sample.ann2), nrow=length(unlist(sample.ann2[1])))
# head(s.ann4[,1:5])
#
# s.ann5 <- as.data.frame(t(s.ann4))
# head(s.ann5[1:5,])
#
# colnames(s.ann5)[1] <- "Tissue"
# colnames(s.ann5)[2] <- "Retina"
# colnames(s.ann5)[3] <- "Well"
#
# head(s.ann3)
# head(s.ann5)
#
#
# # Bind annotations and remove uninformative columns
# s.ann6 <- cbind(s.ann5, s.ann3)
#
# head(s.ann6)
#
# s.ann7 <- within(s.ann6, rm("Well", "Sample", "SeqID"))
#
# head(s.ann7)
#
#
#
# # Bind annotations to count data
# fc.ann2 <- cbind(s.ann7, fc.ann)
#
# # Check that everything's in order
# head(fc.ann2[,1:10])
# tail(fc.ann2[,1:6])
#
#
#
# #### Prepping for statistical analyses ####
# ncol(fc.ann2)
#
# Columns of interest
# coi <- c(6:(ncol(fc.ann2))) # Colums with count data
#
# # Transform counts with log2(x+1)
# fc.log <- fc.ann2
# fc.log[, coi] <- log2(fc.ann2[, coi] + 1)
#
# head(fc.log[,1:10])


# # Adjust batch annotation for Retina 7
# ## NEXT TIME: Try something like "fc4$Platform <- ifelse(fc4$Retina <= 5,"NextSeq", "HiSeq")"

# ## Retina 7 is from the same library as Retina 6, but its samples were prepped differently, so I'm going to label it as its own batch.
#
# ## Add new level to Batch factor
# levels(fc.log$Batch)
#
# fc.log2 <- fc.log
#
# levels(fc.log2$Batch) <- c(levels(fc.log2$Batch),"B004")
#
#
# ## Change Batch factor for Retina 7 samples
# fc.log2 %>% dplyr::filter(Retina == 6 | Retina == 7)
#
# ### Extract Retina 7 data and change Batch factor level
# fc.log2.R7 <- fc.log2 %>% dplyr::filter(Retina == 7) %>% dplyr::mutate(Batch = factor(Batch, labels = c("B003" = "B004")))
#
# fc.log2.R7[,1:10]
#
#
# ### Extract data Retina 1-6 and concatenate with updated Retina 7
# fc.log2.Rs <- fc.log2 %>% dplyr::filter(Retina != 7)
#
# fc.log3 <- rbind(fc.log2.Rs, fc.log2.R7)
#
#
# ### Check
# fc.log3[,1:10] %>% dplyr::filter(Retina == 6 | Retina == 7) # Looks good.
#
#
# # Write out annotated and log-transformed count data for faster loading later
# write.table(
#     x = data.frame(fc.log3, stringsAsFactors = T),
#     file = paste0(out.name, "counts-annotated-log2.csv"),
#     quote = F,
#     sep = ",",
#     row.names = F
# )


# #### IMPORT annotated and log2(x+1) transformed data ####
# fc <- read.csv(file = "output/R_2018-09-04_counts-annotated-log2.csv", header = TRUE, stringsAsFactors = TRUE)
#
# head(fc[,1:6])
# tail(fc[, 1:6])
#
#
# # Set rownames
# rownames(fc) <- fc$Sample
#
#
# # Columns of interest
# coi <- c(6:(ncol(fc))) # Colums with count data
#
#
# #### PCA ####
#
# pca.fc <- prcomp(fc[, coi], center = F, scale. = F)
#
# sink(paste0(out.name, "pcaSummary.txt"))
# summary(pca.fc)
# sink()
#
# plot(pca.fc, type = "lines", main = "PCA of RNA-seq by prcomp")
# title(xlab = "Principle Component")
#
# # Would be better if I found a way to tile through the first 6 PCs by factor of interest (for this, Batch, Tissue, Retina, and Lane).
#
# #ggplot(data=meltedPCA) + geom_bar(aes(x=Var1, y=value, fill=group), stat="identity") + facet_wrap(~Var2) # Interesting, but not useful atm.
#
#
# head(pca$x)
#
# scores <- data.frame(fc, pca.fc$x[,1:9])
#
#
# qplot(x=PC1, y=PC2, data=scores) +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     theme(legend.position="right") +
#     theme_light() +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC1, y=PC2, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Retina)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC1, y=PC2, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Tissue)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC2, y=PC3, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC3, y=PC4, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC3, y=PC4, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Tissue)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC4, y=PC5, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC4, y=PC5, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Tissue)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC5, y=PC6, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC5, y=PC6, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Tissue)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC6, y=PC7, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Batch)), size = 4) +
#     ggtitle("(PCA determined by prcomp)") +
#     coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#
#
# qplot(x=PC6, y=PC7, data=scores) +
#     theme(legend.position="right") +
#     theme_light() +
#     geom_point(aes(colour = factor(Tissue)), size = 4) +
#     ggtitle("(PCA determined by prcomp)")
#
#
# # I'm happy to see tissue-specific differences showing up by PC3. Will go ahead and run DESeq2, since that normalizes for library size and can correct for batch as a covariate.
#
#
#
# #### Hierarchical clustering (dendrogram) ####
# # From https://datascienceplus.com/hierarchical-clustering-in-r/
# # For advanced dendrograms: http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
# ## More info: http://dataaspirant.com/2018/01/08/hierarchical-clustering-r/
#
# ## I'm not terribly interested in this right now, but it's recommended to check that the samples cluster as expected by orthogonal methods, and PCA plus hierarchical clustering into a dendrogram is standard.
#
# # Make data.frame with samples as columns
# head(fc[,1:6])
#
# # Cluster
# clusters <- hclust(dist(fc[,coi]))
#
#
# # To view tiny font in full tree, making as a large PDF for zooming in
# pdf(paste0(out.name, "sampleQC-dendrogram.pdf"), width=40)
#
# ## Plot (from https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html)
# plot(clusters, hang = -1, cex = 1)
#
# ## Close the PDF file's associated graphics device (necessary to finalize the output)
# dev.off()
#
#
# # Lanes definitely cluster together here too. Will merge for easier visualization downstream.
#


#### DESeq2 ####
fc <- read.csv(file = "output/R_2018-09-04_counts-annotated-log2.csv", header = TRUE, stringsAsFactors = TRUE)

rownames(fc) <- fc$Sample


#### Format data for DESeq2 input ####
# Following https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf and https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#summarizedexperiment-input


# DESeq2 wants a matrix of raw counts (not log-transformed!) and a table of annotations

## For counts
fc.orig <- read.table(vFc, header = TRUE)

head(fc.orig[,1:7])

# Set rownames to gene
rownames(fc.orig) <- fc.orig$Geneid

# Clean up column names
names(fc.orig) <- stringr::str_replace(names(fc.orig), "(.*).bam.(.*)", "\\2")

# Extract count data
fc.orig.cts <- fc.orig[,7:ncol(fc.orig)]

head(fc.orig.cts)
ncol(fc.orig.cts)

sapply(fc.orig.cts, class) # Integer, good - DESeq2 requires integer



#### Get annotations as a table of factors ####
head(fc[,1:10])

fc.ann <- fc[,1:5]
head(fc.ann)


# Make ID column with tissue and retina
fc.ann$ID <- factor(paste(fc.ann$Tissue, fc.ann$Retina, sep="_"))


## Later on will merge technical replicates (Batch 1 and 2), so going to give them the same identifier now

### If I wanted to correct for batch, would do that before prepping for DESeq2
### NOTE: In the future, will want to more cleanly distinguish between biological batches (samples collected and prepped under different conditions), libraries, and technical replicates.

fc.ann$Batch[fc.ann$Batch == 'B002'] <- 'B001'

### (This causes as "factor levels were dropped which had no samples" warning when running DESeq2 because the B002 level is now empty.)



# Sort by Sample column to get same order of rows as columns in the count matrix
## DESeq2 checks that the rownames of colData match the column names of countData, which is nice.

## Convert to character
fc.ann$Sample <- as.character(fc.ann$Sample)

head(fc.ann)
sapply(fc.ann, class)

## Sort data frame by Sample column
fc.ann2 <- fc.ann[with(fc.ann, order(Sample)), ]
tail(fc.ann2)



# Convert integer columns to factor
sapply(fc.ann2, class)

fc.ann2$Retina <- factor(fc.ann2$Retina)




#### Make DESeq2 object ####
# (note that ncol(countData) must equal nrow(colData))
# dds <- DESeqDataSetFromMatrix(countData = fc.cts2,
#                               colData = fc.ann2,
#                               design= ~ Tissue + Retina + Batch)


## Hmmm...getting 'Model matrix not full rank' error. Happens when variables can be formed by combinations of other variables, e.g. how Retina fits into Batch. Well, I definitely want to adjust for batch and I'm looking for DE between tissues. From the PCA, major differences between retinas are batch-dependent--retinas are theoretically the same. Differences between the platforms (e.g. GC bias) should be captured by batch as well. So...

dds <- DESeqDataSetFromMatrix(countData = fc.orig.cts,
                              colData = fc.ann2,
                              design = ~ Batch + Tissue)

# (dds is for DESeqDataSet)


#### Filter out low count genes ####
# Only keeping genes with greater than or equal to 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]


#### Collapse technical replicates ####
# DESeq2 recommends collapsing technical replicates = multiple sequencing runs of the same library. In this case means batch 1 and batch 2. See non-unique IDs in:
as.data.frame(colData(dds))


# Collapsing based on ID
ddsC <- collapseReplicates(dds,
                           groupby = dds$ID)

    ## Could record value collapsed with "run = dds$Batch", but collapsing on multiple columns here, so won't bother


# Check
as.data.frame(colData(ddsC))


# Check that collapsed total counts match original total counts
origCts <- rowSums( counts(dds)[ , dds$ID == "D_1" ] )

all( origCts == counts(ddsC)[ ,"D_1" ]) # Should output TRUE



#### Set the control/comparison reference ####
# By default, DESeq2 calculates the FIRST level of a factor as the control (a.k.a. reference level, comparison standard, etc.).

# Making RFZ our "control" since we're interested in comparing it to all other tissues
ddsC$Tissue <- relevel( ddsC$Tissue, "RFZ" )

# Check
as.data.frame(colData(ddsC))
levels(as.data.frame(colData(ddsC))$Tissue)



#### Run DESeq2 ####

# DESeq2 estimates size factors to control for library size, dispersion for each gene, and fits a generalized linear model (glm).

deAnalysis <- DESeq(ddsC)


# Check results
resultsNames(deAnalysis)


res <- results(deAnalysis, name="Tissue_D_vs_RFZ") # Is the same as below
res <- results(deAnalysis, contrast = c("Tissue", "D", "RFZ") )

res
mcols(res, use.names = TRUE)


#### p-values and multiple testing ####
# DESeq2 corrects the pvalue by Benjamini-Hochberg to get padj (a type of FDR, more precisely called BJ-adjusted p values)

sum(res$pvalue < 0.01, na.rm = TRUE) # 637

sum(res$padj < 0.1, na.rm = TRUE) # 281, with ~10% as false positives
sum(res$padj < 0.01, na.rm = TRUE) # 85


# Subsetting to genes that pass FDR threshold and sorting by log2 fold change
resSig <- res[ which(res$padj < 0.1 ), ]

resSig.sorted <- resSig[ order( resSig$log2FoldChange ), ]

head(resSig.sorted) # Strongest downregulation from control
tail(resSig.sorted) # Strongest upregulation from control



#### Log fold change shrinkage for visualization and ranking ####
resLFC <- lfcShrink(deAnalysis, coef="Tissue_D_vs_RFZ", type="apeglm")

    # This is a little slow. Note that https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#summarizedexperiment-input has suggestions for parellelizing at this point if it takes more than 30 seconds to run.

# Check
head(resLFC,4)

mcols(resLFC, use.names = TRUE)


# # Write out data (note that this converts rownames to first column)
# write.csv(
#     x = as.data.frame(resLFC, stringsAsFactors = FALSE),
#     file = paste0(out.name, "deseq2-resultsLFC.csv"),
#     quote = F,
#     row.names = TRUE
# )


# resLFC <- read.csv("output/R_2018-09-05_deseq2-resultsLFC.csv")
# head(resLFC)


#### Diagnostic plots ####

# Plot log2 fold change
plotMA( res, ylim = c(-1, 1) )
    # Hmm,  that's kind of a mess.  On the upside, it's symmetrical vertically. On the downside, it looks like the DESeq2 priors aren't very effective. ("The DESeq2 package incorporates a prior on log2 fold changes, resulting in moderated estimates from genes with low counts and highly variable counts, as can be seen by the narrowing of spread of points on the left side of the plot." -https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf) Not much narrowing here.


plotMA( resLFC, ylim = c(-1, 1) )
    # Plotting the shrunken LFC removes noise associated with changes from low count genes without requiring an arbitrary threshold (like FDR). Much cleaner!



# Plot dispersion estimates
## Red line = fit. Black dots = orig data. Blue dots = adjusted data. Black dots circled in blue (non-optimal color choice) = outliers.
plotDispEsts( deAnalysis, ylim = c(1e-6, 1e1) )
    # Looks reasonable - same shape as in https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf


# Histogram of pvalues
hist( resLFC$pvalue, col="grey" )
    # Looks good, peak at low end suggests some differentially expressed genes


# # Look at pvalues vs. normalized count - mainly for justifying filtering to reduce multiple test correction
# ## create bins using the quantile function
# qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# # "cut" the genes into the bins
# bins <- cut( res$baseMean, qs )
# # rename the levels of the bins using the middle point
# levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# # calculate the ratio of £p£ values less than .01 for each bin
# ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# # plot these ratios
# barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")



#### PCA on DESeq2-normalized counts ####
# Get regularized-log transformed counts from DESeq2 (see https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf#subsection.5.1)
rld <- rlog(deAnalysis)
head(assay(rld))


# Annotate
rld.cts <- t(assay(rld))

head(rld.cts[,1:5])

rld.ann <- rownames(rld.cts)

rld.ann2 <- stringr::str_split(rld.ann, "_")


## Convert annotations (list of lists) into preferred format
rld.ann2 <- matrix(unlist(rld.ann2), nrow=length(unlist(rld.ann2[1])))
head(rld.ann2)

rld.ann3 <- as.data.frame(t(rld.ann2))
head(rld.ann3)


## Name annotation columns
colnames(rld.ann3)[1] <- "Tissue"
colnames(rld.ann3)[2] <- "Retina"


## Add to count data
rld2 <- cbind(rld.ann3, rld.cts)
head(rld2[,1:6])


# DESeq's plotPCA
DESeq2::plotPCA(rld, intgroup = "Retina", ntop = 1000, returnData = FALSE)


# PCA
pca.rld <- prcomp(rld2[, 3:ncol(rld2)], center = F, scale. = F)

sink(paste0(out.name, "pcaSummary.txt"))
summary(pca.rld)
sink()

plot(pca.rld, type = "lines", main = "PCA of RNA-seq by prcomp")
title(xlab = "Principle Components")

head(pca.rld$x)

scores <- data.frame(rld2, pca.rld$x[,1:9])


# Plotting
ggplot(data = scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3)

ggplot(data = scores, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3)

ggplot(data = scores, aes(x = PC3, y = PC4)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3)


ggplot(data = scores, aes(x = PC4, y = PC5)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3)


ggplot(data = scores, aes(x = PC5, y = PC6)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3)



# Plotting 2
ggplot(data = scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3)

ggplot(data = scores, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3)

ggplot(data = scores, aes(x = PC3, y = PC4)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3)


ggplot(data = scores, aes(x = PC4, y = PC5)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3)


ggplot(data = scores, aes(x = PC5, y = PC6)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3)



# Plotting top 1000 most variable genes.  Interesting...
gg <- DESeq2::plotPCA( rld, intgroup = c( "Retina"), returnData = TRUE, ntop = 1000)


ggplot(data = gg, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=name),hjust="inward", vjust="inward", size = 3)




# The good news: Retina 6 and Retina 7 cluster together now!

# The bad news is that Retina 1-5 and 6-7 still cluster separately, and that's driving PC1.

# And oddly, when looking at the top 500-1000 variable genes, PC2 is powered by D7. When looking at all genes, D7 drives PC3. Going hunting...

## Reviewing the FastQC report for D7 shows that R2 maps pretty poorly, but that's also true for D6.  OTOH, Qualimap's BAM QC indicates that D7 has an abnormally high 5'-3' bias compared to all other samples (2.36 instead of ~1!). That's probably what's driving it.

## In PC3 vs. PC4, V7 shows up as a strong driver, and it has the second highest 5'-3' bias (1.60).  Next is T7 with 1.47, then N7 1.38, D3 1.37, V5 1.35, and so on.

#  Hmm. I think that from these results, it may be worth testing running svaseq on a tissue-by-tissue basis and using the SVs in the DESeq2 model to hopefully correct for GC bias (which is probably driving the NextSeq vs. HiSeq separation) and the 5'-3' bias.



#### Normalize with sva ####
# See http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# Data matrix for sva
dat  <- counts(deAnalysis, normalized=TRUE)

# Variable table for sva, independent variable of interest is Tissue
mod <- model.matrix(~as.factor(Tissue), colData(deAnalysis))

# Null model with adjustment variables
mod0 <- model.matrix(~ 1, colData(deAnalysis))

# Testing what SV1-12 would look like
svseq.test <- svaseq(dat, mod, mod0, B = 8, n.sv = 12)

svseq.test$sv


## Checking to see if this captures differences in Tissue (hopefully not) and across Retina (hopefully so)
par(mfrow = c(4, 3), mar = c(3,5,3,1))

for (i in 1:12) {
    stripchart(svseq.test$sv[, i] ~ deAnalysis$Tissue, vertical = TRUE, main = paste0("SV", i))
    abline(h = 0)
}


for (i in 1:12) {
    stripchart(svseq.test$sv[, i] ~ deAnalysis$Retina, vertical = TRUE, main = paste0("SV", i))
    abline(h = 0)
}



# Rerunning svaseq for final number of adjusted SVs
svseq <- svaseq(dat, mod, mod0, B = 8, n.sv = 3)

for (i in 1:4) {
    stripchart(svseq$sv[, i] ~ deAnalysis$Retina, vertical = TRUE, main = paste0("SV", i))
    abline(h = 0)
}

dev.off()


# Add SVA data to previous DESeq2 output
ddssva <- deAnalysis

ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~ SV1 + SV2 + SV3 + Tissue



#### Rerun DESeq2 with SVA data ####

# Applying SVA to remove hidden batch effects
deAnalysis.sva <- DESeq(ddssva)

# # 2 genes didn't converge. Per https://support.bioconductor.org/p/65091/ , I tested increasing the iterations, but 1 gene still didn't converge. For simplicity's sake, will keep DESeq and clean out the remaining 2 non-convergent genes
#
# # Genes to be removed
# deAnalysis.sva.orig[!mcols(deAnalysis.sva.orig)$betaConv]
#
#     # rownames(2): ENSGALG00000043263 ENSGALG00000043688
#
#     # ENSGALG00000043263: "Novel gene", no description on Ensembl (https://useast.ensembl.org/Gallus_gallus/Gene/Summary?db=core;g=ENSGALG00000043263;r=6:17094960-17097738;t=ENSGALT00000067505)
#
#     # ENSGALG00000043688: "Novel gene", no description on Ensembl (https://useast.ensembl.org/Gallus_gallus/Gene/Summary?db=core;g=ENSGALG00000043688;r=13:8592131-8600567;t=ENSGALT00000048782)
#
#     # Yeah, I feel fine removing those.
#
#
#
# # Clean input DDS
# deAnalysis.clean <- deAnalysis.sva.orig[which(mcols(deAnalysis.sva.orig)$betaConv),]
#
#
# # Rerun DESeq2
# deAnalysis.sva <- DESeq(deAnalysis.clean)


# Check results for SVs
resultsNames(deAnalysis.sva) # Looks good




#### MA plot ####
# Shrinking log-fold changes by "normal" instead of apeglm because apeglm doesn't work with SVs

resLFC.sva <- lfcShrink(deAnalysis.sva, coef="Tissue_D_vs_RFZ", type="normal")

head(resLFC.sva)

plotMA(resLFC.sva, ylim = c(-1, 1)) # Not working?  May need to make it manually.




#### DESeq's plotPCA ####

pca.orig <- DESeq2::plotPCA(rld.sva,
                intgroup = "Retina",
                ntop = 1000,
                returnData = TRUE)

# Looks exactly the same as without SVs - I mean, *exactly* the same. After spending a couple hours hunting around, I discovered that's how it's supposed to be - standard PCA is for the unadjusted data.  :p  To visualize adjustment, Michael Love (DESeq2 author) suggested limma's removeBatchEffect.



#### Adjust counts by SVA plus limma's removeBatchEffect for visualization ####
# Need to include {limma} removeBatchEffect to see modifications from SVA by PCA

# "For the DE analysis, adding the batch term accounts for mean shifts in the gene counts, and allows one to isolate the effect of a condition of interest, as you know. Typically, we want to look at the PCA plot of the samples in order to see the natural relationship of the samples to each other, so including the batch effects. However, if you want to see the variation among samples, excluding mean shifts due to batch, you can take the matrix of VST or rlog transformed values, and run limma's removeBatchEffect() function on this matrix, then use that matrix as you would the VST or rlog transformed matrix. This function removes mean shifts which can be accounted for using a batch variable." -https://support.bioconductor.org/p/62954/

# and "Just to be clear, there's an important difference between removing a batch effect and modelling a batch effect. Including the batch in your design formula will model the batch effect in the regression step, which means that the raw data are not modified (so the batch effect is not removed), but instead the regression will estimate the size of the batch effect and subtract it out when performing all other tests. In addition, the model's residual degrees of freedom will be reduced appropriately to reflect the fact that some degrees of freedom were "spent" modelling the batch effects. This is the preferred approach for any method that is capable of using it (this includes DESeq2). You would only remove the batch effect (e.g. using limma's removeBatchEffect function) if you were going to do some kind of downstream analysis that can't model the batch effects, such as training a classifier." -https://support.bioconductor.org/p/76099/

rld.sva <- rlog(deAnalysis.sva, blind = FALSE)


# Make SV covariate matrix (necessary when accounting for multiple SVs)
sva.cov <- cbind(deAnalysis.sva$SV1,
                 deAnalysis.sva$SV2,
                 deAnalysis.sva$SV3)


# Adjust for SV effects on counts with {limma}
rld.cts.sva <- limma::removeBatchEffect(assay(rld.sva),
                                        covariates = sva.cov)

# Transpose and annotate
rld.cts.sva <- t(rld.cts.sva)

rld.ann.sva <- rownames(rld.cts.sva)

rld.ann2.sva <- stringr::str_split(rld.ann.sva, "_")


## Convert annotations (list of lists) into preferred format
rld.ann2.sva <- matrix(unlist(rld.ann2.sva), nrow=length(unlist(rld.ann2.sva[1])))
head(rld.ann2.sva)

rld.ann3.sva <- as.data.frame(t(rld.ann2.sva))
head(rld.ann3.sva)


## Name annotation columns
colnames(rld.ann3.sva)[1] <- "Tissue"
colnames(rld.ann3.sva)[2] <- "Retina"


## Add to count data
rld2.sva <- cbind(rld.ann3.sva, rld.cts.sva)
head(rld2.sva[,1:6])



#### Heatmap ####
# See http://edu.sablab.net/rt2018/scripts4.html for some beautiful heatmaps

# scripts from https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix (will improve on this later)
# # transform to log2(n+1)
# de.log <- normTransform(deAnalysis.sva)
#
# select <- order(rowMeans(counts(deAnalysis.sva, normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
#
# df <- as.data.frame(colData(deAnalysis.sva)[,c("Retina","Tissue")])
#
# # plot
# pheatmap(assay(de.log)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)


# For distances plus map

sampleDists <- dist(rld2.sva[3:ncol(rld2.sva)])


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld2.sva$Retina, rld2.sva$Tissue, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(out.name, "deseq2-heatmap-3svs_", Sys.Date(), ".pdf"))

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dev.off()


# That looks better.  See some clustering by tissue type



#### Heatmap for most variable genes ####
head(rld2.sva[,1:10])

# Subset for count data
rld2.sva.cts <- rld2.sva[, 3:ncol(rld2.sva)]

# Transpose
rld2.sva.cts.t <- t(rld2.sva.cts)

head(rld2.sva.cts.t[,1:10])
class(rld2.sva.cts.t) # Checking for numeric matrix


# Estimate variance for each row
var_genes <- apply(rld2.sva.cts.t, 1, var)

head(var_genes)


# Get gene names for the top most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]

head(select_var)


# Subset logcounts matrix for most variable genes
rld2.sva.var <- rld2.sva.cts.t[select_var,]

dim(rld2.sva.var)

head(rld2.sva.var)


# Transpose count matrix again
rld2.sva.var.t <- t(rld2.sva.var)
head(rld2.sva.var.t)


# Heatmap
sampleDists <- dist(rld2.sva.var.t)

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rld2.sva$Retina, rld2.sva$Tissue, sep="-")

colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(out.name, "deseq2-heatmap-3svs-1000topVar_", Sys.Date(), ".pdf"))

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dev.off()




#### PCA on SVA-adjusted counts ####

pca.rld.sva <- prcomp(rld2.sva[, 3:ncol(rld2.sva)], center = F, scale. = F)

sink(paste0(out.name, "pcaSummary-sva-3svs_",  Sys.Date(), ".txt"))
summary(pca.rld.sva)
sink()


plot(pca.rld.sva, type = "lines", main = "PCA of RNA-seq by prcomp")
title(xlab = "Principle Components")

head(pca.rld.sva$x)

scores.sva <- data.frame(rld2.sva, pca.rld.sva$x[,1:9])


# Get percent variance explained by each PC
head(pca.rld.sva$sdev)

pcVar <- round((((pca.rld.sva$sdev)^2) * 100) / sum(pca.rld.sva$sdev^2), digits = 2)

cumsum((pca.rld.sva$sdev)^2) / sum(pca.rld.sva$sdev^2)

#### PDF ####
pdf(paste0(out.name, "deseq2-pca-3svs_", Sys.Date(), ".pdf"))

# Plotting
ggplot(data = scores.sva, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp, \nPC1 scaled by 0.2 for easier viewing)") +
    coord_fixed(ratio = 0.2, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC1: ",pcVar[1],"% variance")) +
    ylab(paste0("PC2: ",pcVar[2],"% variance"))


ggplot(data = scores.sva, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp, \nPC1 scaled by 0.2 for easier viewing)") +
    coord_fixed(ratio = 0.2, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC1: ",pcVar[1],"% variance")) +
    ylab(paste0("PC2: ",pcVar[2],"% variance"))


dev.off()

# Looks promising, but adjusting for SV1 results in PC1 accounting 99.13% of data. Well, considering the batch/retina confounders for this data, that's probably valid...

## Iterating through adjusting for 1-5 SVs. By SV4, retinas look well combined (with SV1-3, retina 6 still clusters together). Will stay with adjusting for SV1-4. (PC1 for SV1-4 99.44% of variance. With SV1-3, PC1 is 99.3%, so...yeah.)

## That said, will definitely go back and rerun DE analysis on only Retinas 1-5 to compare.

ggplot(data = scores.sva, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC2: ",pcVar[2],"% variance")) +
    ylab(paste0("PC3: ",pcVar[3],"% variance"))


ggplot(data = scores.sva, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC2: ",pcVar[2],"% variance")) +
    ylab(paste0("PC3: ",pcVar[3],"% variance"))


ggplot(data = scores.sva, aes(x = PC3, y = PC4)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC3: ",pcVar[3],"% variance")) +
    ylab(paste0("PC4: ",pcVar[4],"% variance"))


ggplot(data = scores.sva, aes(x = PC3, y = PC4)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp)") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC3: ",pcVar[3],"% variance")) +
    ylab(paste0("PC4: ",pcVar[4],"% variance"))




#### PCA on top most variable genes ####

pca.rld.sva.var <- prcomp(rld2.sva.var.t, center = F, scale. = F)

sink(paste0(out.name, "pcaSummary-sva-3svs-1000topVar_",  Sys.Date(), ".txt"))
summary(pca.rld.sva.var)
sink()


plot(pca.rld.sva.var, type = "lines", main = "PCA of RNA-seq by prcomp")
title(xlab = "Principle Components")

head(pca.rld.sva.var$x)

scores.sva.var <- data.frame(rld2.sva.var.t, pca.rld.sva.var$x[,1:9], rld2.sva[,1:2])


# Get percent variance explained by each PC
head(pca.rld.sva.var$sdev)

pcVar <- round((((pca.rld.sva.var$sdev)^2) * 100) / sum(pca.rld.sva.var$sdev^2), digits = 2)

cumsum((pca.rld.sva.var$sdev)^2) / sum(pca.rld.sva.var$sdev^2)


#### PDF ####
pdf(paste0(out.name, "deseq2-pca-3svs-1000topVar_", Sys.Date(), ".pdf"))

# Plotting

# Orig
percentVar <- round(100 * attr(pca.orig, "percentVar"), digits = 3)

ggplot(data = pca.orig, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(Original DESeq2 PCA determined by prcomp, \nPC1 scaled by 0.2 for easier viewing)") +
    coord_fixed(ratio = 0.2, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=name),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))


# SVA
ggplot(data = scores.sva.var, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Retina)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp, \nPC1 scaled by 0.2 for easier viewing)") +
    coord_fixed(ratio = 0.2, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Tissue),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC1: ",pcVar[1],"% variance")) +
    ylab(paste0("PC2: ",pcVar[2],"% variance"))


ggplot(data = scores.sva.var, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = factor(Tissue)), size = 4) +
    theme(legend.position="right") +
    theme_light() +
    ggtitle("(PCA determined by prcomp, \nPC1 scaled by 0.2 for easier viewing)") +
    coord_fixed(ratio = 0.2, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    geom_text(aes(label=Retina),hjust="inward", vjust="inward", size = 3) +
    xlab(paste0("PC1: ",pcVar[1],"% variance")) +
    ylab(paste0("PC2: ",pcVar[2],"% variance"))


dev.off()



#### Adding gene names ####
# Requires a working internet connection

resLFC$ensembl <- sapply( strsplit( rownames(resLFC), split="\\+" ), "[", 1 )
# DESeq2 assumes the gene names are ensembl, thus this code for generating the "res$ensembl" filter for getBM.

ensembl <- biomaRt::useMart( "ensembl", dataset = "ggallus_gene_ensembl") # Galgal5

genemap <- biomaRt::getBM( attributes = c("ensembl_gene_id",
                                          "entrezgene",
                                          "hgnc_symbol",
                                          "description",
                                          "chromosome_name",
                                          "band",
                                          "strand",
                                          "start_position",
                                          "end_position"),
                           filters = "ensembl_gene_id",
                           values = resLFC$ensembl,
                           mart = ensembl )



# # Write out genemap data faster loading later
# write.table(
#     x = data.frame(genemap, stringsAsFactors = F),
#     file = paste0(out.name, "biomaRt-ensembl-genemap.csv"),
#     quote = F,
#     sep = ",",
#     row.names = F
# )

idx <- match(resLFC$ensembl, genemap$ensembl_gene_id)

resLFC$entrez <- genemap$entrezgene[ idx ]

resLFC$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

resLFC$hgnc_description <- genemap$description[ idx ]




#### p-values and adjusted p-values ####
sum(res.sva$pvalue < 0.01, na.rm = TRUE) # 444

sum(res.sva$padj < 0.1, na.rm = TRUE) # 192, with ~10% as false positives
sum(res.sva$padj < 0.01, na.rm = TRUE) # 66


# Subsetting to genes that pass FDR threshold and sorting by log2 fold change
resSig.sva <- res.sva[ which(res.sva$padj < 0.1 ), ]

resSig.sorted.sva <- resSig.sva[ order( resSig.sva$log2FoldChange ), ]

head(resSig.sorted.sva) # Strongest downregulation from control
tail(resSig.sorted.sva) # Strongest upregulation from control





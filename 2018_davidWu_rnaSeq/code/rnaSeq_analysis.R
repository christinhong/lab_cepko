#### Background ####

# Author: Christin M. Hong
# Objective: Analyzing processed microarray data
# Team: Connie Cepko, Harvard Medical School
# Last updated: 2018-06

# This data was provided as counts in Excel spreadsheets.  I'm assuming it's already been normalized to FPKMs.


#### START INFRASTRUCTURE ####
#### Import libraries ####
# Standard data tidying
library(tidyr)
library(dplyr)
library(stringr)
library(lubridate) # For standardizing and manipulating dates and times
library(httr) # working with urls

# Standard visualization
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(gtable)
library(ggsignif)
library(Cairo) # For visualizing transparency in ggplots and PDFs
library(viridis)
library(RColorBrewer)
library(ezknitr)

# Data analysis
library(pheatmap)
library(dendsort)
# library(cytofkit) # for phenograph, tSNE, etc.  May require R 3.5 to deal with "max DLL" issue
# library(sva) # for ComBat for batch correction


#### Color palettes ####
colors.pg <- brewer.pal(10, "PRGn") # For heatmap


#### Set seed for reproducibility ####
set.seed(200)


#### Set variables ####
# What's the working directory?  Preferably in same location as the data and Rproject.
getwd()

fname <- "dw_microarray_rpe-innateImmunity"

#### Import data ####
d <- read.csv("rnaSeq/output/innate immune genes tissue microarraysgrant ed dw copy.csv")
# I ended up manually editing the sample names in this file because it was the easiest way to generate usable annotations...thankfully there were only 10 samples.  I'm sure there are more elegant ways of doing it, though (what if there were 100s of samples?!).


#### * File name * ####
out.name <-
    paste0("output/",
           fname,
           "_",
           "rnaSeq-analysis_",
           Sys.Date())

print(out.name)


#### END INFRASTRUCTURE ####


#### START ANALYSIS ####

# look at number/distribution of reads, including only complete cases
sapply(d, class)

d <- transform(d, RNAseq = as.integer(RNAseq))

hist(d$RNAseq) # This is odd.  Are these p-values or read counts?

colnames(d)


# Graph log2(x+1) RNAseq values (log = natural log)
hist(t(log2(d[4:11] + 1)))


#### GUESSES ####
head(d)

# Are these raw reads...?  I hope not.  I guess I'll assume that these have at least been scaled to library size, since I don't have the full library to scale with anyway.  (Not sure if it's different for microarray data, but I'm still sure I'm only working with a subset of the data, so assuming that it comes scaled is still reasonable.)

# Will also assume data has already been filtered.

# It looks like there are 2 mouse strains, 3 timepoints, and 2 types of tissue (whole retina vs. RPE?).  There are replicates at the timepoints for the unlabeled tissue type, which I'm going to assume is whole retina.


#### Considering multiple parameters ####
# #### Annotation ####
# # transpose data to annotate by sample (converts all to "character" class because this include numbers and text.  But the character class makes things tricky later on...hmm.)
#### IDEA: Maybe annotate horizontally.  Kind of weird, but should work. ####
# d.t <- t(d)
# head(d.t)
#
# # Exclude n and Unigene row
# d.t2 <- tail(d.t, -2)
# nrow(d.t2)
#
# # rearrange rows
# d.t3 <- d.t2[c(10, 1:9), ]
# head(d.t3)
#
# # annotating
# samples <- stringr::str_split(rownames(d.t3), "_", simplify = TRUE)
#
# samples.df <- as.data.frame(samples)
# samples.df2 <- samples.df
# names(samples.df2)[1] <- "MouseStrain"
# names(samples.df2)[2] <- "Tissue"
# names(samples.df2)[3] <- "Timepoint"
# names(samples.df2)[4] <- "Replicate"
# print(samples.df2)
#
#
# d.full <- cbind(samples.df2, d.t3)


# which ones are most useful for testing differential expression...?
head(d)
colnames(d)
sapply(d,  class)

#### Editing a mislabeled values ####
# duplicate tnfaip8 in gene.name?
dplyr::filter(d, gene.name == "tnfaip8")
# Oh, I see, it's mislabeled in one row---should be tnfaip9 in row with Unigene Mm.31403.  Also:

dplyr::filter(d, grepl("26103", gene.name)) # filtering for partial match with grep

# correcting
d$gene.name <- as.character(d$gene.name)

# subset for correct values (AND to keep values that are correct for all counts)
d.ok <- dplyr::filter(d, Unigene != "Mm.31403" & Unigene != "Mm.45995")

# subset for incorrect (OR to get all errors)
d.error <- dplyr::filter(d, Unigene == "Mm.31403" | Unigene == "Mm.45995")

d.error$gene.name[1] <- "sting"
d.error$gene.name[2] <- "tnfaip9"
d.error

# merge
d.corrected <- rbind(d.ok, d.error)

# check
dplyr::filter(d.corrected, gene.name == "tnfaip8")
tail(d.corrected)


#### tidying ####
# Renaming rows by gene name
d2 <- d.corrected
rownames(d2) <- d.corrected[, 12]

# log2(x+1) scaling
d.val <- log2(d2[4:11] + 1)


#### PCA and k means ####
autoplot(
    kmeans(
        d.val,
        centers = 5,
        iter.max = 100,
        nstart = 20,
        algorithm = "Lloyd"
    ),
    main = paste0("k means clustering"),
    data = d,
    asp = 1
)

flow_pcaPr <- prcomp(d[, c(4:8)], center = F, scale. = F)
summary(flow_pcaPr)

plot(flow_pcaPr, type = "lines", main = "PCA by prcomp")
title(xlab = "Principle Component")



#### Heatmap ####
# From Kam's tutorial: http://slowkow.com/notes/heatmap-tutorial/

## Setting coloring scheme based on quantile breaks so coloring will represent equal proportion of data
quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0 , 1, length.out = n))
    breaks[!duplicated(breaks)]
}

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(d.val, n = 11))


## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(d.val)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(d.val)))


## Plotting sorted heatmap
pheatmap(
    mat = d.val,
    color = inferno(10),
    border_color = NA,
    show_colnames = T,
    show_rownames = T,
    drop_levels = T,
    fontsize = 10,
    main = paste0("Sample vs. log2(x+1) read count with sorted dendrogram\nFrom file: ", fname),
    kmeans_k = NA,
    breaks = flowHM_breaks,
    cluster_cols = flowHM_cluster_cols,
    cluster_rows = flowHM_cluster_rows
)
# OK, sorting makes sense, that's cool.


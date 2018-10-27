# Christin M. Hong
# Last updated: 2018-10
# Harvard Medical School, Connie Cepko Lab

# Formatting filenames from Broad Genomics Platform output from Picard to a human-readable annotation syntax.


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
library(genefilter)

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
        stringsAsFactors = FALSE)

# Color palette
colors.div <- brewer.pal(10, "PRGn") # For heatmap


# Seed for reproducibility
set.seed(200)


#### Set variables ####

# Working directory
wd <- "~/Dropbox/01_Harvard/02_code_github/lab_cepko/2018_seanWang_microglia/"
setwd(wd)



#### START ####

# Goal: "Sample-Well_S##_L###_R1orR2_batchNumber.fastq.gz"


#### Creating annotation table ####
# Import original sample filenames
s.orig <- read.table("sampleFiles.txt")
head(s.orig)


# Split original filenames on periods
s.orig.split <- stringr::str_split(s.orig$V1, "[.]")


# Get lane, barcode, and read data from split names
s.ann1 <- matrix(unlist(s.orig.split), nrow=length(unlist(s.orig.split[1])))

# Check
s.ann1[1:7,1:5]


# Transpose and convert to dataframe
s.ann2 <- as.data.frame(t(s.ann1))
head(s.ann2[1:5,])


# Name annotation columns
colnames(s.ann2)[1] <- "Machine"
colnames(s.ann2)[2] <- "Lane"
colnames(s.ann2)[3] <- "Barcode"
colnames(s.ann2)[4] <- "unmapped"
colnames(s.ann2)[5] <- "Read"
colnames(s.ann2)[6] <- "fastq"
colnames(s.ann2)[7] <- "gz"


# Append original list of files
length(s.orig$V1) # 768
nrow(s.ann2) # 768

s.ann3 <- cbind(s.orig$V1, s.ann2)


# Check that barcodes in original name and Barcode column match
s.ann3[1:10,] # Looks good


# Import data annotations
anno <- read.csv("SK-3JTV_metadata.csv", header = TRUE)
head(anno)
nrow(anno) # 96


# Merge s.ann2 and anno by "Barcode" column
anno2 <- dplyr::left_join(s.ann3, anno, by = "Barcode")
    # Using left join to retain info from original file names
    # "If there are multiple matches between x and y, all combinations of the matches are returned." ---https://dplyr.tidyverse.org/reference/join.html

head(s.ann3)
head(anno)
nrow(anno2) # 768
head(anno2)
tail(anno2)

sapply(anno2, class) # Should all be "character" since I set the default option stringsAsFactors=FALSE



#### Formatting name components ####

# Goal: "Sample-Well_S##_L###_R1orR2_batchNumber.fastq.gz"

# Sort by well
anno3 <- anno2[order(anno2$Well), ]


# Add sample number based on well position (0 padded)
anno3$Well <- as.factor(anno3$Well)

sapply(anno3, class)

anno3$sNum <- as.numeric(anno3$Well)

anno3 <- transform(anno3, sNum = sprintf('S%03d', sNum))


# Format lanes
anno3 <- transform(anno3, Lane = sprintf('L00%s', Lane))


# Format reads
anno3 <- transform(anno3, Read = sprintf('R%s', Read))


# Add batch number and suffix
anno3$suffix <- "001.fastq.gz"


# Append Well to Sample
anno3 <- within(anno3,  SampleW <- paste(Sample, Well, sep="-"))


# Check
anno3[1:12,]
tail(anno3)



### Create new file names ####

# Goal: "Sample-Well_S##_L###_R1orR2_batchNumber.fastq.gz"

anno.final <- within(anno3, newName <- paste(SampleW, sNum, Lane, Read, suffix, sep="_"))

anno.final[1:12,]
tail(anno.final) # Yes!


# Rename original filename column to be bash-friendly
anno.final <- dplyr::rename(anno.final, oldName = s.orig.V1)


### Write out results ####
write.csv(
    x = data.frame(anno.final),
    file = "sampleAnno.csv",
    quote = F,
    row.names = F
)



#### Write out just names ####
anno.names <- anno.final %>% dplyr::select(oldName, newName)

head(anno.names)

write.csv(
    x = data.frame(anno.names),
    file = "sampleNames.csv",
    quote = F,
    row.names = F
)



# Done!

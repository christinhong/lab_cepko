#### Intro ####

# Author: Christin M. Hong
# Objective: Tidying qPCR plate layouts and cycle values
# Last updated: 2018-05
# Team: Cepko, Harvard Medical School


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

# Data analysis
library(XLConnect) # for working with 2007 Excel sheets in R
library(readxl) # for working with newer Excel sheets in R


#### Color palettes ####
colors.pg <- brewer.pal(10, "PRGn") # For heatmap


#### Set seed for reproducibility ####
set.seed(200)


#### Set variables ####
# What's the working directory?  Preferably in same location as the data and Rproject.
getwd()

#### Plate layout ####
fname <- "20180530_2h_plate1.xlsx"

layout <- readxl::read_excel(paste0("layouts/", fname))


#### Results ####
# Load in Workbook with XLConnect to account for older version of Excel
wb <- XLConnect::loadWorkbook("results/20180530 2hours-1 -  Quantification Summary.xlsx")

# Load in worksheet
results <- XLConnect::readWorksheet(wb, sheet=1)


#### * File name * ####
out.name <-
    paste0("output/",
           fname,
           "_",
           "qPCR-processing_",
           Sys.Date())

print(out.name)


#### END INFRASTRUCTURE ####


#### Tidy data ####
# Extract layout
pl1 <- layout[c(2:9), c(2:13)]
print(pl1)

# Convert format to A1-A12, B1-B12, etc. (same as qPCR machine output)
pl1.v <- as.data.frame(c(t(pl1)))

# Split names to label samples
samples <- str_split(pl1.v$`c(t(pl1))`, "_", simplify = TRUE)
samples.df <- as.data.frame(samples)
print(samples.df)

## naming labels
samples.df2 <- samples.df
names(samples.df2)[1] <- "Primer"
names(samples.df2)[2] <- "Group"
names(samples.df2)[3] <- "Mouse"
print(samples.df2)

# Add Ct
plate.df <- cbind(samples.df2, results$Cq)
print(plate.df)

# Export as CSV
write.table(
    x = data.frame(plate.df, stringsAsFactors = T),
    file = paste0(out.name, ".csv"),
    quote = F,
    sep = ",",
    row.names = F
)



#### End tidying ####

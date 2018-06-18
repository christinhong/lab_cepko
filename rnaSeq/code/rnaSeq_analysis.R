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

# Data analysis
library(XLConnect) # for working with 2007 Excel sheets in R
library(readxl) # for working with newer Excel sheets in R
library(edgeR)


#### Color palettes ####
colors.pg <- brewer.pal(10, "PRGn") # For heatmap


#### Set seed for reproducibility ####
set.seed(200)


#### Set variables ####
# What's the working directory?  Preferably in same location as the data and Rproject.
getwd()

fname <- "dw_microarray_rpe-innateImmunity"

#### Import data ####
layout <- readxl::read_excel("output/innate immune genes tissue microarraysgrant ed dw copy.xlsx")


#### * File name * ####
out.name <-
    paste0("output/",
           fname,
           "_",
           "rnaSeq-analysis_",
           Sys.Date())

print(out.name)


#### END INFRASTRUCTURE ####


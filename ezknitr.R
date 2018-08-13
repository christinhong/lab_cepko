library(ezknitr) # I can play with this more.  See https://deanattali.com/blog/ezknitr-package/

setwd("/home/christin/aaa_present_CepkoLab/lab_cepko-code/")

ezknit(file = "rnaSeq/code/rnaSeq_analysis.Rmd",
       out_dir = "rnaSeq/reports",
       fig_dir = "rnaSeq/output/figures",
       verbose = TRUE,
       keep_html = FALSE)


#!/usr/bin/env Rscript

library(funr)

#args = commandArgs(trailingOnly = TRUE)
wd <- getwd()
#print(wd)
#print(dirname(sys.script()))
from_dir <- dirname(sys.script())
to_dir <- wd

# create cazzlink of Rmd file
file.copy(paste0(from_dir, "/report.Rmd"), to_dir, overwrite = TRUE)

rmarkdown::render("report.Rmd", output_format = 'html_document', output_file = 'MToolBox_report.html')

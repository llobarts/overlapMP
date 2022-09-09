#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
name = args[1]

currentpath = getwd()
pathnewfile = paste0(currentpath, "/Results/Reports/Report_overlapping_", name, '.html')
pathrmarkdown = paste0(currentpath, "/Models/psiblast_report.Rmd")
rmarkdown::render(input = pathrmarkdown,
                  output_format = "html_document",
                  output_file = pathnewfile,
                  envir = new.env(),
                  params = list(
                    name = name
                    ))



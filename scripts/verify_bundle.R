#!/usr/bin/env Rscript
# Detect accidental changes to the frozen reference data.
source("R/workflow.R")
m <- csv_read("reference/MD5SUMS.csv")
actual <- unname(tools::md5sum(m$file))
bad <- is.na(actual) | actual != m$md5
if(any(bad)) stop("Missing or modified reference files: ",paste(m$file[bad],collapse=", "))
message("All ",nrow(m)," frozen reference files match their checksums.")

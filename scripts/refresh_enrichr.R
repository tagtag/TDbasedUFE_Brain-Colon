#!/usr/bin/env Rscript
# Optional live analysis. Never overwrites manuscript snapshots.
# Run from the repository root after full analysis, or supply a symbols CSV.
source("R/workflow.R")
need("curl"); need("jsonlite")
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)) args[1] else "results/full/selected_gene_symbols.csv"
d <- csv_read(input)
symbols <- unique(d$gene_symbol[!is.na(d$gene_symbol) & nzchar(d$gene_symbol)])
out <- "results/enrichr_live"
initialize(out)
h <- curl::new_handle()
curl::handle_setform(h,list=paste(symbols,collapse="\n"),
  description="Cognitive aging tensor-selected mouse gene symbols")
response <- curl::curl_fetch_memory("https://maayanlab.cloud/Enrichr/addList",handle=h)
if(response$status_code!=200L) stop("Enrichr submission failed; no retry")
submission <- jsonlite::fromJSON(rawToChar(response$content))
jsonlite::write_json(submission,file.path(out,"submission.json"),auto_unbox=TRUE,pretty=TRUE)
libs <- sub("\\.tsv$","",list.files("reference/enrichr",pattern="\\.tsv$"))
for(lib in libs) {
  message("Enrichr: ",lib)
  url <- paste0("https://maayanlab.cloud/Enrichr/export?userListId=",submission$userListId,
    "&filename=",lib,"&backgroundType=",utils::URLencode(lib,reserved=TRUE))
  download_once(url,file.path(out,paste0(lib,".tsv")))
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

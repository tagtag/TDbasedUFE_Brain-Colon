#!/usr/bin/env Rscript
# Run from the repository root. No automatic installation or retries.
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else "reference"
if(!mode %in% c("reference","download","full"))
  stop("Usage: Rscript run_analysis.R [reference|download|full]")
source("R/workflow.R")
out <- file.path("results",mode)
initialize(out)
if(mode=="download") {
  download_inputs()
} else {
  source("R/correlations.R")
  source("R/enrichment.R")
  source("R/validation.R")
  if(mode=="full") {
    need("rTensor"); need("readxl")
    download_inputs()
    matched <- prepare_samples(out)
    source("R/tensor_analysis.R")
    analysis <- run_tensor_analysis(matched,out)
  } else {
    need("readxl")
    analysis <- list(scores=csv_read("reference/molecular_scores.csv"),
      selected=csv_read("reference/selected_gene_symbols.csv"),
      background=csv_read("reference/background_gene_ids.csv")[[1]])
  }
  compute_correlations(analysis$scores,out)
  ens_enrichment(analysis$selected$gene_id,analysis$background,out)
  comparison <- compare_boruta(analysis$selected,out)
  boruta_excluded_enrichment(analysis$selected,comparison,out)
  enrichr_tables(out)
  validate_results(out,full=mode=="full")
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
message("Completed mode: ",mode,". Outputs: ",out)

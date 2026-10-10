validate_results <- function(out, full=FALSE) {
  checks <- list()
  check <- function(name,ok,detail="") {
    checks[[length(checks)+1L]] <<- data.frame(check=name,passed=isTRUE(ok),detail=detail)
  }
  numerical_match <- function(actual,reference,keys,columns,tolerance=1e-8) {
    key <- function(x) do.call(paste,c(x[keys],sep="::"))
    idx <- match(key(reference),key(actual))
    if(anyNA(idx) || nrow(actual)!=nrow(reference)) return(FALSE)
    all(vapply(columns,function(col) isTRUE(all.equal(as.numeric(actual[[col]][idx]),
      as.numeric(reference[[col]]),tolerance=tolerance,check.attributes=FALSE)),logical(1)))
  }
  corr <- csv_read(file.path(out,"correlations_all_36.csv"))
  ref <- csv_read("reference/correlations.csv")
  check("All 36 Pearson tests", numerical_match(corr,ref,c("molecular_score","behavior"),
    c("n","r","p_value","BH_adjusted_p")))
  adj <- csv_read(file.path(out,"age_adjusted_correlations_all_36.csv"))
  ref <- csv_read("reference/age_adjusted_correlations.csv")
  check("All 36 categorical-age-adjusted tests",numerical_match(adj,ref,
    c("molecular_score","behavior"),c("n","df","partial_r","p_value","BH_adjusted_p")))
  ens <- csv_read(file.path(out,"ENS_all_32_tests.csv"))
  ref <- csv_read("reference/ens/fisher_all_32_tests.csv")
  check("All 32 ENS marker tests",numerical_match(ens,ref,"marker_set",
    c("overlap_n","marker_genes_in_background","fold_enrichment","p_value","BH_adjusted_p")))
  boruta <- csv_read(file.path(out,"Boruta_union_summary.csv"))
  check("Boruta union 506 / overlap 64 / TD-only 1513",
    boruta$Boruta_genes==506 && boruta$overlap==64 && boruta$TD_only==1513)
  comp <- csv_read(file.path(out,"Boruta_overlap_summary.csv"))
  check("All Boruta subsets",numerical_match(comp,csv_read("reference/boruta/overlap_summary.csv"),
    c("source","minimum_runs"),c("features","unique_genes","overlap_with_TD","TD_not_selected")))
  enrich <- csv_read(file.path(out,"TD_only_enrichment_all_terms.csv"))
  ref <- csv_read("reference/boruta/TD_only_local_enrichment.csv")
  check("All Boruta-excluded enrichment terms",numerical_match(enrich,ref,c("library","term"),
    c("universe_genes","query_genes","term_genes","overlap","p_value","BH_adjusted_p")))
  if(full) {
    scores <- csv_read(file.path(out,"molecular_scores.csv"))
    check("Rebuilt molecular scores",numerical_match(scores,csv_read("reference/molecular_scores.csv"),
      "sample_id",names(layer_files),tolerance=1e-6))
    check("75 unique mice",nrow(scores)==75 && !anyDuplicated(scores$sample_id))
  }
  result <- do.call(rbind,checks)
  csv_write(result,file.path(out,"validation.csv"))
  print(result,row.names=FALSE)
  if(!all(result$passed)) stop("Reproduction check failed; inspect validation.csv and outputs.")
  invisible(result)
}

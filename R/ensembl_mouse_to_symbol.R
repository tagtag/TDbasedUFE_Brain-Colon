# Base R only. Reuse the GTF downloaded by test_rna_methylation_overlap.R.
# No downloads, package installation, parallel processing, or retries.
# source("ensembl_mouse_to_symbol.R")
# symbols <- ensembl_mouse_to_symbol(S_gene$selected_ids)

ensembl_mouse_to_symbol <- function(
    gene_ids,
    gtf_file="data/Mus_musculus.GRCm39.116.gtf.gz",
    output_file="selected_gene_symbols.csv") {
  original <- as.character(gene_ids)
  ids <- sub("\\.[0-9]+$","",trimws(original))
  if (anyNA(ids) || any(!grepl("^ENSMUSG[0-9]+$",ids)))
    stop("gene_ids must contain mouse Ensembl gene IDs (ENSMUSG).")
  if (!file.exists(gtf_file)) stop("GTF not found: ",gtf_file)
  con <- if (grepl("\\.gz$",gtf_file)) gzfile(gtf_file,"rt") else file(gtf_file,"rt")
  on.exit(close(con),add=TRUE)
  pieces <- list()
  repeat {
    lines <- readLines(con,n=50000L,warn=FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines,"#") &
      grepl("^[^\t]+\t[^\t]+\tgene\t",lines)]
    if (!length(lines)) next
    attributes <- sub("^(?:[^\t]*\t){8}","",lines,perl=TRUE)
    extract <- function(key) {
      pattern <- paste0('(?:^|;)[[:space:]]*',key,' "([^\"]+)"')
      matches <- regmatches(attributes,regexec(pattern,attributes,perl=TRUE))
      vapply(matches,function(x) if (length(x)>=2L) x[2L] else NA_character_,character(1L))
    }
    genes <- sub("\\.[0-9]+$","",extract("gene_id"))
    symbols <- extract("gene_name")
    keep <- !is.na(genes) & genes %in% ids & !is.na(symbols) & nzchar(symbols)
    if (any(keep)) pieces[[length(pieces)+1L]] <- data.frame(
      gene_id=genes[keep],gene_symbol=symbols[keep],stringsAsFactors=FALSE)
  }
  mapping <- if (length(pieces)) unique(do.call(rbind,pieces)) else
    data.frame(gene_id=character(),gene_symbol=character())
  count <- table(mapping$gene_id)
  ambiguous <- names(count)[count>1L]
  unique_mapping <- mapping[!mapping$gene_id %in% ambiguous,,drop=FALSE]
  symbol <- unique_mapping$gene_symbol[match(ids,unique_mapping$gene_id)]
  result <- data.frame(original_id=original,gene_id=ids,gene_symbol=symbol,
    status=ifelse(ids %in% ambiguous,"ambiguous",
                  ifelse(is.na(symbol),"unmapped","mapped")),stringsAsFactors=FALSE)
  # Preserve input order, repeated input IDs and repeated symbols across genes.
  if (!is.null(output_file)) {
    dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
    write.csv(result,output_file,row.names=FALSE,na="")
  }
  message("Mapped ",sum(result$status=="mapped")," / ",nrow(result)," input IDs.")
  result
}

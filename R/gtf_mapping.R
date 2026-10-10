normalize_mouse_ensembl <- function(ids) sub("\\.[0-9]+$", "", trimws(as.character(ids)))

read_mouse_gtf_mapping <- function(gtf_file) {
  con <- if (grepl("\\.gz$",gtf_file)) gzfile(gtf_file,"rt") else file(gtf_file,"rt")
  on.exit(close(con),add=TRUE)
  pieces <- list()
  repeat {
    lines <- readLines(con,n=50000L,warn=FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines,"#") &
                   grepl("^[^\t]+\t[^\t]+\ttranscript\t",lines)]
    if (!length(lines)) next
    has_ids <- grepl('gene_id "[^\"]+"',lines) &
               grepl('transcript_id "[^\"]+"',lines)
    lines <- lines[has_ids]
    if (!length(lines)) next
    transcript <- sub('.*transcript_id "([^\"]+)".*','\\1',lines)
    gene <- sub('.*gene_id "([^\"]+)".*','\\1',lines)
    pieces[[length(pieces)+1L]] <- data.frame(
      ensembl_transcript_id=normalize_mouse_ensembl(transcript),
      ensembl_gene_id=normalize_mouse_ensembl(gene),stringsAsFactors=FALSE)
  }
  if (!length(pieces)) stop("No transcript/gene annotation records in GTF: ",gtf_file)
  mapping <- unique(do.call(rbind,pieces))
  if (any(!grepl("^ENSMUST[0-9]+$",mapping$ensembl_transcript_id)) ||
      any(!grepl("^ENSMUSG[0-9]+$",mapping$ensembl_gene_id)))
    stop("GTF does not contain the expected mouse Ensembl transcript/gene IDs.")
  mapping
}


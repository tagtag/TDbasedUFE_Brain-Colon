study_annotation <- function(path = "reference/boruta/Ensembl2GeneSymbol.csv") {
  x <- read.csv(path, header = FALSE, col.names = c("gene_id", "transcript_id", "gene_symbol"),
    stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
  unique(x)
}
ens_enrichment <- function(selected_ids, background_ids, out) {
  need("readxl")
  B <- unique(background_ids); S <- unique(selected_ids)
  stopifnot(all(S %in% B))
  anno <- unique(study_annotation()[, c("gene_id", "gene_symbol")])
  anno <- anno[anno$gene_id %in% B & !is.na(anno$gene_symbol), ]
  stopifnot(setequal(anno$gene_id, B), !anyDuplicated(anno$gene_id))
  d <- as.data.frame(readxl::read_excel("reference/ens/Drokhlyansky_Table_S4.xlsx",
    sheet = "Mouse Colon All"))
  f <- d[!is.na(d$padjD) & !is.na(d$coefD) & d$padjD < 0.05 & d$coefD > 0, ]
  csv_write(f, file.path(out, "ENS_filtered_TableS4.csv"))
  split_words <- function(s) strsplit(s, " ", fixed = TRUE)[[1]]
  neuron <- split_words(paste("Tubb3 Elavl4 Ret Phox2b Chrnb4 Eml5 Smpd3 Tagln3",
    "Snap25 Gpr22 Gdap1l1 Stmn3 Chrna3 Scg3 Syt4 Ncan Crmp1 Adcyap1r1 Elavl3 Dlg2 Cacna2d"))
  glia <- split_words(paste("Erbb3 Sox10 Fabp7 Plp1 Gas7 Nid1 Qk Sparc Mest Nfia Wwtr1",
    "Gpm6b Rasa3 Flrt1 Itpripl1 Itga4 Lama4 Postn Ptprz1 Pdpn Col18a1 Nrcam"))
  sets <- list(Published_Neuron = neuron, Published_Glia = glia,
    Published_ENS_union = union(neuron, glia))
  for (category in sort(unique(d$ident)))
    sets[[paste0("S4_", category)]] <- unique(f$gene[f$ident == category])
  stopifnot(length(sets) == 32L)
  rows <- members <- list()
  for (nm in names(sets)) {
    symbols <- unique(sets[[nm]])
    M <- unique(anno$gene_id[anno$gene_symbol %in% symbols])
    a <- length(intersect(S, M)); b <- length(S)-a; c <- length(M)-a
    dd <- length(B)-a-b-c
    tab <- matrix(c(a,b,c,dd), nrow = 2, byrow = TRUE)
    test <- fisher.test(tab, alternative = "greater")
    rows[[nm]] <- data.frame(marker_set = nm, background_n = length(B),
      selected_n = length(S), marker_genes_in_background = length(M), overlap_n = a,
      expected_overlap = length(S)*length(M)/length(B),
      fold_enrichment = if(length(M)) a/(length(S)*length(M)/length(B)) else NA_real_,
      sample_odds_ratio = (a*dd)/(b*c), conditional_odds_ratio = unname(test$estimate),
      p_value = test$p.value, a = a, b = b, c = c, d = dd,
      overlap_symbols = paste(sort(unique(anno$gene_symbol[anno$gene_id %in% intersect(S,M)])), collapse = ";"),
      unmatched_symbols = paste(sort(setdiff(symbols,anno$gene_symbol)), collapse = ";"))
    for (s in symbols) {
      gs <- anno$gene_id[anno$gene_symbol == s]
      if (!length(gs)) gs <- NA_character_
      members[[length(members)+1L]] <- data.frame(marker_set = nm,
        gene_symbol = s, gene_id = gs, selected = gs %in% S)
    }
  }
  result <- do.call(rbind, rows)
  result$BH_adjusted_p <- p.adjust(result$p_value, "BH")
  csv_write(result, file.path(out, "ENS_all_32_tests.csv"))
  csv_write(do.call(rbind,members), file.path(out,"ENS_marker_membership.csv"))
  csv_write(result[match(c("S4_Neuron_1","S4_Neuron_2","S4_Glia","S4_Fibroblast"),result$marker_set), ],
    file.path(out, "Table4_ENS.csv"))
  invisible(result)
}
compare_boruta <- function(selected, out) {
  annotation <- study_annotation()
  # The original annotation lacks this transcript. This frozen, documented
  # Ensembl lookup was used in the manuscript comparison; do not query live REST.
  annotation <- unique(rbind(annotation, data.frame(gene_id = "ENSMUSG00000071984",
    transcript_id = "ENSMUST00000239425", gene_symbol = "Fndc1")))
  tx <- split(annotation$gene_id, annotation$transcript_id)
  records <- summaries <- list(); sets <- list()
  files <- sort(list.files("reference/boruta", pattern = "\\.RDS\\.csv$", full.names = TRUE))
  for (file in files) {
    d <- csv_read(file)
    id <- regmatches(d$feature, regexpr("ENSMUS[GT][0-9]+", d$feature))
    if (length(id) != nrow(d)) stop("Non-Ensembl feature in Boruta source: ", file)
    gid <- vapply(id, function(i) {
      if (startsWith(i,"ENSMUSG")) return(i)
      g <- unique(tx[[i]])
      if (length(g) != 1L) stop("Unmapped or ambiguous Boruta transcript: ", i)
      g
    }, character(1))
    name <- sub("\\.csv\\.RDS\\.csv$", "", basename(file))
    records[[name]] <- data.frame(source = name, feature = d$feature,
      count = as.integer(d$count), gene_id = gid)
    for (threshold in c(1L,10L)) {
      g <- unique(gid[d$count >= threshold])
      key <- paste(name, threshold, sep = "__ge")
      sets[[key]] <- g
      summaries[[key]] <- data.frame(source = name, minimum_runs = threshold,
        features = sum(d$count >= threshold), unique_genes = length(g),
        overlap_with_TD = sum(g %in% selected$gene_id),
        TD_not_selected = sum(!selected$gene_id %in% g))
    }
  }
  boruta <- unique(unlist(sets, use.names = FALSE))
  remaining <- selected[!selected$gene_id %in% boruta, ]
  excluded_symbols <- unique(toupper(annotation$gene_symbol[
    annotation$gene_id %in% boruta & !is.na(annotation$gene_symbol)]))
  csv_write(do.call(rbind,records),file.path(out,"Boruta_mapped_features.csv"))
  csv_write(do.call(rbind,summaries),file.path(out,"Boruta_overlap_summary.csv"))
  csv_write(data.frame(Boruta_genes=length(boruta),TD_genes=nrow(selected),
    overlap=sum(selected$gene_id %in% boruta),TD_only=nrow(remaining)),
    file.path(out,"Boruta_union_summary.csv"))
  csv_write(remaining,file.path(out,"TD_only_genes.csv"))
  csv_write(data.frame(gene_symbol=excluded_symbols),file.path(out,"Boruta_excluded_symbols.csv"))
  list(remaining = remaining, excluded_symbols = excluded_symbols)
}
read_gmt <- function(path) {
  lines <- strsplit(readLines(path,warn=FALSE), "\t", fixed=TRUE)
  sets <- list()
  for (a in lines) {
    if(length(a)<3L) stop("Invalid GMT line in ",path)
    g <- toupper(trimws(a[-c(1,2)])); sets[[a[1]]] <- unique(g[nzchar(g)])
  }
  sets
}
boruta_excluded_enrichment <- function(selected, comparison, out) {
  all_symbols <- unique(toupper(selected$gene_symbol))
  selected_symbols <- unique(toupper(comparison$remaining$gene_symbol))
  stopifnot(!anyNA(all_symbols),!anyNA(selected_symbols))
  rows <- settings <- list()
  libraries <- c("GO_Biological_Process_2026", "PanglaoDB_Augmented_2021")
  for (lib in libraries) {
    sets <- read_gmt(file.path("reference/gene_sets",paste0(lib,".gmt")))
    universe <- setdiff(unique(unlist(sets,use.names=FALSE)),comparison$excluded_symbols)
    Q <- intersect(selected_symbols,universe)
    rs <- lapply(names(sets),function(term) {
      M <- intersect(sets[[term]],universe); overlap <- intersect(Q,M)
      data.frame(library=lib,term=term,universe_genes=length(universe),query_genes=length(Q),
        term_genes=length(M),overlap=length(overlap),
        original_overlap=length(intersect(all_symbols,sets[[term]])),
        p_value=phyper(length(overlap)-1L,length(M),length(universe)-length(M),length(Q),lower.tail=FALSE),
        genes=paste(sort(overlap),collapse=";"))
    })
    rs <- do.call(rbind,rs); rs$BH_adjusted_p <- p.adjust(rs$p_value,"BH")
    rows[[lib]] <- rs
    settings[[lib]] <- data.frame(library=lib,universe=length(universe),query=length(Q),terms=length(sets))
  }
  result <- do.call(rbind,rows)
  csv_write(result,file.path(out,"TD_only_enrichment_all_terms.csv"))
  csv_write(do.call(rbind,settings),file.path(out,"TD_only_enrichment_settings.csv"))
  focus <- c("Modulation of Chemical Synaptic Transmission (GO:0050804)",
    "Chemical Synaptic Transmission (GO:0007268)","Synapse Organization (GO:0050808)",
    "Enteric Neurons","Enteric Glia Cells")
  idx <- match(tolower(focus),tolower(result$term))
  if(anyNA(idx)) stop("Expected Table 7 term missing from frozen GMT")
  table7 <- result[idx, ]; table7$overlap_after_before <- paste(table7$overlap,table7$original_overlap,sep="/")
  csv_write(table7,file.path(out,"Table7_TD_only_enrichment.csv"))
  invisible(result)
}
enrichr_tables <- function(out) {
  libs <- c("Orphanet_Augmented_2021","GTEx_Aging_Signatures_2021")
  for(i in seq_along(libs)) {
    d <- read.delim(file.path("reference/enrichr",paste0(libs[i],".tsv")),
      check.names=FALSE,quote="",comment.char="")
    # Keep the original API export order, as stated in the manuscript captions.
    csv_write(head(d,10),file.path(out,paste0("Table",i,"_",libs[i],".csv")))
  }
}

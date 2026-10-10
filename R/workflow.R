# Shared configuration and input acquisition. All operations are sequential.
ages <- c("2M", "9M", "15M", "24M", "30M")
layer_files <- c(
  hippocampus_RNA = "varianceStabilizedRNASeq_Brain.tsv",
  hippocampus_methylation = "brainMeanPromoterMethylationNoMissingValues.tsv",
  colon_RNA = "varianceStabilizedRNASeq_Colon.tsv",
  colon_methylation = "intestineMeanPromoterMethylationNoMissingValues.tsv")
behavior_file <- "df_allCognitionData_2025-02-17.tsv"
behavior_columns <- c("TrainingPL_LogRegSlope", "TrainingPL_LogRegYIntercept",
  "Probetrial_CognitiveScore", "Retention_CognitiveScore", "TrainingCS_LinRegSlope",
  "TrainingCS_LinRegYIntercept", "InvertedProbetrialPL_Day7",
  "InvertedRetentionPL_BestOf2", "MedianRankCombinedScore")
csv_read <- function(path) read.csv(path, stringsAsFactors = FALSE,
  check.names = FALSE, fileEncoding = "UTF-8-BOM")
csv_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
}
need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Install required package: ", pkg, ". See scripts/install_dependencies.R")
}
initialize <- function(out) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
}
download_once <- function(url, dest) {
  if (file.exists(dest)) return(invisible(dest))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  tmp <- paste0(dest, ".part")
  message("Downloading: ", url)
  status <- download.file(url, tmp, method = "libcurl", mode = "wb")
  if (status != 0L || !file.exists(tmp) || file.info(tmp)$size == 0)
    stop("Download failed; no automatic retry. ", url)
  if (!file.rename(tmp, dest)) stop("Could not finalize download: ", dest)
  invisible(dest)
}
find_input <- function(name) {
  files <- list.files("data/source", recursive = TRUE, full.names = TRUE)
  hits <- files[basename(files) == name]
  if (length(hits) != 1L) stop("Expected exactly one source input: ", name)
  hits
}
download_inputs <- function() {
  archive <- "data/ScriptsAndFilesCognitiveAgingManuscript.zip"
  download_once(paste0("https://zenodo.org/records/15849341/files/",
    "ScriptsAndFilesCognitiveAgingManuscript.zip?download=1"), archive)
  entries <- unzip(archive, list = TRUE)$Name
  if (any(grepl("(^/|(^|/)\\.\\.(/|$)|^[A-Za-z]:)", entries)))
    stop("Unsafe archive paths")
  requested <- c(unname(layer_files), behavior_file, "Ensembl2GeneSymbol.csv",
    "DATA_LICENSE.txt", basename(list.files("reference/boruta", pattern = "RDS.csv$")))
  for (name in requested) {
    hit <- entries[basename(entries) == name]
    if (length(hit) != 1L) stop("Missing or ambiguous archive member: ", name)
    target <- file.path("data/source", hit)
    if (!file.exists(target)) unzip(archive, files = hit, exdir = "data/source")
    if (!file.exists(target)) stop("Extraction failed: ", name)
  }
  download_once(paste0("https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/",
    "Mus_musculus.GRCm39.116.gtf.gz"), "data/Mus_musculus.GRCm39.116.gtf.gz")
  paths <- c(archive, "data/Mus_musculus.GRCm39.116.gtf.gz",
    vapply(requested, find_input, character(1)))
  csv_write(data.frame(file = paths, md5 = unname(tools::md5sum(paths))),
    "data/downloaded_file_checksums.csv")
}
full_header <- function(path) {
  lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(lines) != 2L) stop("Too few input lines: ", path)
  h <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  row <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  if (length(row) == length(h) + 1L) h <- c("feature_id", h)
  else if (length(row) != length(h)) stop("Unexpected TSV header: ", path)
  h <- trimws(h)
  if (anyDuplicated(h)) stop("Duplicate header labels: ", path)
  h
}
prepare_samples <- function(out) {
  ids <- lapply(layer_files, function(f) full_header(find_input(f))[-1])
  h <- full_header(find_input(behavior_file))
  behavior <- read.delim(find_input(behavior_file), header = FALSE, skip = 1,
    col.names = h, colClasses = c("character", rep("numeric", length(h)-1L)),
    row.names = NULL, check.names = FALSE, quote = "", comment.char = "",
    na.strings = c("NA", "NaN", ""))
  names(behavior)[1] <- "sample_id"
  behavior$sample_id <- trimws(behavior$sample_id)
  stopifnot(!anyDuplicated(behavior$sample_id),
    setequal(names(behavior)[-1], behavior_columns))
  ids$behavior <- behavior$sample_id
  common <- Reduce(intersect, ids)
  complete <- apply(is.finite(as.matrix(behavior[, behavior_columns])), 1, all)
  eligible <- setdiff(intersect(common, behavior$sample_id[complete]), "30M_16N")
  # Use the frozen publication map, not a new random sample on each run.
  sm <- csv_read("reference/sample_id_map.csv")
  stopifnot(nrow(sm) == 75L, !anyDuplicated(sm$sample_id),
    all(sm$sample_id %in% eligible), all(table(sm$age) == 15L))
  count <- function(x) vapply(ages, function(a) sum(startsWith(x, paste0(a,"_"))), integer(1))
  t5 <- do.call(rbind, lapply(names(ids), function(nm)
    data.frame(dataset = nm, t(count(ids[[nm]])), total = length(ids[[nm]]), check.names = FALSE)))
  t6 <- data.frame(age = ages, common = count(common), eligible = count(eligible),
    selected = count(sm$sample_id))
  t6 <- rbind(t6, data.frame(age = "Total", common = length(common),
    eligible = length(eligible), selected = nrow(sm)))
  csv_write(t5, file.path(out, "Table5_source_samples.csv"))
  csv_write(t6, file.path(out, "Table6_sample_selection.csv"))
  matched <- cbind(sm, behavior[match(sm$sample_id, behavior$sample_id), behavior_columns])
  csv_write(matched, file.path(out, "sample_behavior_matched.csv"))
  csv_write(sm, file.path(out, "sample_id_map.csv"))
  stopifnot(identical(as.integer(t6$common), c(15L,16L,16L,17L,17L,81L)),
    identical(as.integer(t6$eligible), c(15L,15L,16L,17L,15L,78L)))
  matched
}

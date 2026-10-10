# Common-gene tensor: tissue x modality x age x replicate x gene.
# Methylation: arithmetic mean of signed/raw promoter values per ENSMUSG gene.
# RNA: published VST values. No preliminary SVD, no component sign alignment.
# Reuse the existing 75 sample IDs and transcript->gene mapping, without downloads.
# Default normalization: center/scale across common genes for each tissue/modality/sample.
# All operations are sequential. No package installation or retries.
# Usage:
# source("build_cognitive_gene_tensor.R")
# G <- build_cognitive_gene_tensor()
# Z_gene <- G$Z
# HOSVD_gene <- G$HOSVD

build_cognitive_gene_tensor <- function(
    work_dir="data/source",
    sample_map_file="results/full/sample_id_map.csv",
    mapping_file="results/full/transcript_to_gene_map.csv",
    output_dir="results/full/tensor",
    input_normalization=c("sample_zscore","none")) {
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1",BLIS_NUM_THREADS="1")
  if (requireNamespace("RhpcBLASctl",quietly=TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  if (!requireNamespace("rTensor",quietly=TRUE)) stop("Package rTensor is required.")
  input_normalization <- match.arg(input_normalization)
  normalize_id <- function(x) sub("\\.[0-9]+$","",trimws(as.character(x)))
  ages <- c("2M","9M","15M","24M","30M")
  tissues <- c("hippocampus","colon")
  modalities <- c("RNA","methylation")
  sm <- read.csv(sample_map_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!all(c("age","replicate","sample_id") %in% names(sm)))
    stop("Sample map needs age, replicate, sample_id.")
  if (nrow(sm)!=75L || anyNA(sm[,c("age","replicate","sample_id")]) ||
      anyDuplicated(sm$sample_id) || any(!sm$age %in% ages) ||
      any(!sm$replicate %in% 1:15)) stop("Invalid 75-individual sample map.")
  grid <- expand.grid(age=ages,replicate=1:15,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  key <- paste(sm$age,sm$replicate,sep="::")
  if (anyDuplicated(key)) stop("Duplicate age/replicate pair.")
  idx <- match(paste(grid$age,grid$replicate,sep="::"),key)
  if (anyNA(idx)) stop("Missing age/replicate pair.")
  sm <- sm[idx,,drop=FALSE]
  sm$sample_index <- seq_len(75L)
  sample_ids <- sm$sample_id
  if (any(!startsWith(sample_ids,paste0(sm$age,"_")))) stop("Sample IDs and age labels differ.")

  mapping <- read.csv(mapping_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!all(c("ensembl_transcript_id","ensembl_gene_id") %in% names(mapping)))
    stop("Mapping CSV needs ensembl_transcript_id, ensembl_gene_id.")
  map <- unique(data.frame(transcript_id=normalize_id(mapping$ensembl_transcript_id),
                          gene_id=normalize_id(mapping$ensembl_gene_id)))
  valid <- !is.na(map$transcript_id) & !is.na(map$gene_id) &
    grepl("^ENSMUST[0-9]+$",map$transcript_id) & grepl("^ENSMUSG[0-9]+$",map$gene_id)
  map <- map[valid,,drop=FALSE]
  multiplicity <- table(map$transcript_id)
  ambiguous <- names(multiplicity)[multiplicity>1L]
  map <- map[!map$transcript_id %in% ambiguous,,drop=FALSE]

  filenames <- c(hippocampus_RNA="varianceStabilizedRNASeq_Brain.tsv",
    colon_RNA="varianceStabilizedRNASeq_Colon.tsv",
    hippocampus_methylation="brainMeanPromoterMethylationNoMissingValues.tsv",
    colon_methylation="intestineMeanPromoterMethylationNoMissingValues.tsv")
  available <- list.files(work_dir,recursive=TRUE,full.names=TRUE)
  paths <- vapply(filenames,function(name) {
    hits <- available[basename(available)==name]
    if (length(hits)!=1L) stop("Expected one extracted input file: ",name)
    hits
  },character(1L))
  read_layer <- function(path) {
    lines <- readLines(path,n=2L,warn=FALSE)
    if (length(lines)!=2L) stop("Too few lines: ",path)
    h <- strsplit(lines[1],"\t",fixed=TRUE)[[1]]
    first <- strsplit(lines[2],"\t",fixed=TRUE)[[1]]
    if (length(first)==length(h)+1L) h <- c("feature_id",h)
    else if (length(first)!=length(h)) stop("Unexpected header: ",path)
    h <- trimws(h)
    if (anyDuplicated(h) || anyNA(match(sample_ids,h[-1]))) stop("Invalid sample columns: ",path)
    classes <- rep("NULL",length(h)); classes[1] <- "character"
    classes[match(sample_ids,h)] <- "numeric"
    tab <- read.delim(path,header=FALSE,skip=1L,col.names=h,colClasses=classes,
      row.names=NULL,check.names=FALSE,quote="",comment.char="",na.strings=c("NA","NaN",""))
    ids <- normalize_id(tab[[1]])
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid feature IDs: ",path)
    x <- as.matrix(tab[,sample_ids,drop=FALSE]); rownames(x) <- ids
    x
  }
  layers <- list(); statuses <- list(); counts <- list()
  for (layer in names(filenames)) {
    message("Reading and converting: ",layer)
    x <- read_layer(paths[[layer]])
    if (endsWith(layer,"_methylation")) {
      if (any(!grepl("^ENSMUST[0-9]+$",rownames(x)))) stop("Expected ENSMUST: ",layer)
      gene <- map$gene_id[match(rownames(x),map$transcript_id)]
      statuses[[layer]] <- data.frame(layer=layer,transcript_id=rownames(x),gene_id=gene,
        status=ifelse(is.na(gene),"unmapped_or_ambiguous","mapped"))
      ok <- !is.na(gene)
      if (!any(ok)) stop("No mapped genes: ",layer)
      # No transcript selection before averaging; no na.rm: incomplete genes excluded later.
      sums <- rowsum(x[ok,,drop=FALSE],group=gene[ok],reorder=TRUE,na.rm=FALSE)
      n <- table(gene[ok])
      x <- sweep(sums,1L,as.numeric(n[rownames(sums)]),"/")
      counts[[layer]] <- data.frame(layer=layer,gene_id=rownames(x),
        n_transcripts=as.integer(n[rownames(x)]))
    } else if (any(!grepl("^ENSMUSG[0-9]+$",rownames(x)))) stop("Expected ENSMUSG: ",layer)
    layers[[layer]] <- x
  }
  genes_before_filter <- sort(Reduce(intersect,lapply(layers,rownames)))
  if (length(genes_before_filter)<2L) stop("Too few genes shared by all four layers.")
  finite <- rep(TRUE,length(genes_before_filter))
  for (layer in names(layers))
    finite <- finite & rowSums(!is.finite(layers[[layer]][genes_before_filter,,drop=FALSE]))==0L
  genes <- genes_before_filter[finite]
  if (length(genes)<2L) stop("Too few finite common genes.")
  tensor <- array(NA_real_,dim=c(2L,2L,5L,15L,length(genes)),
    dimnames=list(tissue=tissues,modality=modalities,age=ages,
                  replicate=paste0("rep",1:15),gene=genes))
  normalization <- list()
  for (tissue in tissues) for (modality in modalities) {
    layer <- paste(tissue,modality,sep="_")
    # Rows = age (fastest), then replicate; columns = common genes.
    m <- t(layers[[layer]][genes,sample_ids,drop=FALSE])
    center <- rep(0,nrow(m)); scale <- rep(1,nrow(m))
    if (input_normalization=="sample_zscore") {
      center <- rowMeans(m); scale <- apply(m,1L,sd)
      if (any(!is.finite(scale) | scale<=0)) stop("Cannot standardize: ",layer)
      m <- sweep(sweep(m,1L,center,"-"),1L,scale,"/")
    }
    tensor[tissue,modality,,,] <- array(m,dim=c(5L,15L,length(genes)))
    normalization[[layer]] <- data.frame(layer=layer,sample_id=sample_ids,
                                         center=center,scale=scale)
  }
  rm(layers); gc(verbose=FALSE)
  Z <- list(tensor=tensor,genes=genes,sample_map=sm,
    settings=list(input_normalization=input_normalization,
      methylation_aggregation="arithmetic mean per gene, before intersection and normalization",
      mapping_file=mapping_file,sample_map_file=sample_map_file,
      genes_common_before_finite_filter=length(genes_before_filter),
      genes_removed_nonfinite=sum(!finite)),
    normalization=do.call(rbind,normalization))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(Z,file.path(output_dir,"cognitive_gene_tensor.rds"))
  write.csv(sm,file.path(output_dir,"sample_id_map.csv"),row.names=FALSE)
  write.csv(data.frame(gene_id=genes),file.path(output_dir,"common_genes.csv"),row.names=FALSE)
  write.csv(do.call(rbind,statuses),file.path(output_dir,"methylation_mapping_status.csv"),row.names=FALSE)
  write.csv(do.call(rbind,counts),file.path(output_dir,"methylation_gene_transcript_counts.csv"),row.names=FALSE)
  write.csv(Z$normalization,file.path(output_dir,"normalization.csv"),row.names=FALSE)
  # Gene unfolding has N rows and 300 columns: at most 300 nonzero directions.
  # Retain ALL possible nonzero directions, avoiding an unnecessary N x N basis.
  ranks <- c(2L,2L,5L,15L,min(length(genes),300L))
  message("HOSVD of tensor: ",paste(dim(tensor),collapse=" x "))
  HOSVD <- rTensor::hosvd(rTensor::as.tensor(tensor),ranks=ranks)
  rownames(HOSVD$U[[5]]) <- genes
  saveRDS(HOSVD,file.path(output_dir,"cognitive_gene_HOSVD.rds"))
  message("Saved tensor and HOSVD in ",output_dir)
  invisible(list(Z=Z,HOSVD=HOSVD))
}

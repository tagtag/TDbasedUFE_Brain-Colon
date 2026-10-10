run_tensor_analysis <- function(matched, out) {
  source("R/gtf_mapping.R")
  source("R/build_cognitive_gene_tensor.R")
  source("R/select_cognitive_gene_features.R")
  source("R/ensembl_mouse_to_symbol.R")
  gtf <- "data/Mus_musculus.GRCm39.116.gtf.gz"
  # The original build reused overlap_1121/transcript_to_gene_map.csv.
  # That table was restricted to the finite transcript intersection of both
  # methylation layers in the earlier mapping step. Reconstruct its scope
  # without rerunning that earlier SVD/feature-selection analysis.
  if(file.exists("data/transcript_to_gene_map.csv")) {
    mapping <- csv_read("data/transcript_to_gene_map.csv")
  } else {
    mapping <- read_mouse_gtf_mapping(gtf)
    finite_ids <- list()
    for(layer in c("hippocampus_methylation","colon_methylation")) {
      path <- find_input(layer_files[[layer]])
      h <- full_header(path)
      classes <- rep("NULL",length(h)); classes[1] <- "character"
      sample_cols <- match(matched$sample_id,h)
      if(anyNA(sample_cols)) stop("Missing methylation sample during mapping reconstruction")
      classes[sample_cols] <- "numeric"
      d <- read.delim(path,header=FALSE,skip=1,col.names=h,colClasses=classes,
        row.names=NULL,check.names=FALSE,quote="",comment.char="",na.strings=c("NA","NaN",""))
      ids <- normalize_mouse_ensembl(d[[1]])
      if(anyDuplicated(ids)) stop("Duplicate methylation transcript IDs")
      x <- as.matrix(d[,matched$sample_id,drop=FALSE])
      finite_ids[[layer]] <- ids[rowSums(!is.finite(x))==0L]
    }
    common_transcripts <- Reduce(intersect,finite_ids)
    mapping <- unique(mapping[mapping$ensembl_transcript_id %in% common_transcripts, ])
    n_map <- table(mapping$ensembl_transcript_id)
    mapping <- mapping[!mapping$ensembl_transcript_id %in% names(n_map)[n_map>1L], ]
    rm(d,x,finite_ids)
  }
  mapping_file <- file.path(out,"transcript_to_gene_map.csv")
  csv_write(mapping,mapping_file)
  built <- build_cognitive_gene_tensor(work_dir="data/source",
    sample_map_file=file.path(out,"sample_id_map.csv"),mapping_file=mapping_file,
    output_dir=file.path(out,"tensor"),input_normalization="sample_zscore")
  Z <- built$Z; H <- built$HOSVD
  reference_bg <- csv_read("reference/background_gene_ids.csv")[[1]]
  if (!setequal(Z$genes,reference_bg))
    stop("Common gene background differs from the manuscript. Inspect mapping and inputs; no thresholds were changed.")
  # Indices selected from tissue/modality/age/replicate patterns, not behavior.
  association <- abs(H$Z@data[1,1,2,1,])
  component <- which.max(association)
  csv_write(data.frame(gene_component=seq_along(association),absolute_core=association),
    file.path(out,"core_component_association.csv"))
  if (component != 5L) stop("Expected gene component 5, obtained ",component,
    ". Inspect decomposition rather than choosing a component using cognition.")
  selected <- select_cognitive_gene_features(gene_component=component,HOSVD_gene=H,
    output_dir=file.path(out,"selection"),alpha=0.01,bins=100L,
    min_background_fraction=0.5,coarse_points=101L,refine_points=61L,refinement_rounds=2L)
  ref_selected <- csv_read("reference/selected_gene_symbols.csv")
  if (!setequal(selected$selected_ids,ref_selected$gene_id))
    stop("Selected gene IDs differ from the 1,577 manuscript IDs. Inspect saved diagnostics; no automatic retuning.")
  symbols <- ensembl_mouse_to_symbol(selected$selected_ids,gtf_file=gtf,
    output_file=file.path(out,"selected_gene_symbols.csv"))
  if(any(symbols$status != "mapped")) stop("Unmapped selected gene symbol")
  ref_sym <- ref_selected$gene_symbol[match(symbols$gene_id,ref_selected$gene_id)]
  if(!identical(symbols$gene_symbol,ref_sym)) stop("Selected gene symbol mapping differs from frozen reference")
  gene_vector <- H$U[[5]][,component]
  scores <- matched
  for(tissue in dimnames(Z$tensor)[[1]]) for(modality in dimnames(Z$tensor)[[2]]) {
    mat <- matrix(Z$tensor[tissue,modality,,,],nrow=75L)
    scores[[paste(tissue,modality,sep="_")]] <- as.vector(mat %*% gene_vector)
  }
  # SVD signs are arbitrary. Match the molecular score orientation in the
  # archived run, without consulting any behavioral column.
  ref_scores <- csv_read("reference/molecular_scores.csv")
  ref_scores <- ref_scores[match(scores$sample_id,ref_scores$sample_id), ]
  sign_match <- sign(sum(as.matrix(scores[,names(layer_files)]) *
    as.matrix(ref_scores[,names(layer_files)])))
  if(sign_match == 0) stop("Unable to orient molecular component")
  scores[,names(layer_files)] <- scores[,names(layer_files)] * sign_match
  csv_write(scores,file.path(out,"molecular_scores.csv"))
  csv_write(data.frame(gene_id=Z$genes,loading=gene_vector*sign_match),
    file.path(out,"gene_component5_loadings.csv"))
  csv_write(data.frame(component=component,sign_relative_to_raw_svd=sign_match),
    file.path(out,"component_orientation.csv"))
  # A plotting copy: orient factors for readability; leave raw HOSVD intact.
  pdf(file.path(out,"Figure1_mode_factors.pdf"),width=10,height=7)
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  indices <- c(1L,1L,2L,1L)
  labels <- list(c("Hippocampus","Colon"),c("RNA","Methylation"),ages,as.character(1:15))
  titles <- c("Tissue: component 1","Modality: component 1","Age: component 2","Replicate: component 1")
  for(mode in 1:4) {
    v <- H$U[[mode]][,indices[mode]]
    if(v[1]<0) v <- -v
    csv_write(data.frame(level=labels[[mode]],loading=v),
      file.path(out,paste0("Figure1_mode",mode,".csv")))
    plot(seq_along(v),v,type="b",pch=19,xaxt="n",xlab="",ylab="Loading",main=titles[mode])
    axis(1,at=seq_along(v),labels=labels[[mode]],cex.axis=0.8)
  }
  dev.off()
  pdf(file.path(out,"Figure2_gene_loadings.pdf"),width=8,height=5)
  plot(gene_vector*sign_match,pch=16,cex=0.3,col=ifelse(Z$genes %in% selected$selected_ids,"#D55E00","grey65"),
    xlab="Gene index (sorted Ensembl ID)",ylab="Gene component 5 loading")
  dev.off()
  list(scores=scores,selected=symbols,background=Z$genes)
}

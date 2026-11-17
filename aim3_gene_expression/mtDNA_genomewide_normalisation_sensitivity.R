#!/usr/bin/env Rscript
# Sensitivity analysis: personalised mitochondrial counts with nuclear-derived
# STAR whole-transcriptome size factors, compared with mtDNA-only DESeq2 factors.
# Run: Rscript mtDNA_genomewide_normalisation_sensitivity.R
suppressPackageStartupMessages(library(DESeq2))

ROOT <- "/home/dayna/Aim_3"
EXPR <- file.path(ROOT, "mtdna_expression_matrix.tsv")
PE_META <- file.path(ROOT, "pe_metadata_unified.tsv")
SE_META <- file.path(ROOT, "se_metadata_unified.tsv")
GTF <- file.path(ROOT, "chrM_genes.gtf")
OUT <- file.path(ROOT, "normalisation_sensitivity")
STRAND_COL <- 4L  # R 1-indexed; STAR reverse-stranded column 4
MIN_NUCLEAR_SAMPLES <- 5L
MIN_NUCLEAR_COUNT <- 10L
# IMPORTANT: leave FALSE until original SE alignment settings are confirmed.
ALLOW_SE <- TRUE

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(EXPR), file.exists(PE_META), file.exists(SE_META), file.exists(GTF))
expr <- read.delim(EXPR, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(all(c("gene_id", "gene_name", "gene_type") %in% names(expr)))
gene_info <- expr[,c("gene_id","gene_name","gene_type")]
mat <- as.matrix(expr[,setdiff(names(expr),names(gene_info)),drop=FALSE])
rownames(mat) <- expr$gene_id
storage.mode(mat) <- "integer"
stopifnot(nrow(mat)==37L, !anyDuplicated(rownames(mat)))
mt_ids <- unique(expr$gene_id)
# Use all available nuclear genes in original STAR output, not personalised mt genes.
read_star <- function(sample, folder) {
  f <- file.path(ROOT, folder, sample, paste0(sample,"_ReadsPerGene.out.tab"))
  if (!file.exists(f)) return(NULL)
  d <- read.delim(f, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                  check.names=FALSE, quote="", comment.char="")
  if (ncol(d)!=4L || nrow(d)<10000L || !identical(as.character(d[1:4,1]),
       c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")))
    stop("Unexpected STAR format: ",f)
  v <- as.integer(d[[STRAND_COL]])
  names(v) <- d[[1]]
  v[-seq_len(4L)]
}
run_cohort <- function(label, meta_path, aligned_dir) {
  message("\n=== ",label," ===")
  meta <- read.delim(meta_path, stringsAsFactors=FALSE)
  stopifnot(all(c("sample_id","haplogroup","age","sex") %in% names(meta)),
            !anyDuplicated(meta$sample_id))
  ids <- intersect(meta$sample_id, colnames(mat))
  ids <- ids[colSums(is.na(mat[,ids,drop=FALSE]))==0L]
  raw <- lapply(ids, read_star, folder=aligned_dir)
  names(raw) <- ids
  missing <- names(raw)[vapply(raw,is.null,logical(1))]
  if(length(missing)) {
    writeLines(missing,file.path(OUT,paste0(label,"_missing_original_STAR.txt")))
    warning(label,": missing original STAR files for ",length(missing)," samples; see output")
  }
  raw <- raw[!vapply(raw,is.null,logical(1))]
  if(length(raw)<10L) stop("Too few matched original STAR files for ",label)
  common_genes <- Reduce(intersect,lapply(raw,names))
  if(length(common_genes)<10000L) stop("Original STAR gene IDs inconsistent")
  # Exclude mitochondrial genes using exact IDs and version-free ENSG roots.
  unversion <- function(x) sub("\\.[0-9]+$","",x)
  nuclear <- common_genes[!unversion(common_genes) %in% unversion(mt_ids)]
  # Also remove any genes whose name is absent from the reference matrix? No:
  # original STAR annotation is genome-wide and all non-mt genes are eligible.
  gmat <- vapply(raw,function(x) as.integer(x[nuclear]),integer(length(nuclear)))
  rownames(gmat) <- nuclear
  if(anyNA(gmat)) stop("NA in genome-wide count matrix")
  keep <- rowSums(gmat>=MIN_NUCLEAR_COUNT)>=MIN_NUCLEAR_SAMPLES
  gmat <- gmat[keep,,drop=FALSE]
  if(nrow(gmat)<1000L) stop("Too few expressed nuclear genes")
  # DESeq2 poscounts handles zeros; estimate from nuclear genes ONLY.
  sf <- estimateSizeFactorsForMatrix(gmat,type="poscounts")
  sf <- sf/exp(mean(log(sf)))  # geometric mean = 1
  stopifnot(all(is.finite(sf)),all(sf>0))
  ids <- colnames(gmat)
  metadata <- meta[match(ids,meta$sample_id),,drop=FALSE]
  rownames(metadata) <- ids
  metadata$L_status <- factor(ifelse(startsWith(metadata$haplogroup,"L"),"L","non_L"),
                              levels=c("non_L","L"))
  metadata$L_subgroup <- factor(ifelse(startsWith(metadata$haplogroup,"L"),
                                    substr(metadata$haplogroup,1,2),"non_L"))
  metadata$L_subgroup <- relevel(metadata$L_subgroup,ref="non_L")
  write.table(data.frame(sample_id=ids,cohort=label,nuclear_size_factor=sf,
                         nuclear_library_counts=colSums(gmat),
                         mt_library_counts=colSums(mat[,ids,drop=FALSE])),
              file.path(OUT,paste0(label,"_size_factors.tsv")),
              sep="\t",quote=FALSE,row.names=FALSE)
  message(label,": ",length(ids)," samples; ",nrow(gmat)," expressed nuclear genes")
  rm(gmat,raw); gc(verbose=FALSE)

  results_all <- list()
  for (comparison in c("L0_vs_nonL","L_vs_nonL")) {
    if(comparison=="L0_vs_nonL") {
      sel <- metadata$L_subgroup %in% c("L0","non_L")
      group <- droplevels(metadata$L_subgroup[sel])
      design <- ~ group
      contrast <- c("group","L0","non_L")
    } else {
      sel <- rep(TRUE,nrow(metadata))
      group <- droplevels(metadata$L_status[sel])
      design <- ~ group
      contrast <- c("group","L","non_L")
    }
    sel_ids <- rownames(metadata)[sel]
    if(any(table(group)<5)) stop("Insufficient group sizes: ",comparison)
    count_subset <- mat[,sel_ids,drop=FALSE]
    cd <- data.frame(group=group,row.names=sel_ids)
    for(mode in c("mtDNA_only","nuclear_reference")) {
      dds <- DESeqDataSetFromMatrix(countData=count_subset,colData=cd,design=design)
      if(mode=="nuclear_reference") {
        sf_sel <- sf[match(sel_ids,ids)]
        # Recenter within exact comparison subset for interpretable count scale.
        sizeFactors(dds) <- sf_sel/exp(mean(log(sf_sel)))
      }
      # Existing analysis uses default mtDNA-only size factors.
      dds <- DESeq(dds,quiet=TRUE)
      res <- as.data.frame(results(dds,contrast=contrast))
      res$gene_id <- rownames(res)
      res <- merge(gene_info,res,by="gene_id",sort=FALSE)
      res$cohort <- label; res$comparison <- comparison; res$normalisation <- mode
      results_all[[paste(comparison,mode,sep="_")]] <- res
      write.table(res,file.path(OUT,paste(label,comparison,mode,"results.tsv",sep="_")),
                  sep="\t",quote=FALSE,row.names=FALSE)
      if(mode=="mtDNA_only") {
        mt_sf <- sizeFactors(dds)
        write.table(data.frame(sample_id=names(mt_sf),mtDNA_size_factor=mt_sf),
                    file.path(OUT,paste(label,comparison,"mtDNA_size_factors.tsv",sep="_")),
                    sep="\t",quote=FALSE,row.names=FALSE)
      }
      message(label," ",comparison," ",mode,": ",sum(res$padj<0.05,na.rm=TRUE)," FDR-significant genes")
    }
  }
  do.call(rbind,results_all)
}
pe <- run_cohort("GC6_74",PE_META,"aligned_PE")
se <- if(ALLOW_SE) run_cohort("TANDEM",SE_META,"aligned_SE") else NULL
all <- rbind(pe,se)
write.table(all,file.path(OUT,"all_sensitivity_results.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
# Side-by-side matched comparison and independent replication by FDR and direction.
key <- c("gene_name","cohort","comparison")
wide <- reshape(all[,c(key,"normalisation","log2FoldChange","pvalue","padj")],
                idvar=key,timevar="normalisation",direction="wide")
write.table(wide,file.path(OUT,"original_vs_nuclear_reference.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
if(ALLOW_SE) {
  for(mode in c("mtDNA_only","nuclear_reference")) for(comp in c("L0_vs_nonL","L_vs_nonL")) {
    d <- all[all$normalisation==mode & all$comparison==comp,]
    a <- d[d$cohort=="GC6_74",c("gene_name","log2FoldChange","padj")]
    b <- d[d$cohort=="TANDEM",c("gene_name","log2FoldChange","padj")]
    m <- merge(a,b,by="gene_name",suffixes=c("_GC6_74","_TANDEM"))
    m$replicated <- !is.na(m$padj_GC6_74) & !is.na(m$padj_TANDEM) &
      m$padj_GC6_74<0.05 & m$padj_TANDEM<0.05 &
      sign(m$log2FoldChange_GC6_74)==sign(m$log2FoldChange_TANDEM)
    write.table(m,file.path(OUT,paste("replication",comp,mode,"tsv",sep=".")),
                sep="\t",quote=FALSE,row.names=FALSE)
    message("Replication ",comp," ",mode,": ",sum(m$replicated)," genes")
  }
}
message("\nDone. Outputs in ",OUT)


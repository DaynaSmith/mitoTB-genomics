#!/usr/bin/env Rscript
# FINAL nuclear eQTL pipeline: PE + SE, sex REQUIRED in every model.
# Discovery, replication, haplogroup sensitivity, QC, thesis tables in ONE run.
# Run from directory containing the four input TSV files.
# Requires: MatrixEQTL, DESeq2. Does not modify original analyses.
# Note: nuclear eQTL only; does NOT rerun mtDNA-variant eQTLs.

suppressPackageStartupMessages({
  if (!requireNamespace('MatrixEQTL', quietly=TRUE)) stop('Install MatrixEQTL first')
  if (!requireNamespace('DESeq2', quietly=TRUE)) stop('Install DESeq2 first')
  library(MatrixEQTL); library(DESeq2)
})
EXPR_FILE <- 'mtdna_expression_matrix.tsv'
GENO_FILE <- 'ndna_strictQC_noNUMT_genotype_matrix.tsv'
META_FILES <- c(PE='pe_metadata_unified.tsv', SE='se_metadata_unified.tsv')
OUTDIR <- 'FINAL_nuclear_eQTL_SEX_INCLUDED'
MIN_MAC <- 10L; MIN_PRESENCE <- 0.90; N_PCS <- 5L
ALPHA <- 0.05
# QC is recalculated on the sex-complete cohort to match the analysed samples.
# PCs are also recalculated on that cohort (not inherited from 405-sample PE PCA).
for (f in c(EXPR_FILE, GENO_FILE, META_FILES)) if (!file.exists(f)) stop('Missing input: ', f)
if (dir.exists(OUTDIR)) stop('Output directory already exists: ', OUTDIR,
                             '. Rename/move it before a fresh run to avoid mixing results.')
dir.create(OUTDIR, recursive=TRUE)
logcon <- file(file.path(OUTDIR,'run_log.txt'),open='wt')
logmsg <- function(...) { s <- paste0(...,collapse=''); cat(s,'\n'); writeLines(s,logcon); flush(logcon) }
write_tsv <- function(x, filename) write.table(x,file.path(OUTDIR,filename),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
read_tsv <- function(x, ...) read.delim(x,check.names=FALSE,stringsAsFactors=FALSE,...)

logmsg('Loading expression counts and strict-QC NUMT-excluded nuclear genotype matrix')
expr <- read_tsv(EXPR_FILE)
if (!all(c('gene_id','gene_name','gene_type') %in% names(expr))) stop('Expression file missing gene annotations')
if (anyDuplicated(expr$gene_id)) stop('Duplicate expression gene IDs')
count_cols <- setdiff(names(expr),c('gene_id','gene_name','gene_type'))
raw_counts <- as.matrix(expr[,count_cols,drop=FALSE]); rownames(raw_counts) <- expr$gene_id
suppressWarnings(storage.mode(raw_counts) <- 'numeric')
if (any(!is.na(raw_counts) & (raw_counts < 0 | raw_counts != floor(raw_counts)))) stop('Expression counts not nonnegative integers')
# Original normalization: exclude any expression sample with missing gene counts,
# normalize globally, and log2(normalized count+1).
raw_counts <- raw_counts[,colSums(is.na(raw_counts))==0,drop=FALSE]
storage.mode(raw_counts) <- 'integer'
dds <- DESeqDataSetFromMatrix(countData=raw_counts,
                              colData=data.frame(row.names=colnames(raw_counts)),design=~1)
dds <- estimateSizeFactors(dds)
norm_expr <- log2(DESeq2::counts(dds,normalized=TRUE)+1)
rm(dds)
geno <- as.matrix(read_tsv(GENO_FILE,row.names=1))
suppressWarnings(storage.mode(geno) <- 'numeric')
if (anyDuplicated(rownames(geno)) || anyDuplicated(colnames(geno))) stop('Duplicate genotype SNP/sample identifiers')
if (any(!is.na(geno) & !(geno %in% 0:2))) stop('Genotype data must be 0/1/2 dosage or NA')
logmsg('Input: ',nrow(geno),' SNPs; ',ncol(geno),' genotype samples; ',ncol(norm_expr),' expression samples')

parse_hap <- function(x) {
  hap <- trimws(as.character(x)); hap[is.na(hap) | hap=='' | toupper(hap) %in% c('NA','UNKNOWN','UNASSIGNED')] <- NA_character_
  # Valid L0-L6 vs non-L. Do not classify unknown as non-L.
  isL <- grepl('^L[0-6]',hap)
  isNonL <- !is.na(hap) & !isL
  lstatus <- ifelse(is.na(hap),NA_character_,ifelse(isL,'L','non-L'))
  major <- ifelse(is.na(hap),NA_character_,ifelse(isL,substr(hap,1,2),'non-L'))
  major[major %in% c('L4','L5','L6')] <- 'Other L'
  list(L=factor(lstatus,levels=c('non-L','L')),
       major=factor(major,levels=c('non-L','L0','L1','L2','L3','Other L')))
}

cohorts <- list()
for (label in names(META_FILES)) {
  logmsg('\n===== ',label,' : cohort preparation =====')
  meta <- read_tsv(META_FILES[[label]])
  if (!all(c('sample_id','sex','haplogroup') %in% names(meta))) stop(label,': metadata needs sample_id, sex, haplogroup')
  if (anyDuplicated(meta$sample_id)) stop(label,': duplicate metadata sample IDs')
  samples <- Reduce(intersect,list(as.character(meta$sample_id),colnames(geno),colnames(norm_expr)))
  m <- meta[match(samples,meta$sample_id),,drop=FALSE]
  sex <- trimws(as.character(m$sex)); sex[is.na(sex) | sex==''] <- NA_character_
  if (any(!is.na(sex) & !sex %in% c('Female','Male'))) stop(label,': unexpected sex categories: ',paste(unique(sex),collapse=', '))
  logmsg(label,': overlap before sex filter = ',length(samples),'; missing sex = ',sum(is.na(sex)))
  keep <- !is.na(sex)
  samples <- samples[keep]; m <- m[keep,,drop=FALSE]; sex <- sex[keep]
  if (length(samples) < 15 || length(unique(sex)) != 2) stop(label,': insufficient sex-complete data')
  m$sex_model <- as.numeric(factor(sex,levels=c('Female','Male'))) - 1L
  hap <- parse_hap(m$haplogroup)
  m$L_status <- hap$L; m$major_hap_model <- hap$major
  logmsg(label,': final sex-complete samples = ',length(samples))
  logmsg(label,': major haplogroups: ',paste(names(table(m$major_hap_model,useNA='ifany')),table(m$major_hap_model,useNA='ifany'),collapse='; '))
  # Gene presence from sex-complete analysis cohort, BEFORE normalization.
  gene_presence <- rowMeans(!is.na(raw_counts[,samples,drop=FALSE]) & raw_counts[,samples,drop=FALSE]>0)
  pass_genes <- names(gene_presence)[gene_presence >= MIN_PRESENCE]
  # MAC on exactly the sex-complete analysed samples.
  g <- geno[,samples,drop=FALSE]
  allele_sum <- rowSums(g,na.rm=TRUE)
  observed <- rowSums(!is.na(g))
  mac <- pmin(allele_sum,2*observed-allele_sum)
  pass_snps <- rownames(g)[observed>0 & mac>=MIN_MAC]
  logmsg(label,': SNPs MAC>=10 = ',length(pass_snps),'; genes presence>=0.90 = ',length(pass_genes))
  if (length(pass_snps)<N_PCS || !length(pass_genes)) stop(label,': no testable data')
  # PCA exactly follows original algorithm, on the current sex-complete sample set.
  pca_g <- g[pass_snps,,drop=FALSE]
  pca_g <- pca_g[rowSums(is.na(pca_g))==0,,drop=FALSE]
  v <- apply(pca_g,1,var)
  pca_g <- pca_g[is.finite(v)&v>0,,drop=FALSE]
  if (nrow(pca_g)<N_PCS) stop(label,': insufficient SNPs for PCA')
  if (nrow(pca_g)>50000L) {set.seed(42);pca_g<-pca_g[sample(nrow(pca_g),50000L),,drop=FALSE]}
  pca <- prcomp(t(pca_g),center=TRUE,scale.=TRUE)
  pcs <- pca$x[,seq_len(N_PCS),drop=FALSE]
  colnames(pcs) <- paste0('genoPC',seq_len(N_PCS))
  write_tsv(data.frame(sample_id=rownames(pcs),pcs,check.names=FALSE),paste0(label,'_genotype_pcs.tsv'))
  write_tsv(data.frame(sample_id=samples,sex=sex,haplogroup=m$haplogroup,
                       L_status=m$L_status,major_hap_model=m$major_hap_model),paste0(label,'_analysis_samples.tsv'))
  write_tsv(data.frame(snps=pass_snps),paste0(label,'_MAC10_pass_snps.tsv'))
  write_tsv(data.frame(gene=pass_genes),paste0(label,'_presence_pass_genes.tsv'))
  cohorts[[label]] <- list(samples=samples,meta=m,pcs=pcs,pass_snps=pass_snps,pass_genes=pass_genes)
  logmsg(label,': PCA on ',nrow(pca_g),' SNPs x ',ncol(pca_g),' samples')
}

run_discovery <- function(label) {
  c <- cohorts[[label]]; samples <- c$samples
  g <- geno[c$pass_snps,samples,drop=FALSE]
  y <- norm_expr[c$pass_genes,samples,drop=FALSE]
  cv <- rbind(sex=c$meta$sex_model,t(c$pcs[samples,,drop=FALSE]))
  stopifnot(identical(colnames(cv),samples),!anyNA(cv),nrow(cv)==6L)
  tests <- nrow(g)*nrow(y); threshold <- ALPHA/tests
  logmsg(label,': running MatrixEQTL modelLINEAR; n=',length(samples),
         '; SNPs=',nrow(g),'; genes=',nrow(y),'; tests=',tests,
         '; Bonferroni threshold=',format(threshold,scientific=TRUE),'; sex covariate=YES')
  snps <- SlicedData$new();snps$CreateFromMatrix(g);snps$ResliceCombined(sliceSize=5000)
  gene <- SlicedData$new();gene$CreateFromMatrix(y);gene$ResliceCombined(sliceSize=500)
  cvrt <- SlicedData$new();cvrt$CreateFromMatrix(cv)
  # Significant results are retained in memory; full results go to file.
  outfile <- file.path(OUTDIR,paste0(label,'_ALL_nuclear_eQTL_tests.tsv'))
  me <- Matrix_eQTL_engine(snps=snps,gene=gene,cvrt=cvrt,output_file_name=outfile,
                           pvOutputThreshold=1, useModel=modelLINEAR,errorCovariance=numeric(),
                           verbose=FALSE,pvalue.hist=FALSE)
  sig <- me$all$eqtls
  sig <- sig[is.finite(sig$pvalue) & sig$pvalue < threshold,,drop=FALSE]
  sig$bonferroni_sig <- TRUE
  sig$bonferroni_threshold <- threshold
  sig <- sig[order(sig$pvalue),,drop=FALSE]
  write_tsv(sig,paste0(label,'_SIGNIFICANT_nuclear_eQTL.tsv'))
  logmsg(label,': Bonferroni significant SNP-gene pairs = ',nrow(sig))
  list(sig=sig,threshold=threshold,tests=tests)
}
discovery <- lapply(names(cohorts),run_discovery); names(discovery)<-names(cohorts)

# Targeted evaluation of every PE significant pair in SE, even if not
# significant in SE. Use SE discovery-wide threshold for replication.
pe_sig <- discovery$PE$sig; se_sig <- discovery$SE$sig
replication <- data.frame()
if (nrow(pe_sig)) {
  se <- cohorts$SE
  pairs <- pe_sig[pe_sig$snps %in% se$pass_snps & pe_sig$gene %in% se$pass_genes,,drop=FALSE]
  logmsg('PE significant pairs evaluable in SE under SE QC: ',nrow(pairs))
  if (nrow(pairs)) {
    replication <- do.call(rbind,lapply(seq_len(nrow(pairs)),function(i) {
      p <- pairs[i,];s <- se$samples
      df <- data.frame(expression=as.numeric(norm_expr[p$gene,s]),
                       genotype=as.numeric(geno[p$snps,s]),sex=se$meta$sex_model,
                       se$pcs[s,,drop=FALSE],check.names=FALSE)
      df <- df[complete.cases(df),,drop=FALSE]
      fit <- lm(expression~genotype+sex+genoPC1+genoPC2+genoPC3+genoPC4+genoPC5,data=df)
      z <- summary(fit)$coefficients
      if (!'genotype' %in% rownames(z)) stop('Non-estimable replication genotype: ',p$snps)
      data.frame(gene=p$gene,snps=p$snps,PE_beta=p$beta,PE_p=p$pvalue,
                 SE_n=nrow(df),SE_beta=z['genotype','Estimate'],SE_p=z['genotype','Pr(>|t|)'],
                 same_direction=sign(p$beta)==sign(z['genotype','Estimate']),
                 replicated=sign(p$beta)==sign(z['genotype','Estimate']) &
                   z['genotype','Pr(>|t|)']<discovery$SE$threshold)
    }))
  }
}
write_tsv(replication,'PE_significant_pairs_evaluated_in_SE.tsv')

# Sensitivity for each cohort's newly discovered significant pairs.
# Identical sample sets across baseline, L/non-L, and major-haplogroup
# within each pair. Sex is REQUIRED in all three models.
fit_sensitivity <- function(label) {
  c <- cohorts[[label]]; sig <- discovery[[label]]$sig
  if (!nrow(sig)) {logmsg(label,': no significant pairs for sensitivity');return(data.frame())}
  out <- lapply(seq_len(nrow(sig)),function(i) {
    p <- sig[i,];s <- c$samples
    df <- data.frame(expression=as.numeric(norm_expr[p$gene,s]),
                     genotype=as.numeric(geno[p$snps,s]),sex=c$meta$sex_model,
                     c$pcs[s,,drop=FALSE],L_status=c$meta$L_status,
                     major_hap_model=c$meta$major_hap_model,check.names=FALSE)
    # Complete cases for ALL 3 models; no dropping sex.
    df <- droplevels(df[complete.cases(df),,drop=FALSE])
    if (nrow(df)<15 || length(unique(df$genotype))<2) stop(label,': too few sensitivity samples for ',p$snps,' / ',p$gene)
    f0 <- expression~genotype+sex+genoPC1+genoPC2+genoPC3+genoPC4+genoPC5
    fL <- update(f0,.~.+L_status)
    fM <- update(f0,.~.+major_hap_model)
    extract <- function(f) {
      m <- lm(f,data=df); co <- summary(m)$coefficients
      if (!'genotype' %in% rownames(co)) stop('Non-estimable genotype effect')
      c(beta=co['genotype','Estimate'],p=co['genotype','Pr(>|t|)'])
    }
    b <- extract(f0);l <- extract(fL);m <- extract(fM)
    data.frame(gene=p$gene,snps=p$snps,n=nrow(df),
               discovery_beta=p$beta,discovery_p=p$pvalue,
               baseline_beta=b['beta'],baseline_p=b['p'],
               L_beta=l['beta'],L_p=l['p'],major_beta=m['beta'],major_p=m['p'],
               reduction_abs_beta_L_pct=if (b['beta']!=0) 100*(abs(b['beta'])-abs(l['beta']))/abs(b['beta']) else NA_real_,
               reduction_abs_beta_major_pct=if (b['beta']!=0) 100*(abs(b['beta'])-abs(m['beta']))/abs(b['beta']) else NA_real_,
               major_sign_flip=sign(b['beta'])!=sign(m['beta']),
               discovery_threshold=discovery[[label]]$threshold)
  })
  ans <- do.call(rbind,out)
  write_tsv(ans,paste0(label,'_haplogroup_sensitivity_ALL_significant.tsv'))
  logmsg(label,': sensitivity pairs = ',nrow(ans),
         '; attenuated after major haplogroup = ',sum(ans$reduction_abs_beta_major_pct>0,na.rm=TRUE),
         '; baseline max abs beta discrepancy = ',max(abs(ans$baseline_beta-ans$discovery_beta),na.rm=TRUE))
  ans
}
sensitivity <- lapply(names(cohorts),fit_sensitivity); names(sensitivity)<-names(cohorts)

# One combined, publication-ready long-format table, without assuming
# significant counts or a fixed set of genes/variants.
summary <- do.call(rbind,lapply(names(cohorts),function(label) {
  c <- cohorts[[label]];d <- discovery[[label]];s <- sensitivity[[label]]
  data.frame(cohort=label,n=length(c$samples),sex_covariate='Included',
             n_snps=length(c$pass_snps),n_genes=length(c$pass_genes),
             n_tests=d$tests,bonferroni_threshold=d$threshold,
             n_significant=nrow(d$sig),n_attenuated_major=if(nrow(s))sum(s$reduction_abs_beta_major_pct>0) else 0)
}))
write_tsv(summary,'COHORT_SUMMARY_FOR_RESULTS.tsv')
write_tsv(data.frame(cohort='PE',sensitivity$PE), 'TABLE_4_16_PE_sex_adjusted.tsv')
write_tsv(data.frame(cohort='SE',sensitivity$SE), 'TABLE_TANDEM_SE_sex_adjusted.tsv')
logmsg('\n=== FINISHED: ALL MODELS INCLUDE SEX ===')
logmsg('Results directory: ',OUTDIR)
print(summary)
close(logcon)


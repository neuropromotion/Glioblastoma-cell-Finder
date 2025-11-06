library(Seurat) 

get_chromosome_means <- function(counts,
                                 path_to_mapped_genes = "/mnt/jack-5/amismailov/miRNA_study/GTF/gene_mapping.csv",
                                 min_sum_log = 100,
                                 scale_factor = 1000,
                                 remove_hla = TRUE) {
  # 1) Базовые фильтры и нормализация
  row_sums <- rowSums(counts)
  counts <- counts[row_sums > 0, , drop = FALSE]
  
  norm_counts <- t(t(counts) / colSums(counts)) * scale_factor
  log_counts <- log2(norm_counts + 1)
  
  gene_filter <- rowSums(log_counts) >= min_sum_log
  log_counts_filtered <- log_counts[gene_filter, , drop = FALSE]
  
  if (remove_hla) {
    hla_genes <- grep("^HLA", rownames(log_counts_filtered), value = TRUE)
    if (length(hla_genes) > 0) {
      log_counts_filtered <- log_counts_filtered[setdiff(rownames(log_counts_filtered), hla_genes), , drop = FALSE]
    }
  }
  
  # 2) Загрузка обновлённого маппинга (CSV)
  gene_chromosome_map <- read.csv(path_to_mapped_genes, stringsAsFactors = FALSE)
  
  # Приводим к ожидаемым именам
  col_symbol <- "hgnc_symbol"
  col_ens <- "ensembl_id"
  col_chr <- "chromosome"
  
  gene_chromosome_map[[col_ens]] <- as.character(gene_chromosome_map[[col_ens]])
  gene_chromosome_map[[col_symbol]] <- as.character(gene_chromosome_map[[col_symbol]])
  gene_chromosome_map[[col_chr]] <- as.character(gene_chromosome_map[[col_chr]])
  
  # Снимаем префикс "chr"
  gene_chromosome_map[[col_chr]] <- sub("^chr", "", gene_chromosome_map[[col_chr]])
  
  # 3) HGNC → ENS
  hgnc_in_matrix <- rownames(log_counts_filtered)
  idx <- match(hgnc_in_matrix, gene_chromosome_map[[col_symbol]])
  ens_mapped <- gene_chromosome_map[[col_ens]][idx]
  
  keep <- !is.na(ens_mapped) & nzchar(ens_mapped)
  log_counts_filtered <- log_counts_filtered[keep, , drop = FALSE]
  ens_mapped <- ens_mapped[keep]
  
  rownames(log_counts_filtered) <- ens_mapped
  
  # Убираем дубликаты ENS
  if (any(duplicated(ens_mapped))) {
    ens_levels <- unique(ens_mapped)
    agg_list <- lapply(ens_levels, function(e) {
      m <- log_counts_filtered[rownames(log_counts_filtered) == e, , drop = FALSE]
      if (nrow(m) == 1) m else matrix(colMeans(m), nrow = 1, dimnames = list(e, colnames(m)))
    })
    log_counts_filtered <- do.call(rbind, agg_list)
  }
  
  # 4) ENS → Chromosome
  map_sub <- gene_chromosome_map[gene_chromosome_map[[col_ens]] %in% rownames(log_counts_filtered),
                                 c(col_ens, col_chr), drop = FALSE]
  
  ens2chr <- unique(map_sub)
  
  # оставляем только аутосомы
  autosomes <- as.character(1:22)
  ens2chr <- ens2chr[ens2chr[[col_chr]] %in% autosomes, , drop = FALSE]
  
  # 5) Средние по хромосомам
  chr_list <- setNames(vector("list", length(autosomes)), autosomes)
  for (chrom in autosomes) {
    ens_on_chr <- ens2chr[[col_ens]][ens2chr[[col_chr]] == chrom]
    ens_on_chr <- intersect(ens_on_chr, rownames(log_counts_filtered))
    if (length(ens_on_chr) > 0) {
      chr_list[[chrom]] <- colMeans(log_counts_filtered[ens_on_chr, , drop = FALSE])
    } else {
      chr_list[[chrom]] <- rep(NA_real_, ncol(log_counts_filtered))
      names(chr_list[[chrom]]) <- colnames(log_counts_filtered)
    }
  }
  
  chromosome_means <- do.call(cbind, chr_list)
  colnames(chromosome_means) <- paste0("Chr", autosomes)
  rownames(chromosome_means) <- colnames(log_counts_filtered)
  
  as.data.frame(chromosome_means, check.names = FALSE)
}



#-----------------CHECK------------------
# counts <- GetAssayData(seurat.object, assay = "RNA", layer = "counts")
# chromosome_means = get_chromosome_means(counts)
# 
# chromosome_means

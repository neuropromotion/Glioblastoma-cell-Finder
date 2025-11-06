library(Seurat) 

get_chromosome_means <- function(counts, path_to_mapped_genes='path_to/mapped_genes.txt'){
  cat('Start getting chromosome_means...')
  row_sums <- rowSums(counts)
  counts <- counts[row_sums > 0, ]
  norm_counts <- t(t(counts) / colSums(counts) * 1000)  
  log_counts <- log2(norm_counts + 1) 
  gene_filter <- rowSums(log_counts) >= 100
  log_counts_filtered <- log_counts[gene_filter, ]
  # Remove HLA-genes
  hla_genes <- grep("HLA", rownames(log_counts_filtered), value = TRUE)
  log_counts_filtered <- log_counts_filtered[!(rownames(log_counts_filtered) %in% hla_genes), ]
  gene_chromosome_map <- read.table(path_to_mapped_genes, sep = '\t', header = T)  # Файл с генами и их хромосомами
  chromosome_means <- sapply(unique(gene_chromosome_map$Chromosome.scaffold.name), function(chrom) {
    genes <- intersect(
      gene_chromosome_map$HGNC.symbol[gene_chromosome_map$Chromosome.scaffold.name == chrom],
      rownames(log_counts_filtered)  
    )
    
    if (length(genes) > 0) {
      colMeans(log_counts_filtered[genes, , drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, ncol(log_counts_filtered))  
    }
  })
  chromosome_means <- as.data.frame(chromosome_means)
  colnames(chromosome_means) <- paste0("Chr", rep(1:22))
  rownames(chromosome_means) <- rownames(t(log_counts_filtered))
  return(chromosome_means)
}


#-----------------CHECK------------------
# counts <- GetAssayData(seurat.object, assay = "RNA", layer = "counts")
# chromosome_means = get_chromosome_means(counts)
# 
# chromosome_means





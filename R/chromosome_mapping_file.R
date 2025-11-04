gtf <- "path_to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz"
out <- "gene_mapping.txt"
 
# Читать GTF и оставить только feature == "gene"
dt <- fread(cmd = paste("zcat", shQuote(gtf)),
            sep = "\t", header = FALSE, quote = "",
            col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))

dt <- dt[feature == "gene"]

# Извлечение атрибутов: gene_id (ENS), gene_name (HGNC)
extract_attr <- function(x, key) {
  m <- str_match(x, paste0(key, " \"([^\"]+)\";"))
  m[,2]
}

dt[, ENS := extract_attr(attribute, "gene_id")]
dt[, HGNC.symbol := extract_attr(attribute, "gene_name")]

# Привести хромосому, убрать префикс chr, но ничего не фильтровать
dt[, Chromosome.scaffold.name := sub("^chr", "", seqname, ignore.case = TRUE)]

# Оставить нужные поля, НЕ удалять строки с пустым HGNC
res <- dt[, .(ENS, HGNC.symbol, Chromosome.scaffold.name)]

# Удаляем только строки без ENS (некорректные)
res <- res[!(is.na(ENS) | ENS == "")]

# Сохранить как txt (tab-delimited)
fwrite(res, file = out, sep = "\t", quote = FALSE)

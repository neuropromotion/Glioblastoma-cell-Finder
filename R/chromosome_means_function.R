library(httr2)
library(jsonlite)
library(progress)

.extract_ensembl_gene <- function(ens) {
  if (is.null(ens)) return(NA_character_)
  if (is.list(ens)) {
    # ens может быть списком (list of lists) или именованным списком (dict)
    if (!is.null(ens$gene)) {
      return(as.character(ens$gene))
    } else if (length(ens) > 0 && !is.null(ens[[1]]$gene)) {
      return(as.character(ens[[1]]$gene))
    }
  }
  NA_character_
}

# Поиск по синонимам для одного символа
.search_by_synonyms <- function(symbol, verbose = FALSE) {
  # 1) получить alias/other_names и ensembl
  q <- request("https://mygene.info/v3/query") |>
    req_url_query(q = symbol,
                  species = "human",
                  fields = "symbol,alias,other_names,ensembl.gene",
                  size = 1) |>
    req_perform()
  
  res <- q |> resp_body_string() |> fromJSON(simplifyVector = FALSE)
  if (is.null(res$hits) || length(res$hits) == 0) return(NA_character_)
  hit <- res$hits[[1]]
  
  # если ensembl уже есть — вернуть его
  ens <- .extract_ensembl_gene(hit$ensembl)
  if (!is.na(ens)) return(ens)
  
  # собрать синонимы
  aliases <- character(0)
  if (!is.null(hit$alias)) {
    if (is.character(hit$alias)) aliases <- c(aliases, hit$alias)
  }
  if (!is.null(hit$other_names)) {
    if (is.character(hit$other_names)) aliases <- c(aliases, hit$other_names)
  }
  aliases <- unique(setdiff(aliases, symbol))
  if (length(aliases) == 0) return(NA_character_)
  if (verbose) message(sprintf("  %s → синонимы: %s",
                               symbol,
                               paste(utils::head(aliases, 5), collapse = ", ")))
  
  # 2) искать по каждому синониму точечным запросом
  for (syn in aliases) {
    qs <- request("https://mygene.info/v3/query") |>
      req_url_query(q = syn,
                    scopes = "symbol",
                    fields = "ensembl.gene",
                    species = "human",
                    size = 1) |>
      req_perform()
    
    rs <- qs |> resp_body_string() |> fromJSON(simplifyVector = FALSE)
    if (!is.null(rs$hits) && length(rs$hits) > 0) {
      ens2 <- .extract_ensembl_gene(rs$hits[[1]]$ensembl)
      if (!is.na(ens2)) return(ens2)
    }
  }
  NA_character_
}

# Основная функция: вектор символов -> именованный вектор Ensembl ID
find_ens <- function(list_of_genes, use_synonyms = TRUE, verbose = FALSE, batch_size = 1000) {
  stopifnot(is.character(list_of_genes))
  genes <- unique(list_of_genes)
  n <- length(genes)
  
  # 1) batch querymany: scopes=symbol, fields=ensembl.gene
  # будем отправлять пачками, чтобы не упереться в лимиты URL
  results <- setNames(rep(NA_character_, n), genes)
  
  if (n > 0) {
    pb <- progress_bar$new(format = "Primary search [:bar] :percent :current/:total",
                           total = ceiling(n / batch_size), clear = FALSE, width = 70)
    
    for (i in seq(1, n, by = batch_size)) {
      idx <- i:min(i + batch_size - 1, n)
      payload <- list(
        q = unname(genes[idx]),
        scopes = "symbol",
        fields = "ensembl.gene",
        species = "human"
      )
      req <- request("https://mygene.info/v3/query") |>
        req_method("POST") |>
        req_headers("Content-Type" = "application/json") |>
        req_body_json(payload)
      
      resp <- req_perform(req)
      dat  <- resp |> resp_body_string() |> fromJSON(simplifyVector = FALSE)
      
      # dat — список ответов по каждому запросу
      if (is.list(dat) && length(dat) == length(idx)) {
        for (k in seq_along(idx)) {
          entry <- dat[[k]]
          sym <- entry$query %||% genes[idx[k]]
          ens <- .extract_ensembl_gene(entry$ensembl)
          nf  <- isTRUE(entry$notfound)
          if (!is.na(ens) && !nf) {
            results[sym] <- ens
          } else {
            # пока оставляем NA — попадёт в фазу синонимов
          }
        }
      }
      pb$tick()
    }
  }
  
  # 2) поиск по синонимам для оставшихся
  if (use_synonyms) {
    missing <- names(results)[is.na(results)]
    if (length(missing) > 0) {
      if (verbose) message(sprintf("\n🔍 Поиск по синонимам для %d генов…", length(missing)))
      pb2 <- progress_bar$new(format = "Synonym search  [:bar] :percent :current/:total",
                              total = length(missing), clear = FALSE, width = 70)
      for (g in missing) {
        ans <- tryCatch(.search_by_synonyms(g, verbose = verbose),
                        error = function(e) NA_character_)
        if (!is.na(ans)) results[g] <- ans
        pb2$tick()
      }
    }
  }
  
  # 3) отчёт
  found <- sum(!is.na(results))
  total <- length(results)
  message(sprintf("\n📊 Результаты:\n  Найдено: %d/%d (%.1f%%)\n  Не найдено: %d",
                  found, total, 100 * found / total, total - found))
  
  # вернуть только найденные (без NA)
  results_clean <- results[!is.na(results)]
  return(results_clean)
}
get_chromosome_means <- function(counts,
                                 path_to_mapped_genes = "path_to/gene_mapping.txt",
                                 min_sum_log = 100,
                                 scale_factor = 1000,
                                 remove_hla = TRUE) {
  # 1) Базовые фильтры и нормализация
  row_sums <- rowSums(counts)
  counts <- counts[row_sums > 0, , drop = FALSE]
  
  # CPM-like с кастомным scale_factor и лог
  norm_counts <- t(t(counts) / colSums(counts)) * scale_factor
  log_counts <- log2(norm_counts + 1)
  
  # Фильтр по сумме лог-экспрессий
  gene_filter <- rowSums(log_counts) >= min_sum_log
  log_counts_filtered <- log_counts[gene_filter, , drop = FALSE]
  
  # Опционально убрать HLA
  if (remove_hla) {
    hla_genes <- grep("^HLA", rownames(log_counts_filtered), value = TRUE)
    if (length(hla_genes) > 0) {
      log_counts_filtered <- log_counts_filtered[setdiff(rownames(log_counts_filtered), hla_genes), , drop = FALSE]
    }
  }
  
  # 2) Загрузка маппинга генов
  # Ожидаемые колонки: HGNC.symbol, ENS, Chromosome.scaffold.name
  gene_chromosome_map <- read.table(path_to_mapped_genes, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Приведём имена колонок к стандартным (если вдруг они чуть отличаются)
  nm <- names(gene_chromosome_map)
  # Попробуем автоматически найти подходящие колонки
  col_symbol <- nm[grepl("^hgnc(\\.symbol)?$|^symbol$", tolower(nm))][1]
  if (is.na(col_symbol)) col_symbol <- "HGNC.symbol"
  col_ens <- nm[grepl("^ens(\\.gene(id)?)?$|ensembl", tolower(nm))][1]
  if (is.na(col_ens)) col_ens <- "ENS"
  col_chr <- nm[grepl("chromosome|scaffold", tolower(nm))][1]
  if (is.na(col_chr)) col_chr <- "Chromosome.scaffold.name"
  
  stopifnot(all(c(col_symbol, col_ens, col_chr) %in% names(gene_chromosome_map)))
  
  # Удалим пустые/NA ENS, приведём к character
  gene_chromosome_map[[col_symbol]] <- as.character(gene_chromosome_map[[col_symbol]])
  gene_chromosome_map[[col_ens]]    <- as.character(gene_chromosome_map[[col_ens]])
  gene_chromosome_map[[col_chr]]    <- as.character(gene_chromosome_map[[col_chr]])
  gene_chromosome_map <- gene_chromosome_map[!is.na(gene_chromosome_map[[col_ens]]) & nzchar(gene_chromosome_map[[col_ens]]), , drop = FALSE]
  
  # 3) Сопоставление HGNC -> ENS для строк log_counts_filtered
  hgnc_in_matrix <- rownames(log_counts_filtered)
  # Берём первый матч по символу (если дубликаты — выбираем первый)
  idx <- match(hgnc_in_matrix, gene_chromosome_map[[col_symbol]])
  ens_mapped <- gene_chromosome_map[[col_ens]][idx]
  
  # Есть случаи, когда символ не нашёлся — выкинем такие строки
  keep <- !is.na(ens_mapped) & nzchar(ens_mapped)
  log_counts_filtered <- log_counts_filtered[keep, , drop = FALSE]
  ens_mapped <- ens_mapped[keep]
  
  # Если остались дубликаты ENS (разные символы к одному ENS) — агрегируем по ENS как среднее
  # Сначала заменим rownames на ENS
  rownames(log_counts_filtered) <- ens_mapped
  
  # Аггрегация по дубликатам ENS (среднее по строкам с одним ENS)
  # Быстро через tapply-индексацию
  if (any(duplicated(ens_mapped))) {
    ens_levels <- unique(ens_mapped)
    # усреднить по каждой группе ENS
    agg_list <- lapply(ens_levels, function(e) {
      m <- log_counts_filtered[rownames(log_counts_filtered) == e, , drop = FALSE]
      if (nrow(m) == 1) {
        m
      } else {
        matrix(colMeans(m, na.rm = TRUE), nrow = 1, dimnames = list(e, colnames(m)))
      }
    })
    log_counts_filtered <- do.call(rbind, agg_list)
  }
  
  # 4) Подготовим множество хромосом и ENS->CHR словарь
  # Оставим только ENS, которые реально присутствуют в матрице
  map_sub <- gene_chromosome_map[gene_chromosome_map[[col_ens]] %in% rownames(log_counts_filtered), c(col_ens, col_chr), drop = FALSE]
  
  # На случай дубликатов ENS->CHR (бывает редкая неоднозначность): берём наиболее часто встречающуюся хромосому
  if (any(duplicated(map_sub[[col_ens]]))) {
    map_sub <- aggregate(setNames(rep(1, nrow(map_sub)), NULL),
                         by = list(ENS = map_sub[[col_ens]], CHR = map_sub[[col_chr]]), FUN = sum)
    # для каждого ENS выбрать CHR с макс. счётом
    best_chr <- tapply(map_sub$x, map_sub$ENS, function(v) {
      chrs <- map_sub$CHR[map_sub$ENS == names(v)[1]]
      chrs[which.max(v)]
    })
    ens2chr <- data.frame(ENS = names(best_chr), CHR = unname(best_chr), stringsAsFactors = FALSE)
  } else {
    ens2chr <- data.frame(ENS = map_sub[[col_ens]], CHR = map_sub[[col_chr]], stringsAsFactors = FALSE)
  }
  
  # Оставим только аутосомы 1..22
  autosomes <- as.character(1:22)
  ens2chr <- ens2chr[ens2chr$CHR %in% autosomes, , drop = FALSE]
  
  # 5) Усреднение по аутосомам
  # Для каждой хромосомы берём ENS из ens2chr и усредняем соответствующие строки
  chr_list <- setNames(vector("list", length(autosomes)), autosomes)
  for (chrom in autosomes) {
    ens_on_chr <- ens2chr$ENS[ens2chr$CHR == chrom]
    ens_on_chr <- intersect(ens_on_chr, rownames(log_counts_filtered))
    if (length(ens_on_chr) > 0) {
      chr_list[[chrom]] <- colMeans(log_counts_filtered[ens_on_chr, , drop = FALSE], na.rm = TRUE)
    } else {
      chr_list[[chrom]] <- rep(NA_real_, ncol(log_counts_filtered))
      names(chr_list[[chrom]]) <- colnames(log_counts_filtered)
    }
  }
  
  chromosome_means <- do.call(cbind, chr_list)
  colnames(chromosome_means) <- paste0("Chr", autosomes)
  rownames(chromosome_means) <- colnames(log_counts_filtered)  # клетки по строкам
  
  # Вернём data.frame
  chromosome_means <- as.data.frame(chromosome_means, check.names = FALSE)
  return(chromosome_means)
}

#-----------------CHECK------------------
# genes <- c("ELPMA", "DCX", "ADGRG1", "AGT")  # FOO123 — заведомо фальшивый
# res <- find_ens(genes, use_synonyms = TRUE, verbose = TRUE)
# unname(res)
# 
# counts <- GetAssayData(merged_filtered, assay = "RNA", layer = "counts")
# chromosome_means = get_chromosome_means(counts)
# 




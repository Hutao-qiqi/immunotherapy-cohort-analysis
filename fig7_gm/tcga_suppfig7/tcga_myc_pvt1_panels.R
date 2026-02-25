suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(hexbin)
  library(patchwork)
  library(cowplot)
  library(httr2)
  library(jsonlite)
  library(digest)
  library(pROC)
  library(ggpubr)
  library(ggrepel)
})

canonical_tcga <- function(x) {
  x <- as.character(x)
  m <- regexpr("TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  hit <- m != -1
  out[hit] <- substr(x[hit], m[hit], m[hit] + attr(m, "match.length")[hit] - 1)
  out
}

tcga_sample_type <- function(sample_id) {
  sample_id <- as.character(sample_id)
  m <- regexec("TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-([0-9]{2})$", sample_id, perl = TRUE)
  regm <- regmatches(sample_id, m)
  out <- rep(NA_character_, length(sample_id))
  hit <- lengths(regm) == 2
  out[hit] <- vapply(regm[hit], function(z) z[[2]], character(1))
  out
}

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(argv)) {
      out[[substring(key, 3)]] <- TRUE
      break
    }
    val <- argv[[i + 1]]
    if (startsWith(val, "--")) {
      out[[substring(key, 3)]] <- TRUE
      i <- i + 1
      next
    }
    out[[substring(key, 3)]] <- val
    i <- i + 2
  }
  out
}

fmt_p <- function(p) {
  if (is.na(p)) return("P=NA")
  if (p < 1e-4) return("P<0.0001")
  sprintf("P=%.3g", p)
}

safe_ggsave <- function(path, plot, ..., alt_suffix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    {
      ggplot2::ggsave(filename = path, plot = plot, ...)
      path
    },
    error = function(e) {
      alt <- sub("\\.pdf$", paste0("_", alt_suffix, ".pdf"), path, ignore.case = TRUE)
      message("ggsave failed for ", path, " (", conditionMessage(e), "). Writing to: ", alt)
      ggplot2::ggsave(filename = alt, plot = plot, ...)
      alt
    }
  )
}

ppv_metrics_from_scores <- function(y, score, direction = "<") {
  ok <- !is.na(y) & !is.na(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(y) < 10 || length(unique(y)) < 2) return(NULL)

  roc_obj <- pROC::roc(y, score, quiet = TRUE, direction = direction)
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    transpose = FALSE,
    ret = c("threshold", "sensitivity", "specificity")
  )

  thr <- as.numeric(best$threshold[[1]])
  sens <- as.numeric(best$sensitivity[[1]])
  spec <- as.numeric(best$specificity[[1]])
  auc <- as.numeric(pROC::auc(roc_obj))

  pred_pos <- if (direction == "<") score >= thr else score <= thr
  tp <- sum(pred_pos & y == 1L)
  fp <- sum(pred_pos & y == 0L)
  tn <- sum((!pred_pos) & y == 0L)
  fn <- sum((!pred_pos) & y == 1L)

  ppv <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  npv <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
  acc <- (tp + tn) / (tp + fp + tn + fn)
  prev <- mean(y == 1L)

  list(
    threshold = thr,
    auc = auc,
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = npv,
    accuracy = acc,
    prevalence = prev,
    tp = tp, fp = fp, tn = tn, fn = fn,
    n = length(y)
  )
}

plot_roc_from_y_long <- function(dt_long, title, direction = "<") {
  dt <- copy(dt_long)
  dt <- dt[!is.na(y) & !is.na(expr)]

  rocs <- dt[, {
    ok <- !is.na(expr) & !is.na(y)
    if (sum(ok) < 10 || length(unique(y[ok])) < 2) return(list(roc = NULL, auc = NA_real_))
    r <- pROC::roc(y[ok], expr[ok], quiet = TRUE, direction = direction)
    list(roc = list(r), auc = as.numeric(pROC::auc(r)))
  }, by = gene]

  curves <- rbindlist(lapply(seq_len(nrow(rocs)), function(i) {
    g <- rocs$gene[[i]]
    r_item <- rocs$roc[[i]]
    if (is.null(r_item)) return(NULL)
    if (inherits(r_item, "roc")) {
      r <- r_item
    } else if (is.list(r_item) && length(r_item) >= 1 && inherits(r_item[[1]], "roc")) {
      r <- r_item[[1]]
    } else {
      return(NULL)
    }
    data.table(gene = g, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }), fill = TRUE)

  auc_labels <- rocs[, .(gene, auc)]
  auc_labels[, label := sprintf("%s AUC=%.3f", gene, auc)]

  ggplot(curves, aes(x = fpr, y = tpr, color = gene)) +
    geom_line(linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    coord_equal() +
    labs(x = "False positive rate", y = "True positive rate", title = title) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    guides(color = guide_legend(title = NULL)) +
    annotate(
      "text",
      x = 0.58,
      y = 0.18,
      label = paste(auc_labels$label, collapse = "\n"),
      hjust = 0,
      size = 3.3
    )
}

gdc_map_aliquot_to_sample <- function(aliquot_ids, cache_path = NULL, batch_size = 50, sleep_s = 0.15) {
  aliquot_ids <- unique(as.character(aliquot_ids))
  aliquot_ids <- aliquot_ids[!is.na(aliquot_ids) & aliquot_ids != ""]

  cached <- data.table(aliquot_id = character(), sample_submitter_id = character())
  missing <- aliquot_ids
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached0 <- fread(cache_path)
    if (all(c("aliquot_id", "sample_submitter_id") %in% names(cached0))) {
      cached <- unique(cached0[, .(aliquot_id = as.character(aliquot_id), sample_submitter_id = as.character(sample_submitter_id))])
      missing <- setdiff(aliquot_ids, cached$aliquot_id)
    }
  }
  if (length(missing) == 0) return(cached)

  uuid_re <- "^[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}$"
  missing <- missing[grepl(uuid_re, missing, perl = TRUE)]
  if (length(missing) == 0) return(cached)

  map_rows <- list()
  for (i in seq(1, length(missing), by = batch_size)) {
    batch <- missing[i:min(i + batch_size - 1, length(missing))]
    filt <- jsonlite::toJSON(
      list(op = "in", content = list(field = "samples.portions.analytes.aliquots.aliquot_id", value = as.list(batch))),
      auto_unbox = TRUE
    )
    req <- request("https://api.gdc.cancer.gov/cases") |>
      req_url_query(
        filters = filt,
        fields = "samples.submitter_id,samples.portions.analytes.aliquots.aliquot_id",
        size = "2000"
      )
    resp <- req_perform(req)
    body <- resp_body_json(resp, simplifyVector = FALSE)
    hits <- body$data$hits
    if (length(hits) > 0) {
      for (h in hits) {
        samples <- h$samples
        if (length(samples) == 0) next
        for (s in samples) {
          sid <- s$submitter_id
          portions <- s$portions
          if (length(portions) == 0) next
          for (p in portions) {
            analytes <- p$analytes
            if (length(analytes) == 0) next
            for (a in analytes) {
              aliquots <- a$aliquots
              if (length(aliquots) == 0) next
              for (q in aliquots) {
                aid <- q$aliquot_id
                if (is.null(aid) || is.null(sid)) next
                map_rows[[length(map_rows) + 1]] <- list(aliquot_id = as.character(aid), sample_submitter_id = as.character(sid))
              }
            }
          }
        }
      }
    }
    Sys.sleep(sleep_s)
  }

  fresh <- if (length(map_rows) == 0) {
    data.table(aliquot_id = character(), sample_submitter_id = character())
  } else {
    unique(rbindlist(map_rows))
  }
  out <- unique(rbindlist(list(cached, fresh), fill = TRUE))
  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    fwrite(out, cache_path, sep = "\t")
  }
  out
}

read_gene_coords <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("genes", envir = e, inherits = FALSE)) {
    stop("Gene coords RDA missing `genes`: ", path)
  }
  genes <- as.data.table(get("genes", envir = e))
  setnames(genes, c("GeneSymbol", "Chr", "Start", "End"))
  genes[, GeneSymbol := as.character(GeneSymbol)]
  genes[, Chr := as.integer(Chr)]
  genes[, Start := as.integer(Start)]
  genes[, End := as.integer(End)]
  genes[, Pos := as.integer(floor((Start + End) / 2))]
  genes
}

discover_gistic_files <- function(cnv_root) {
  all_files <- list.files(cnv_root, pattern = "all_data_by_genes\\.txt$", recursive = TRUE, full.names = TRUE)
  thr_files <- list.files(cnv_root, pattern = "all_thresholded\\.by_genes\\.txt$", recursive = TRUE, full.names = TRUE)

  extract_code <- function(p) {
    parts <- strsplit(normalizePath(p, winslash = "\\", mustWork = FALSE), "\\\\")[[1]]
    cnv_dir <- NA_character_
    for (seg in parts) {
      if (grepl("cnv", seg, ignore.case = TRUE)) {
        cnv_dir <- seg
        break
      }
    }
    toupper(sub("^([A-Za-z0-9]+).*$", "\\1", cnv_dir))
  }

  best_by_code <- function(files) {
    if (length(files) == 0) return(data.table(cancer = character(), path = character(), ncol = integer()))
    dt <- rbindlist(lapply(files, function(p) {
      hdr <- tryCatch(readLines(p, n = 1, warn = FALSE), error = function(e) NA_character_)
      ncol <- if (length(hdr) == 1 && !is.na(hdr)) length(strsplit(hdr, "\t", fixed = TRUE)[[1]]) else NA_integer_
      data.table(cancer = extract_code(p), path = p, ncol = as.integer(ncol))
    }), fill = TRUE)
    dt <- dt[!is.na(cancer) & cancer != "" & !is.na(ncol)]
    setorder(dt, cancer, -ncol, path)
    dt[, .SD[1], by = cancer]
  }

  list(
    all = best_by_code(all_files),
    thr = best_by_code(thr_files)
  )
}

read_gistic_gene_vector <- function(path, gene, cache_gdc_map = NULL) {
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) != 1) stop("Failed to read header: ", path)
  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (length(cols) < 5) stop("Unexpected GISTIC format: ", path)
  sample_cols <- cols[4:length(cols)]
  sample <- canonical_tcga(sample_cols)

  uuid_re <- "^[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}$"
  need_map <- is.na(sample) & grepl(uuid_re, sample_cols, perl = TRUE)
  if (any(need_map)) {
    map_dt <- gdc_map_aliquot_to_sample(sample_cols[need_map], cache_path = cache_gdc_map)
    if (nrow(map_dt) > 0) {
      setkey(map_dt, aliquot_id)
      mapped <- map_dt[J(sample_cols[need_map]), sample_submitter_id]
      sample[need_map] <- canonical_tcga(mapped)
    }
  }

  pref <- paste0(gene, "\t")
  found <- NULL
  repeat {
    chunk <- readLines(con, n = 20000, warn = FALSE)
    if (length(chunk) == 0) break
    idx <- which(startsWith(chunk, pref))
    if (length(idx) > 0) {
      found <- chunk[[idx[[1]]]]
      break
    }
  }
  if (is.null(found)) return(data.table(sample = character(), value = numeric()))
  fields <- strsplit(found, "\t", fixed = TRUE)[[1]]
  vals <- suppressWarnings(as.numeric(fields[4:length(fields)]))
  dt <- data.table(sample = sample, value = vals)
  dt <- dt[!is.na(sample) & !is.na(value)]
  dt <- dt[, .(value = mean(value)), by = sample]
  dt
}

extract_cnv_from_segments <- function(seg_dt, chr, pos, sample_col, chr_col, start_col, end_col, value_col) {
  dt <- copy(as.data.table(seg_dt))
  setnames(dt, c(sample_col, chr_col, start_col, end_col, value_col), c("Sample", "Chr", "Start", "End", "Value"), skip_absent = TRUE)
  dt[, Sample := canonical_tcga(Sample)]
  dt <- dt[!is.na(Sample)]
  dt <- dt[as.integer(Chr) == as.integer(chr) & as.integer(Start) <= as.integer(pos) & as.integer(End) >= as.integer(pos)]
  if (nrow(dt) == 0) return(data.table(sample = character(), value = numeric()))
  dt[, .(value = mean(as.numeric(Value), na.rm = TRUE)), by = Sample][, .(sample = Sample, value)]
}

xena_post <- function(host, query) {
  r <- request(paste0(host, "/data/")) |>
    req_body_raw(query) |>
    req_headers("Content-Type" = "text/plain") |>
    req_perform()
  resp_body_string(r)
}

xena_dataset_probe_values <- function(host, dataset, samples, probes, query_path = file.path("xena_queries", "datasetProbeValues.xq")) {
  qfn <- paste(readLines(query_path, warn = FALSE), collapse = "\n")
  q <- paste0("(", qfn, " ",
              jsonlite::toJSON(dataset, auto_unbox = TRUE), " ",
              jsonlite::toJSON(as.list(samples), auto_unbox = TRUE), " ",
              jsonlite::toJSON(as.list(probes), auto_unbox = TRUE),
              ")")
  txt <- xena_post(host, q)
  x <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  values <- x[[2]]
  dt <- data.table(sample = samples)
  for (i in seq_along(probes)) {
    v <- suppressWarnings(as.numeric(values[[i]]))
    dt[[probes[[i]]]] <- v
  }
  dt
}

xena_expr_fetch_cached <- function(host, dataset, samples, probes, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key_hash <- digest::digest(list(dataset = dataset, samples = samples, probes = probes), algo = "xxhash64")
  key <- paste0("xena_expr_", gsub("[^A-Za-z0-9]+", "_", dataset), "_", key_hash, ".tsv")
  cache_path <- file.path(cache_dir, key)
  if (file.exists(cache_path)) {
    dt <- fread(cache_path)
    if (all(c("sample", probes) %in% names(dt))) return(dt[, c("sample", probes), with = FALSE])
  }
  dt <- xena_dataset_probe_values(host, dataset, samples, probes)
  fwrite(dt, cache_path, sep = "\t")
  dt
}

chunked <- function(x, size) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / size))
}

xena_expr_fetch_cached_batched <- function(host, dataset, samples, probes, cache_dir, batch_size = 50) {
  parts <- chunked(probes, batch_size)
  if (length(parts) == 0) return(data.table(sample = samples))
  out <- data.table(sample = samples)
  for (p in parts) {
    dt <- xena_expr_fetch_cached(host, dataset, samples, p, cache_dir = cache_dir)
    out <- merge(out, dt, by = "sample", all.x = TRUE, all.y = FALSE)
  }
  out
}

plot_panel_a_hex <- function(df_long, method) {
  stats <- df_long[, {
    ct <- suppressWarnings(cor.test(cnv, expr, method = method, exact = FALSE))
    .(rho = unname(ct$estimate[[1]]), p = ct$p.value, n = .N)
  }, by = gene]
  stats[, label := sprintf("Rho=%.2f, %s, n=%d", rho, vapply(p, fmt_p, character(1)), n)]

  ggplot(df_long, aes(x = cnv, y = expr)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(name = "Count", trans = "log10") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, color = "white") +
    facet_wrap(~gene, scales = "free_y") +
    geom_text(
      data = stats,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.15, size = 3.2
    ) +
    labs(x = "CNV (segment mean)", y = "mRNA expression (log2 scale)", title = "Panel A: Continuous CNV–mRNA correlation") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
}

plot_panel_b_violin <- function(df_cat) {
  df_cat[, gistic := factor(gistic, levels = c("Loss", "Diploid", "Gain", "Amp"))]
  ggplot(df_cat, aes(x = gistic, y = expr, fill = gistic)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.2) +
    geom_boxplot(width = 0.14, outlier.size = 0.3, alpha = 0.6) +
    facet_wrap(~gene, scales = "free_y") +
    ggpubr::stat_compare_means(
      comparisons = list(c("Diploid", "Amp")),
      method = "wilcox.test",
      label = "p.signif",
      size = 4
    ) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = "GISTIC state", y = "mRNA expression", title = "Panel B: Categorical validation (Diploid vs Amp)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
}

plot_panel_d_bubble <- function(df_cor) {
  df <- copy(df_cor)
  df[, neglog10p := -log10(pmax(p, 1e-300))]
  order_dt <- df[, .(order_score = mean(r, na.rm = TRUE)), by = cancer][order(-order_score)]
  df[, cancer := factor(cancer, levels = order_dt$cancer)]

  ggplot(df, aes(x = gene, y = cancer)) +
    geom_point(aes(size = neglog10p, color = r), alpha = 0.95) +
    scale_color_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, name = "Spearman R") +
    scale_size_continuous(name = expression(-log[10](p))) +
    labs(x = NULL, y = "TCGA cancer type", title = "Panel D: Pan-cancer consistency") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_panel_e_roc <- function(df_cat) {
  dt <- df_cat[gistic %in% c("Loss", "Diploid", "Gain", "Amp")]
  dt[, y := fifelse(gistic %in% c("Gain", "Amp"), 1L, 0L)]

  plot_roc_from_y_long(
    dt_long = dt[, .(sample, gene, expr, y)],
    title = "Panel E: Predictive power (mRNA -> CNV Gain/Amp)",
    direction = "<"
  )
}

extract_cytoband_genes <- function(path, cytoband_prefix = "8q24") {
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  readLines(con, n = 1, warn = FALSE)
  genes <- character()
  repeat {
    chunk <- readLines(con, n = 20000, warn = FALSE)
    if (length(chunk) == 0) break
    for (ln in chunk) {
      tab1 <- regexpr("\t", ln, fixed = TRUE)[[1]]
      if (tab1 <= 1) next
      gene <- substr(ln, 1, tab1 - 1)
      rest <- substr(ln, tab1 + 1, nchar(ln))
      tab2 <- regexpr("\t", rest, fixed = TRUE)[[1]]
      if (tab2 <= 1) next
      rest2 <- substr(rest, tab2 + 1, nchar(rest))
      tab3 <- regexpr("\t", rest2, fixed = TRUE)[[1]]
      if (tab3 <= 1) next
      cyt <- substr(rest2, 1, tab3 - 1)
      if (startsWith(cyt, cytoband_prefix)) genes <- c(genes, gene)
    }
  }
  unique(genes)
}

read_gistic_gene_vectors <- function(path, genes, cache_gdc_map = NULL) {
  genes <- unique(as.character(genes))
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1, warn = FALSE)
  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  sample_cols <- cols[4:length(cols)]
  sample <- canonical_tcga(sample_cols)

  uuid_re <- "^[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}$"
  need_map <- is.na(sample) & grepl(uuid_re, sample_cols, perl = TRUE)
  if (any(need_map)) {
    map_dt <- gdc_map_aliquot_to_sample(sample_cols[need_map], cache_path = cache_gdc_map)
    if (nrow(map_dt) > 0) {
      setkey(map_dt, aliquot_id)
      mapped <- map_dt[J(sample_cols[need_map]), sample_submitter_id]
      sample[need_map] <- canonical_tcga(mapped)
    }
  }

  wanted <- setNames(as.list(rep(NA_character_, length(genes))), genes)
  found <- setNames(rep(FALSE, length(genes)), genes)
  gene_prefixes <- paste0(genes, "\t")

  repeat {
    chunk <- readLines(con, n = 25000, warn = FALSE)
    if (length(chunk) == 0) break
    for (j in seq_along(genes)) {
      g <- genes[[j]]
      if (found[[g]]) next
      pref <- gene_prefixes[[j]]
      idx <- which(startsWith(chunk, pref))
      if (length(idx) > 0) {
        wanted[[g]] <- chunk[[idx[[1]]]]
        found[[g]] <- TRUE
      }
    }
    if (all(found)) break
  }

  out <- list()
  for (g in genes) {
    line <- wanted[[g]]
    if (is.na(line)) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    vals <- suppressWarnings(as.numeric(fields[4:length(fields)]))
    dt <- data.table(sample = sample, value = vals)
    dt <- dt[!is.na(sample) & !is.na(value)]
    dt <- dt[, .(value = mean(value)), by = sample]
    out[[g]] <- dt
  }
  out
}

plot_panel_c_rank <- function(df_rank) {
  df <- copy(df_rank)
  df[, is_key := gene %in% c("MYC", "PVT1")]
  df[, gene_key := paste0(gene, "___", cancer)]
  lvl <- df[, .(gene_key = gene_key[order(r)]), by = cancer]$gene_key
  df[, gene_key := factor(gene_key, levels = rev(lvl))]

  key_df <- df[is_key == TRUE]

  ggplot(df, aes(x = r, y = gene_key)) +
    geom_segment(aes(x = 0, xend = r, yend = gene_key, color = is_key), linewidth = 0.55, alpha = 0.9) +
    geom_point(aes(color = is_key), size = ifelse(df$is_key, 2.6, 1.8), alpha = 0.95) +
    ggrepel::geom_label_repel(
      data = key_df,
      aes(label = gene),
      color = "#D7301F",
      fill = "white",
      label.size = 0.25,
      size = 3.6,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "#D7301F",
      segment.size = 0.5,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "grey55"), guide = "none") +
    facet_wrap(~cancer, scales = "free_y") +
    labs(x = "Spearman R (CNV vs mRNA)", y = NULL, title = "Panel C: 8q24 specificity (SKCM & BLCA)") +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_panel_c_schematic <- function(genes, gene_coords, title = "8q24 locus") {
  coords <- gene_coords[GeneSymbol %in% genes & Chr == 8]
  if (nrow(coords) == 0) coords <- gene_coords[GeneSymbol %in% c("PVT1", "MYC") & Chr == 8]
  coords <- coords[!is.na(Start) & !is.na(End)]
  coords[, mid := as.numeric((Start + End) / 2)]
  coords <- coords[order(mid)]
  coords[, is_key := GeneSymbol %in% c("PVT1", "MYC")]

  if (nrow(coords) == 0) {
    return(ggplot() + theme_void() + labs(title = title))
  }

  pad <- 5e5
  y_min <- min(coords$Start) - pad
  y_max <- max(coords$End) + pad

  ticks <- coords[, .(gene = GeneSymbol, y = mid, is_key)]

  ggplot() +
    geom_segment(aes(x = 0.5, xend = 0.5, y = y_min, yend = y_max), linewidth = 6, color = "grey85") +
    geom_segment(
      data = ticks,
      aes(x = 0.30, xend = 0.70, y = y, yend = y, color = is_key),
      linewidth = ifelse(ticks$is_key, 1.15, 0.35),
      alpha = 0.95,
      lineend = "round"
    ) +
    geom_text(
      data = ticks[is_key == TRUE],
      aes(x = 0.74, y = y, label = gene),
      hjust = 0,
      color = "#D7301F",
      size = 3.6,
      fontface = "bold"
    ) +
    scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "black"), guide = "none") +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    xlim(0, 1) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 6))
    )
}

plot_panel_c_combined <- function(df_rank, gene_coords) {
  genes <- unique(df_rank$gene)
  p_left <- plot_panel_c_schematic(genes, gene_coords, title = "8q24")
  p_right <- plot_panel_c_rank(df_rank)
  p_left + p_right + plot_layout(widths = c(1.05, 3.0))
}

summarize_panelc_pancancer <- function(per_cancer_dt, method = "median") {
  dt <- copy(per_cancer_dt)
  dt <- dt[!is.na(r) & !is.na(n) & n >= 10]
  if (nrow(dt) == 0) return(data.table(gene = character(), r = numeric(), cancers = integer()))
  if (method == "median") {
    dt[, .(r = median(r, na.rm = TRUE), cancers = .N), by = gene]
  } else if (method == "mean_fisherz") {
    dt[, {
      z <- atanh(pmin(pmax(r, -0.999999), 0.999999))
      w <- sqrt(n - 3)
      zbar <- sum(z * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
      .(r = tanh(zbar), cancers = .N)
    }, by = gene]
  } else {
    stop("Unknown summary method: ", method)
  }
}

plot_panel_c_pancancer_rank <- function(df_sum) {
  df <- copy(df_sum)
  df <- df[order(r)]
  df[, is_key := gene %in% c("MYC", "PVT1")]
  df[, gene_key := paste0(gene, "___", seq_len(.N))]
  df[, gene_key := factor(gene_key, levels = df$gene_key)]
  key_df <- df[is_key == TRUE]

  ggplot(df, aes(x = r, y = gene_key)) +
    geom_segment(aes(x = 0, xend = r, yend = gene_key, color = is_key), linewidth = 0.55, alpha = 0.9) +
    geom_point(aes(color = is_key), size = ifelse(df$is_key, 2.6, 1.7), alpha = 0.95) +
    ggrepel::geom_label_repel(
      data = key_df,
      aes(label = gene),
      color = "#D7301F",
      fill = "white",
      label.size = 0.25,
      size = 3.8,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "#D7301F",
      segment.size = 0.5,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "grey55"), guide = "none") +
    labs(x = "Pan-cancer Spearman R (CNV vs mRNA)", y = NULL, title = "Panel C (pan-cancer): 8q24 specificity") +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_panel_c_pancancer_combined <- function(df_sum, gene_coords) {
  genes <- unique(df_sum$gene)
  p_left <- plot_panel_c_schematic(genes, gene_coords, title = "8q24")
  p_right <- plot_panel_c_pancancer_rank(df_sum)
  p_left + p_right + plot_layout(widths = c(1.05, 3.0))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  cnv_root <- if (!is.null(args[["cnv-root"]])) args[["cnv-root"]] else "."
  outdir <- if (!is.null(args[["outdir"]])) args[["outdir"]] else file.path("results", "tcga_panels")
  meta_path <- if (!is.null(args[["meta"]])) args[["meta"]] else file.path("免疫队列", "fig7", "GM", "data", "TCGA", "Survival_SupplementalTable_S1_20171025_xena_sp")
  gene_coords_rda <- if (!is.null(args[["gene-coords-rda"]])) args[["gene-coords-rda"]] else file.path("GBM_CNV!", "genes_GR.rda")

  host <- if (!is.null(args[["xena-host"]])) args[["xena-host"]] else "https://legacy.xenahubs.net"
  expr_dataset <- if (!is.null(args[["xena-expr-dataset"]])) args[["xena-expr-dataset"]] else "TCGA.PANCAN.sampleMap/HiSeqV2"

  genes <- if (!is.null(args[["genes"]])) strsplit(args[["genes"]], ",", fixed = TRUE)[[1]] else c("MYC", "PVT1")
  genes <- unique(trimws(genes))
  method <- if (!is.null(args[["method"]])) tolower(args[["method"]]) else "spearman"
  if (!method %in% c("spearman", "pearson")) stop("`--method` must be spearman or pearson")

  sample_types <- if (!is.null(args[["sample-types"]])) strsplit(args[["sample-types"]], ",", fixed = TRUE)[[1]] else c("01", "03")
  sample_types <- unique(trimws(sample_types))

  panelc_cancers <- if (!is.null(args[["panelc-cancers"]])) strsplit(args[["panelc-cancers"]], ",", fixed = TRUE)[[1]] else c("SKCM", "BLCA")
  panelc_cancers <- unique(trimws(panelc_cancers))
  cytoband_prefix <- if (!is.null(args[["cytoband-prefix"]])) args[["cytoband-prefix"]] else "8q24"
  panelc_top <- if (!is.null(args[["panelc-top"]])) as.integer(args[["panelc-top"]]) else 0L
  if (is.na(panelc_top) || panelc_top < 0) stop("`--panelc-top` must be >= 0 (0 means no limit)")
  panelc_scope <- tolower(args[["panelc-scope"]] %||% "pair") # pair|pancancer
  panelc_summary <- tolower(args[["panelc-summary"]] %||% "median") # median|mean_fisherz
  if (!panelc_scope %in% c("pair", "pancancer")) stop("`--panelc-scope` must be pair or pancancer")
  if (!panelc_summary %in% c("median", "mean_fisherz")) stop("`--panelc-summary` must be median or mean_fisherz")
  only_panelc <- isTRUE(args[["only-panelc"]])
  panelc_resume <- !identical(tolower(args[["panelc-resume"]] %||% "true"), "false")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(outdir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  message("[1/8] Loading meta: ", meta_path)
  meta <- fread(meta_path)
  meta <- meta[, .(sample = canonical_tcga(sample), cancer = as.character(`cancer type abbreviation`))]
  meta <- meta[!is.na(sample) & tcga_sample_type(sample) %in% sample_types]
  meta <- unique(meta)

  message("[2/8] Discovering CNV files under: ", cnv_root)
  gistic <- discover_gistic_files(cnv_root)
  gene_coords <- read_gene_coords(gene_coords_rda)
  gene_coords <- gene_coords[GeneSymbol %in% genes]

  if (only_panelc) {
    message("Running Panel C only (--only-panelc).")
    message("[8/8] Panel C (8q24 ranking)")
  } else {
  message("[3/8] Building continuous CNV table")
  gdc_cache_all <- file.path(cache_dir, "gdc_aliquot_to_sample.tsv")

  cnv_rows <- list()
  for (cancer_code in sort(unique(meta$cancer))) {
    cancer_samples <- unique(meta[cancer == cancer_code, sample])
    if (length(cancer_samples) < 3) next

    if (cancer_code %in% gistic$all$cancer) {
      path <- gistic$all[cancer == cancer_code, path][[1]]
      for (g in genes) {
        v <- read_gistic_gene_vector(path, g, cache_gdc_map = gdc_cache_all)
        if (nrow(v) == 0) next
        dt <- merge(data.table(sample = cancer_samples), v, by = "sample", all = FALSE)
        if (nrow(dt) == 0) next
        cnv_rows[[length(cnv_rows) + 1]] <- dt[, .(sample, cancer = cancer_code, gene = g, cnv = value)]
      }
    } else if (cancer_code == "GBM") {
      gbm_rda <- file.path(cnv_root, "GBM_CNV!", "GBM_CNV_download.rda")
      if (!file.exists(gbm_rda)) next
      e <- new.env(parent = emptyenv())
      load(gbm_rda, envir = e)
      if (!exists("data", envir = e, inherits = FALSE)) next
      seg <- as.data.table(get("data", envir = e))
      for (g in genes) {
        row <- gene_coords[GeneSymbol == g]
        if (nrow(row) != 1) next
        v <- extract_cnv_from_segments(seg, row$Chr, row$Pos, "Sample", "Chromosome", "Start", "End", "Segment_Mean")
        if (nrow(v) == 0) next
        dt <- merge(data.table(sample = cancer_samples), v, by = "sample", all = FALSE)
        cnv_rows[[length(cnv_rows) + 1]] <- dt[, .(sample, cancer = cancer_code, gene = g, cnv = value)]
      }
    } else if (cancer_code == "OV") {
      ov_path <- file.path(cnv_root, "OV_CNV!", "OV_CNV.txt")
      if (!file.exists(ov_path)) next
      ov_seg <- fread(ov_path)
      if (!("GDC_Aliquot" %in% names(ov_seg))) next
      map_dt <- gdc_map_aliquot_to_sample(ov_seg$GDC_Aliquot, cache_path = file.path(cache_dir, "ov_aliquot_to_sample.tsv"))
      setkey(map_dt, aliquot_id)
      ov_seg[, Sample := map_dt[J(GDC_Aliquot), sample_submitter_id]]
      ov_seg <- ov_seg[!is.na(Sample)]
      for (g in genes) {
        row <- gene_coords[GeneSymbol == g]
        if (nrow(row) != 1) next
        v <- extract_cnv_from_segments(ov_seg, row$Chr, row$Pos, "Sample", "Chromosome", "Start", "End", "Segment_Mean")
        if (nrow(v) == 0) next
        dt <- merge(data.table(sample = cancer_samples), v, by = "sample", all = FALSE)
        cnv_rows[[length(cnv_rows) + 1]] <- dt[, .(sample, cancer = cancer_code, gene = g, cnv = value)]
      }
    }
  }
  cnv_cont <- rbindlist(cnv_rows, fill = TRUE)
  fwrite(cnv_cont, file.path(outdir, "cnv_continuous.tsv"), sep = "\t")

  message("[4/8] Fetching expression from Xena (cached): ", expr_dataset)
  samples_all <- sort(unique(cnv_cont$sample))
  expr_wide <- xena_expr_fetch_cached(host, expr_dataset, samples_all, genes, cache_dir = cache_dir)
  expr_long <- melt(expr_wide, id.vars = "sample", variable.name = "gene", value.name = "expr")
  expr_long <- expr_long[!is.na(expr)]

  pan <- merge(cnv_cont, expr_long, by = c("sample", "gene"), all = FALSE)
  pan <- merge(pan, meta, by = c("sample", "cancer"), all = FALSE)
  fwrite(pan, file.path(outdir, "pan_cnv_expr_long.tsv"), sep = "\t")

  message("[5/8] Panel A (hexbin scatter + lm + rho/p)")
  pA <- plot_panel_a_hex(pan, method = method)
  safe_ggsave(file.path(outdir, "PanelA_hex_scatter.pdf"), pA, width = 10.5, height = 4.2, limitsize = FALSE)

  message("[6/8] Panel D (33 cancer × 2 genes bubble)")
  cor_dt <- pan[, {
    ct <- suppressWarnings(cor.test(cnv, expr, method = method, exact = FALSE))
    .(r = unname(ct$estimate[[1]]), p = ct$p.value, n = .N)
  }, by = .(cancer, gene)]
  fwrite(cor_dt, file.path(outdir, "per_cancer_cor.tsv"), sep = "\t")
  pD <- plot_panel_d_bubble(cor_dt)
  safe_ggsave(file.path(outdir, "PanelD_pancancer_bubble.pdf"), pD, width = 6.2, height = 8.8, limitsize = FALSE)

  message("[7/8] Panel B/E (GISTIC categories + ROC)")
  cat_rows <- list()
  for (cancer_code in sort(unique(meta$cancer))) {
    if (!cancer_code %in% gistic$thr$cancer) next
    path <- gistic$thr[cancer == cancer_code, path][[1]]
    cancer_samples <- unique(meta[cancer == cancer_code, sample])
    for (g in genes) {
      v <- read_gistic_gene_vector(path, g, cache_gdc_map = gdc_cache_all)
      if (nrow(v) == 0) next
      dt <- merge(data.table(sample = cancer_samples), v, by = "sample", all = FALSE)
      if (nrow(dt) == 0) next
      dt[, gistic := fifelse(value <= -1, "Loss",
                             fifelse(value < 0.5, "Diploid",
                                     fifelse(value < 1.5, "Gain", "Amp")))]
      cat_rows[[length(cat_rows) + 1]] <- dt[, .(sample, cancer = cancer_code, gene = g, gistic)]
    }
  }
  gistic_cat <- rbindlist(cat_rows, fill = TRUE)
  gistic_cat <- merge(gistic_cat, expr_long, by = c("sample", "gene"), all = FALSE)
  fwrite(gistic_cat, file.path(outdir, "pan_gistic_expr.tsv"), sep = "\t")

  pB <- plot_panel_b_violin(gistic_cat)
  safe_ggsave(file.path(outdir, "PanelB_violin_box.pdf"), pB, width = 9.8, height = 4.2, limitsize = FALSE)

  # Add combined predictor curve for Panel E: MYC_PVT1 = mean(MYC, PVT1)
  if (all(c("MYC", "PVT1") %in% unique(gistic_cat$gene))) {
    wide <- dcast(gistic_cat[gene %in% c("MYC", "PVT1")], sample + cancer + gistic ~ gene, value.var = "expr")
    if (all(c("MYC", "PVT1") %in% names(wide))) {
      wide[, expr := rowMeans(.SD, na.rm = FALSE), .SDcols = c("MYC", "PVT1")]
      comb <- wide[!is.na(expr), .(sample, cancer, gistic, gene = "MYC_PVT1", expr)]
      gistic_cat <- rbindlist(list(gistic_cat, comb), use.names = TRUE, fill = TRUE)
    }
  }

  pE <- plot_panel_e_roc(gistic_cat)
  safe_ggsave(file.path(outdir, "PanelE_ROC.pdf"), pE, width = 5.6, height = 5.0, limitsize = FALSE)

  # Additional ROC: only MP (MYC & PVT1 both Gain/Amp) as positive class
  if (all(c("MYC", "PVT1") %in% unique(gistic_cat$gene))) {
    mp_label <- dcast(
      gistic_cat[gene %in% c("MYC", "PVT1"), .(sample, gene, gistic)],
      sample ~ gene,
      value.var = "gistic"
    )
    if (all(c("MYC", "PVT1") %in% names(mp_label))) {
      mp_label <- mp_label[!is.na(MYC) & !is.na(PVT1)]
      mp_label[, y := as.integer(MYC %in% c("Gain", "Amp") & PVT1 %in% c("Gain", "Amp"))]
      mp_long <- merge(
        gistic_cat[gene %in% c("MYC", "PVT1", "MYC_PVT1"), .(sample, gene, expr)],
        mp_label[, .(sample, y)],
        by = "sample",
        all = FALSE
      )
      pE_mp <- plot_roc_from_y_long(
        dt_long = mp_long,
        title = "Panel E (MP+ only): mRNA -> (MYC & PVT1 both Gain/Amp)",
        direction = "<"
      )
      safe_ggsave(file.path(outdir, "PanelE_ROC_MPonly.pdf"), pE_mp, width = 5.8, height = 5.0, limitsize = FALSE)

      ppv_mp <- mp_long[, {
        m <- ppv_metrics_from_scores(y = y, score = expr, direction = "<")
        if (is.null(m)) return(NULL)
        as.list(m)
      }, by = gene]
      fwrite(ppv_mp, file.path(outdir, "PanelE_PPV_metrics_MPonly.tsv"), sep = "\t")

      # Even stricter positive class: MYC & PVT1 both "Amp"
      mp_label[, y_amp := as.integer(MYC == "Amp" & PVT1 == "Amp")]
      mp_long_amp <- merge(
        gistic_cat[gene %in% c("MYC", "PVT1", "MYC_PVT1"), .(sample, gene, expr)],
        mp_label[, .(sample, y = y_amp)],
        by = "sample",
        all = FALSE
      )
      pE_mp_amp <- plot_roc_from_y_long(
        dt_long = mp_long_amp,
        title = "Panel E (MP Amp only): mRNA -> (MYC & PVT1 both Amp)",
        direction = "<"
      ) +
        scale_color_manual(
          values = c(MYC_PVT1 = "#2E7D32", PVT1 = "#7E57C2", MYC = "#64B5F6"),
          breaks = c("MYC_PVT1", "PVT1", "MYC")
        )
      safe_ggsave(file.path(outdir, "PanelE_ROC_MPAmpOnly.pdf"), pE_mp_amp, width = 5.8, height = 5.0, limitsize = FALSE)

      ppv_mp_amp <- mp_long_amp[, {
        m <- ppv_metrics_from_scores(y = y, score = expr, direction = "<")
        if (is.null(m)) return(NULL)
        as.list(m)
      }, by = gene]
      fwrite(ppv_mp_amp, file.path(outdir, "PanelE_PPV_metrics_MPAmpOnly.tsv"), sep = "\t")
    }
  }

  ppv_dt <- gistic_cat[, {
    y <- fifelse(gistic %in% c("Gain", "Amp"), 1L, 0L)
    m <- ppv_metrics_from_scores(y = y, score = expr, direction = "<")
    if (is.null(m)) return(NULL)
    as.list(m)
  }, by = gene]
  fwrite(ppv_dt, file.path(outdir, "PanelE_PPV_metrics.tsv"), sep = "\t")

  message("[8/8] Panel C (8q24 ranking)")
  }

  gdc_cache_all <- file.path(cache_dir, "gdc_aliquot_to_sample.tsv")
  rank_rows <- list()
  # Choose a stable gene list for 8q24
  gene_list_path <- NULL
  if (any(gistic$all$cancer %in% c("BRCA", "BLCA", "SKCM"))) {
    pref <- c("BRCA", "BLCA", "SKCM")
    gene_list_path <- gistic$all[cancer %in% pref][order(match(cancer, pref))]$path[[1]]
  } else if (nrow(gistic$all) > 0) {
    gene_list_path <- gistic$all$path[[1]]
  }
  if (is.null(gene_list_path)) stop("No GISTIC all_data_by_genes.txt found for Panel C gene list.")

  genes8 <- extract_cytoband_genes(gene_list_path, cytoband_prefix = cytoband_prefix)
  genes8 <- unique(c(genes8, genes))
  genes8 <- genes8[!is.na(genes8) & genes8 != ""]
  if (length(genes8) == 0) stop("No genes found for cytoband prefix: ", cytoband_prefix)

  if (panelc_scope == "pair") {
    for (cancer_code in panelc_cancers) {
      if (!cancer_code %in% gistic$all$cancer) next
      path <- gistic$all[cancer == cancer_code, path][[1]]
      cancer_samples <- unique(meta[cancer == cancer_code, sample])
      if (length(cancer_samples) < 10) next

      cnv_map <- read_gistic_gene_vectors(path, genes8, cache_gdc_map = gdc_cache_all)
      if (length(cnv_map) == 0) next

      expr8_wide <- xena_expr_fetch_cached_batched(host, expr_dataset, cancer_samples, genes8, cache_dir = cache_dir, batch_size = 50)
      expr8_long <- melt(expr8_wide, id.vars = "sample", variable.name = "gene", value.name = "expr")
      expr8_long <- expr8_long[!is.na(expr)]

      for (g in genes8) {
        cnv_dt <- cnv_map[[g]]
        if (is.null(cnv_dt) || nrow(cnv_dt) < 10) next
        merged <- merge(cnv_dt, expr8_long[gene == g, .(sample, expr)], by = "sample", all = FALSE)
        if (nrow(merged) < 10) next
        ct <- suppressWarnings(cor.test(merged$value, merged$expr, method = method, exact = FALSE))
        rank_rows[[length(rank_rows) + 1]] <- data.table(
          cancer = cancer_code,
          gene = g,
          r = as.numeric(unname(ct$estimate[[1]])),
          p = as.numeric(ct$p.value),
          n = as.integer(nrow(merged))
        )
      }
    }
    rank_dt <- rbindlist(rank_rows, fill = TRUE)
    fwrite(rank_dt, file.path(outdir, "panelC_8q24_rank.tsv"), sep = "\t")

    rank_plot <- if (panelc_top == 0L) rank_dt else rank_dt[order(-r), head(.SD, panelc_top), by = cancer]
    pC <- plot_panel_c_combined(rank_plot, gene_coords)
    safe_ggsave(file.path(outdir, "PanelC_8q24_rank.pdf"), pC, width = 10.5, height = 6.6, limitsize = FALSE)
  } else {
    # pan-cancer: compute per-cancer r for all cancers, then summarize across cancers
    per_cancer_path <- file.path(outdir, "panelC_8q24_rank_per_cancer.tsv")
    done_cancers <- character()
    if (panelc_resume && file.exists(per_cancer_path)) {
      old <- fread(per_cancer_path)
      done_cancers <- unique(old$cancer)
      message("Panel C resume: found existing per-cancer results for ", length(done_cancers), " cancers.")
    }

    for (cancer_code in sort(unique(meta$cancer))) {
      if (!cancer_code %in% gistic$all$cancer) next
      if (panelc_resume && cancer_code %in% done_cancers) next
      path <- gistic$all[cancer == cancer_code, path][[1]]
      cancer_samples <- unique(meta[cancer == cancer_code, sample])
      if (length(cancer_samples) < 10) next

      message("  - Panel C (pan-cancer): ", cancer_code, " n=", length(cancer_samples))
      cnv_map <- read_gistic_gene_vectors(path, genes8, cache_gdc_map = gdc_cache_all)
      if (length(cnv_map) == 0) next

      expr8_wide <- xena_expr_fetch_cached_batched(host, expr_dataset, cancer_samples, genes8, cache_dir = cache_dir, batch_size = 50)
      expr8_long <- melt(expr8_wide, id.vars = "sample", variable.name = "gene", value.name = "expr")
      expr8_long <- expr8_long[!is.na(expr)]

      for (g in genes8) {
        cnv_dt <- cnv_map[[g]]
        if (is.null(cnv_dt) || nrow(cnv_dt) < 10) next
        merged <- merge(cnv_dt, expr8_long[gene == g, .(sample, expr)], by = "sample", all = FALSE)
        if (nrow(merged) < 10) next
        ct <- suppressWarnings(cor.test(merged$value, merged$expr, method = method, exact = FALSE))
        rank_rows[[length(rank_rows) + 1]] <- data.table(
          cancer = cancer_code,
          gene = g,
          r = as.numeric(unname(ct$estimate[[1]])),
          p = as.numeric(ct$p.value),
          n = as.integer(nrow(merged))
        )
      }

      if (length(rank_rows) > 0) {
        chunk_dt <- rbindlist(rank_rows, fill = TRUE)
        fwrite(chunk_dt, per_cancer_path, sep = "\t", append = file.exists(per_cancer_path))
        rank_rows <- list()
      }
    }
    per_cancer <- if (file.exists(per_cancer_path)) fread(per_cancer_path) else data.table()

    sum_dt <- summarize_panelc_pancancer(per_cancer, method = panelc_summary)
    setorder(sum_dt, -r)
    fwrite(sum_dt, file.path(outdir, "panelC_8q24_rank_pancancer.tsv"), sep = "\t")

    plot_dt <- if (panelc_top == 0L) sum_dt else head(sum_dt, panelc_top)
    pC <- plot_panel_c_pancancer_combined(plot_dt, gene_coords)
    safe_ggsave(file.path(outdir, "PanelC_8q24_rank_pancancer.pdf"), pC, width = 10.5, height = 6.6, limitsize = FALSE)
  }

  if (!only_panelc) {
    fig <- (pA / pB / pC / pD / pE) + plot_layout(heights = c(1.1, 1.0, 1.8, 2.3, 1.2))
    safe_ggsave(file.path(outdir, "Figure_PanCancer_CNV_mRNA_panels.pdf"), fig, width = 11.0, height = 22.0, limitsize = FALSE)
  }

  message("Done. Outputs in: ", outdir)
}

if (sys.nframe() == 0) {
  main()
}

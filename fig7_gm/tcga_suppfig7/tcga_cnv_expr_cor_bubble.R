suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(httr2)
  library(jsonlite)
  library(digest)
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
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_text_file <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("<0.0001")
  sprintf("%.3g", p)
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

read_gene_coords <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("genes", envir = e, inherits = FALSE)) {
    stop("Gene coords RDA does not contain object named `genes`: ", path)
  }
  genes <- get("genes", envir = e, inherits = FALSE)
  genes <- as.data.table(genes)
  setnames(genes, c("GeneSymbol", "Chr", "Start", "End"))
  genes[, GeneSymbol := as.character(GeneSymbol)]
  genes[, Chr := as.integer(Chr)]
  genes[, Start := as.integer(Start)]
  genes[, End := as.integer(End)]
  genes[, Pos := as.integer(floor((Start + End) / 2))]
  genes
}

## Xena expression fetch (use legacy hub for HiSeqV2)
xena_post <- function(host, query) {
  r <- request(paste0(host, "/data/")) |>
    req_body_raw(query) |>
    req_headers("Content-Type" = "text/plain") |>
    req_perform()
  resp_body_string(r)
}

xena_dataset_probe_values <- function(host, dataset, samples, probes, query_path) {
  qfn <- read_text_file(query_path)
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
    dt[[probes[[i]]]] <- suppressWarnings(as.numeric(values[[i]]))
  }
  dt
}

xena_expr_fetch_cached <- function(host, dataset, samples, probes, cache_dir, query_path) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key_hash <- digest::digest(list(dataset = dataset, samples = samples, probes = probes), algo = "xxhash64")
  key <- paste0("xena_expr_", gsub("[^A-Za-z0-9]+", "_", dataset), "_", key_hash, ".tsv")
  cache_path <- file.path(cache_dir, key)
  if (file.exists(cache_path)) {
    dt <- fread(cache_path)
    if (all(c("sample", probes) %in% names(dt))) return(dt[, c("sample", probes), with = FALSE])
  }
  dt <- xena_dataset_probe_values(host, dataset, samples, probes, query_path = query_path)
  fwrite(dt, cache_path, sep = "\t")
  dt
}

discover_gistic_all_data_files <- function(cnv_root) {
  files <- list.files(
    cnv_root,
    pattern = "all_data_by_genes\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(data.table(cancer = character(), path = character(), ncol = integer()))
  }

  dt <- rbindlist(lapply(files, function(p) {
    parts <- strsplit(normalizePath(p, winslash = "\\", mustWork = FALSE), "\\\\", fixed = FALSE)[[1]]
    cnv_dir <- NA_character_
    for (seg in parts) {
      if (grepl("cnv", seg, ignore.case = TRUE)) {
        cnv_dir <- seg
        break
      }
    }
    cancer <- toupper(sub("^([A-Za-z0-9]+).*$", "\\1", cnv_dir))
    hdr <- tryCatch(readLines(p, n = 1, warn = FALSE), error = function(e) NA_character_)
    ncol <- if (length(hdr) == 1 && !is.na(hdr)) length(strsplit(hdr, "\t", fixed = TRUE)[[1]]) else NA_integer_
    data.table(cancer = cancer, path = p, ncol = as.integer(ncol))
  }), fill = TRUE)

  dt <- dt[!is.na(cancer) & cancer != "" & !is.na(ncol)]
  if (nrow(dt) == 0) {
    return(dt)
  }
  setorder(dt, cancer, -ncol, path)
  dt <- dt[, .SD[1], by = cancer]
  dt[]
}

read_gistic_genes_from_all_data <- function(path, genes, gdc_cache_path = NULL) {
  genes <- unique(as.character(genes))
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) != 1) stop("Failed to read header: ", path)
  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (length(cols) < 5) stop("Unexpected all_data_by_genes format: ", path)
  sample_cols <- cols[4:length(cols)]
  sample <- canonical_tcga(sample_cols)

  uuid_re <- "^[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}$"
  need_map <- is.na(sample) & grepl(uuid_re, sample_cols, perl = TRUE)
  if (any(need_map)) {
    map_dt <- gdc_map_aliquot_to_sample(sample_cols[need_map], cache_path = gdc_cache_path)
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
    chunk <- readLines(con, n = 20000, warn = FALSE)
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

  res <- lapply(genes, function(g) {
    line <- wanted[[g]]
    if (is.na(line)) return(NULL)
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    vals <- suppressWarnings(as.numeric(fields[4:length(fields)]))
    dt <- data.table(sample = sample, cnv = vals)
    dt <- dt[!is.na(sample) & !is.na(cnv)]
    dt <- dt[, .(cnv = mean(cnv)), by = sample]
    dt
  })
  names(res) <- genes
  res
}

combine_expr_mean <- function(expr_by_gene, genes, combined_name = "MYC_PVT1") {
  genes <- unique(as.character(genes))
  genes <- genes[genes %in% names(expr_by_gene)]
  if (length(genes) < 2) return(expr_by_gene)

  dts <- lapply(genes, function(g) {
    dt <- copy(expr_by_gene[[g]])
    setnames(dt, "expr", g)
    dt
  })
  wide <- Reduce(function(a, b) merge(a, b, by = "sample", all = FALSE), dts)
  if (nrow(wide) == 0) return(expr_by_gene)
  wide[, expr := rowMeans(.SD, na.rm = FALSE), .SDcols = genes]
  expr_by_gene[[combined_name]] <- wide[!is.na(expr), .(sample, expr)]
  expr_by_gene
}

combine_cnv_mean <- function(cnv_by_gene, genes, combined_name = "MYC_PVT1") {
  genes <- unique(as.character(genes))
  genes <- genes[genes %in% names(cnv_by_gene)]
  if (length(genes) < 2) return(cnv_by_gene)

  dts <- lapply(genes, function(g) {
    dt <- copy(cnv_by_gene[[g]])
    setnames(dt, "cnv", g)
    dt
  })
  wide <- Reduce(function(a, b) merge(a, b, by = "sample", all = FALSE), dts)
  if (nrow(wide) == 0) return(cnv_by_gene)
  wide[, cnv := rowMeans(.SD, na.rm = FALSE), .SDcols = genes]
  cnv_by_gene[[combined_name]] <- wide[!is.na(cnv), .(sample, cnv)]
  cnv_by_gene
}

extract_cnv_from_segments <- function(seg_dt, chr, pos, sample_col, chr_col, start_col, end_col, value_col) {
  dt <- copy(seg_dt)
  setDT(dt)
  setnames(dt, c(sample_col, chr_col, start_col, end_col, value_col), c("Sample", "Chr", "Start", "End", "Value"), skip_absent = TRUE)
  dt[, Sample := canonical_tcga(Sample)]
  dt <- dt[!is.na(Sample)]

  dt <- dt[as.integer(Chr) == as.integer(chr) & as.integer(Start) <= as.integer(pos) & as.integer(End) >= as.integer(pos)]
  if (nrow(dt) == 0) {
    return(data.table(sample = character(), cnv = numeric()))
  }
  dt[, .(cnv = mean(as.numeric(Value), na.rm = TRUE)), by = Sample][, .(sample = Sample, cnv)]
}

gdc_map_aliquot_to_sample <- function(aliquot_ids, cache_path = NULL, batch_size = 50, sleep_s = 0.2) {
  suppressPackageStartupMessages({
    library(httr2)
    library(jsonlite)
  })

  aliquot_ids <- unique(as.character(aliquot_ids))
  aliquot_ids <- aliquot_ids[!is.na(aliquot_ids) & aliquot_ids != ""]

  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- fread(cache_path)
    if (all(c("aliquot_id", "sample_submitter_id") %in% names(cached))) {
      cached <- unique(cached[, .(aliquot_id = as.character(aliquot_id), sample_submitter_id = as.character(sample_submitter_id))])
      missing <- setdiff(aliquot_ids, cached$aliquot_id)
      if (length(missing) == 0) return(cached)
    } else {
      cached <- data.table(aliquot_id = character(), sample_submitter_id = character())
      missing <- aliquot_ids
    }
  } else {
    cached <- data.table(aliquot_id = character(), sample_submitter_id = character())
    missing <- aliquot_ids
  }

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
    if (length(hits) == 0) {
      Sys.sleep(sleep_s)
      next
    }

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

plot_bubble_rank <- function(df, out_path, title = NULL, color_by = "p", method = "spearman") {
  df <- copy(df)
  df[, neglog10p := -log10(pmax(p, 1e-300))]
  df[, cancer := factor(cancer, levels = df[order(r), cancer])]

  color_by <- tolower(color_by %||% "p")
  if (!color_by %in% c("p", "r")) stop("plot_bubble_rank: `color_by` must be 'p' or 'r'")

  if (color_by == "r") {
    method_label <- if (tolower(method) == "spearman") "Spearman \u03c1" else "Pearson r"
    breaks_r <- c(-1, -0.5, 0, 0.5, 1)

    p <- ggplot(df, aes(x = r, y = cancer, size = n, color = r)) +
      geom_point(alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
      scale_color_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        limits = c(-1, 1),
        breaks = breaks_r,
        labels = sprintf("%.1f", breaks_r),
        oob = scales::squish,
        name = method_label,
        guide = guide_colorbar(
          barheight = grid::unit(55, "mm"),
          barwidth = grid::unit(7, "mm"),
          ticks.colour = "white",
          ticks.linewidth = 0.6,
          frame.colour = NA
        )
      ) +
      scale_size_continuous(name = "n") +
      labs(x = "Correlation (r)", y = "TCGA cancer type", title = title) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      )
  } else {
    p <- ggplot(df, aes(x = r, y = cancer, size = n, color = neglog10p)) +
      geom_point(alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
      scale_color_viridis_c(name = expression(-log[10](p))) +
      scale_size_continuous(name = "n") +
      labs(x = "Correlation (r)", y = "TCGA cancer type", title = title) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, face = "bold")
      )
  }

  safe_ggsave(out_path, p, width = 8.5, height = max(5, 0.23 * nrow(df) + 1.5), limitsize = FALSE)
}

plot_bubble_matrix <- function(df, out_path, title = NULL, method = "spearman") {
  df <- copy(df)
  df[, cancer := factor(cancer, levels = df[order(r), cancer])]
  df[, gene := factor(gene, levels = unique(gene))]

  method_label <- if (tolower(method) == "spearman") "Spearman ρ" else "Pearson r"
  breaks_r <- c(-1, -0.5, 0, 0.5, 1)

  p <- ggplot(df, aes(x = gene, y = cancer)) +
    geom_point(aes(size = n, color = r), alpha = 0.95) +
    scale_color_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = breaks_r,
      labels = sprintf("%.1f", breaks_r),
      oob = scales::squish,
      name = method_label,
      guide = guide_colorbar(
        barheight = grid::unit(60, "mm"),
        barwidth = grid::unit(7, "mm"),
        ticks.colour = "white",
        ticks.linewidth = 0.6,
        frame.colour = NA
      )
    ) +
    scale_size_continuous(name = "n") +
    labs(x = NULL, y = "TCGA cancer type", title = title) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  safe_ggsave(
    out_path,
    p,
    width = max(4.2, 1.8 + 0.6 * length(levels(df$gene))),
    height = max(6, 0.23 * nrow(df) + 1.8),
    limitsize = FALSE
  )
}

plot_bubble <- function(df, out_path, title = NULL, style = "rank", method = "spearman", color_by = "p") {
  style <- tolower(style %||% "rank")
  if (style == "matrix") {
    plot_bubble_matrix(df, out_path, title = title, method = method)
  } else {
    plot_bubble_rank(df, out_path, title = title, color_by = color_by, method = method)
  }
}

plot_scatter <- function(df, out_path, title = NULL, subtitle = NULL) {
  p <- ggplot(df, aes(x = cnv, y = expr)) +
    geom_point(alpha = 0.55, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "#2C7FB8") +
    labs(x = "CNV (segment mean)", y = "mRNA expression", title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  safe_ggsave(out_path, p, width = 6.5, height = 5.0, limitsize = FALSE)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  cnv_root <- if (!is.null(args[["cnv-root"]])) args[["cnv-root"]] else "."
  expr_rdata <- if (!is.null(args[["expr-rdata"]])) args[["expr-rdata"]] else file.path("免疫队列", "fig7", "GM", "data", "TCGA", "bulkExpMatrix.Rdata")
  meta_path <- if (!is.null(args[["meta"]])) args[["meta"]] else file.path("免疫队列", "fig7", "GM", "data", "TCGA", "Survival_SupplementalTable_S1_20171025_xena_sp")
  gene_coords_rda <- if (!is.null(args[["gene-coords-rda"]])) args[["gene-coords-rda"]] else file.path("GBM_CNV!", "genes_GR.rda")
  outdir <- if (!is.null(args[["outdir"]])) args[["outdir"]] else file.path("results", "cnv_expr_cor")
  genes <- if (!is.null(args[["genes"]])) strsplit(args[["genes"]], ",", fixed = TRUE)[[1]] else c("MYC", "PVT1")
  genes <- unique(trimws(genes))
  method <- if (!is.null(args[["method"]])) tolower(args[["method"]]) else "spearman"
  expr_transform <- if (!is.null(args[["expr-transform"]])) tolower(args[["expr-transform"]]) else "none"
  do_scatter <- isTRUE(args[["scatter"]])
  sample_types <- if (!is.null(args[["sample-types"]])) strsplit(args[["sample-types"]], ",", fixed = TRUE)[[1]] else c("01", "03")
  sample_types <- unique(trimws(sample_types))
  if (any(!grepl("^[0-9]{2}$", sample_types))) {
    stop("`--sample-types` must be a comma-separated list of 2-digit codes, e.g. 01,03")
  }

  expr_source <- tolower(args[["expr-source"]] %||% "auto") # auto|local|xena
  xena_host <- args[["xena-host"]] %||% "https://legacy.xenahubs.net"
  xena_dataset <- args[["xena-dataset"]] %||% "TCGA.PANCAN.sampleMap/HiSeqV2"
  xena_query <- args[["xena-query"]] %||% file.path("xena_queries", "datasetProbeValues.xq")
  combine <- tolower(args[["combine"]] %||% "none") # none|mean
  combine_name <- args[["combine-name"]] %||% "MYC_PVT1"
  combined_only <- isTRUE(args[["combined-only"]])
  bubble_out <- args[["bubble-out"]] %||% NULL
  bubble_style <- tolower(args[["bubble-style"]] %||% "rank") # rank|matrix
  bubble_color <- tolower(args[["bubble-color"]] %||% if (bubble_style == "rank" && method == "pearson") "r" else "p") # p|r

  if (!method %in% c("spearman", "pearson")) {
    stop("`--method` must be spearman or pearson")
  }
  if (!expr_transform %in% c("none", "log2p1")) {
    stop("`--expr-transform` must be none or log2p1")
  }
  if (!expr_source %in% c("auto", "local", "xena")) {
    stop("`--expr-source` must be auto, local, or xena")
  }
  if (!combine %in% c("none", "mean")) {
    stop("`--combine` must be none or mean")
  }
  if (!bubble_style %in% c("rank", "matrix")) {
    stop("`--bubble-style` must be rank or matrix")
  }
  if (!bubble_color %in% c("p", "r")) {
    stop("`--bubble-color` must be p or r")
  }
  if (expr_source %in% c("xena", "auto") && !file.exists(xena_query)) {
    stop("Missing Xena query file: ", xena_query)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("[1/6] Loading gene coordinates: ", gene_coords_rda)
  gene_coords <- read_gene_coords(gene_coords_rda)
  gene_coords <- gene_coords[GeneSymbol %in% genes]
  if (nrow(gene_coords) == 0) stop("No requested genes found in gene coords: ", paste(genes, collapse = ", "))
  missing_coords <- setdiff(genes, gene_coords$GeneSymbol)
  if (length(missing_coords) > 0) {
    warning("Missing gene coords for: ", paste(missing_coords, collapse = ", "))
  }

  expr_by_gene <- list()
  if (expr_source %in% c("local", "auto")) {
    message("[2/6] Loading expression matrix: ", expr_rdata)
    e <- new.env(parent = emptyenv())
    load(expr_rdata, envir = e)
    expr_name <- if (exists("bulkExpMatrix", envir = e, inherits = FALSE)) {
      "bulkExpMatrix"
    } else {
      objs <- ls(envir = e)
      mats <- objs[vapply(objs, function(nm) is.matrix(get(nm, envir = e)), logical(1))]
      if (length(mats) == 0) stop("No matrix object found in: ", expr_rdata)
      mats[[1]]
    }
    expr_mat <- get(expr_name, envir = e)

    expr_samples <- canonical_tcga(colnames(expr_mat))
    expr_keep <- !is.na(expr_samples) & tcga_sample_type(expr_samples) %in% sample_types

    for (g in genes) {
      if (!g %in% rownames(expr_mat)) {
        warning("Expression matrix missing gene: ", g)
        next
      }
      v <- as.numeric(expr_mat[g, ])
      if (expr_transform == "log2p1") v <- log2(v + 1)
      dt <- data.table(sample = expr_samples, expr = v)
      dt <- dt[expr_keep & !is.na(expr)]
      dt <- dt[, .(expr = mean(expr)), by = sample]
      expr_by_gene[[g]] <- dt
    }
  } else {
    message("[2/6] Expression source: Xena (skip local Rdata)")
  }

  message("[3/6] Loading sample metadata: ", meta_path)
  meta <- fread(meta_path)
  if (!all(c("sample", "cancer type abbreviation") %in% names(meta))) {
    stop("Meta file missing required columns: sample, cancer type abbreviation")
  }
  meta <- meta[, .(sample = canonical_tcga(sample), cancer = as.character(`cancer type abbreviation`))]
  meta <- meta[!is.na(sample) & tcga_sample_type(sample) %in% sample_types]
  meta <- unique(meta)

  if (combined_only && combine != "mean") {
    stop("`--combined-only` requires `--combine mean`.")
  }

  if (expr_source %in% c("xena", "auto")) {
    missing_genes <- setdiff(genes, names(expr_by_gene))
    if (expr_source == "xena" || length(missing_genes) > 0) {
      cache_dir <- file.path(outdir, "cache")
      samples_all <- sort(unique(meta$sample))
      xdt <- xena_expr_fetch_cached(
        host = xena_host,
        dataset = xena_dataset,
        samples = samples_all,
        probes = genes,
        cache_dir = cache_dir,
        query_path = xena_query
      )
      for (g in genes) {
        if (!g %in% names(xdt)) next
        dt <- xdt[, .(sample = canonical_tcga(sample), expr = as.numeric(get(g)))]
        dt <- dt[!is.na(sample) & !is.na(expr)]
        dt <- dt[, .(expr = mean(expr)), by = sample]
        expr_by_gene[[g]] <- dt
      }
    }
  }

  if (length(expr_by_gene) == 0) {
    stop("No expression data available for requested genes. Try `--expr-source xena`.")
  }
  if (combine == "mean") {
    expr_by_gene <- combine_expr_mean(expr_by_gene, genes, combined_name = combine_name)
  }

  analysis_genes <- if (combined_only) combine_name else unique(c(genes, if (combine == "mean") combine_name else character()))

  message("[4/6] Discovering CNV inputs under: ", cnv_root)
  gistic_files <- discover_gistic_all_data_files(cnv_root)
  cnv_sources <- unique(c(gistic_files$cancer, "GBM", "OV"))

  results <- list()
  skipped <- list()

  for (cancer_code in sort(unique(meta$cancer))) {
    if (!cancer_code %in% cnv_sources) next

    cancer_samples <- meta[cancer == cancer_code, sample]
    if (length(cancer_samples) < 3) next

    for (g in analysis_genes) {
      if (is.null(expr_by_gene[[g]])) next
      overlap_n <- sum(cancer_samples %in% expr_by_gene[[g]]$sample)
      if (overlap_n < 3) {
        skipped[[length(skipped) + 1]] <- data.table(
          cancer = cancer_code,
          gene = g,
          reason = "insufficient_expression_overlap",
          n_expr_overlap = as.integer(overlap_n)
        )
      }
    }

    cnv_by_gene <- list()

    if (cancer_code %in% gistic_files$cancer) {
      path <- gistic_files[cancer == cancer_code, path][[1]]
      message("  - CNV: ", cancer_code, " (GISTIC all_data_by_genes): ", path)
      cnv_by_gene <- read_gistic_genes_from_all_data(
        path,
        genes,
        gdc_cache_path = file.path(outdir, "cache", "gdc_aliquot_to_sample.tsv")
      )
    } else if (cancer_code == "GBM") {
      gbm_rda <- file.path(cnv_root, "GBM_CNV!", "GBM_CNV_download.rda")
      if (!file.exists(gbm_rda)) {
        warning("Missing GBM CNV RDA: ", gbm_rda, " (skipping GBM)")
        next
      }
      message("  - CNV: GBM (segments): ", gbm_rda)
      gbm_env <- new.env(parent = emptyenv())
      load(gbm_rda, envir = gbm_env)
      if (!exists("data", envir = gbm_env, inherits = FALSE)) {
        warning("GBM RDA missing object `data`: ", gbm_rda, " (skipping GBM)")
        next
      }
      seg <- as.data.table(get("data", envir = gbm_env))
      for (g in genes) {
        row <- gene_coords[GeneSymbol == g]
        if (nrow(row) != 1) next
        cnv_by_gene[[g]] <- extract_cnv_from_segments(
          seg_dt = seg,
          chr = row$Chr,
          pos = row$Pos,
          sample_col = "Sample",
          chr_col = "Chromosome",
          start_col = "Start",
          end_col = "End",
          value_col = "Segment_Mean"
        )
      }
    } else if (cancer_code == "OV") {
      ov_path <- file.path(cnv_root, "OV_CNV!", "OV_CNV.txt")
      if (!file.exists(ov_path)) {
        warning("Missing OV CNV file: ", ov_path, " (skipping OV)")
        next
      }
      message("  - CNV: OV (segments + GDC mapping): ", ov_path)
      ov_seg <- fread(ov_path)
      if (!all(c("GDC_Aliquot", "Chromosome", "Start", "End", "Segment_Mean") %in% names(ov_seg))) {
        warning("OV CNV file format unexpected: ", ov_path, " (skipping OV)")
        next
      }
      cache_path <- file.path(outdir, "cache", "ov_aliquot_to_sample.tsv")
      map_dt <- gdc_map_aliquot_to_sample(ov_seg$GDC_Aliquot, cache_path = cache_path)
      setkey(map_dt, aliquot_id)
      ov_seg[, Sample := map_dt[J(GDC_Aliquot), sample_submitter_id]]
      ov_seg <- ov_seg[!is.na(Sample)]
      for (g in genes) {
        row <- gene_coords[GeneSymbol == g]
        if (nrow(row) != 1) next
        cnv_by_gene[[g]] <- extract_cnv_from_segments(
          seg_dt = ov_seg,
          chr = row$Chr,
          pos = row$Pos,
          sample_col = "Sample",
          chr_col = "Chromosome",
          start_col = "Start",
          end_col = "End",
          value_col = "Segment_Mean"
        )
      }
    }

    if (combine == "mean" && combine_name %in% analysis_genes) {
      cnv_by_gene <- combine_cnv_mean(cnv_by_gene, genes, combined_name = combine_name)
    }

    for (g in analysis_genes) {
      if (is.null(expr_by_gene[[g]])) next
      cnv_dt <- cnv_by_gene[[g]]
      if (is.null(cnv_dt) || nrow(cnv_dt) == 0) next

      merged <- merge(
        data.table(sample = cancer_samples),
        merge(cnv_dt, expr_by_gene[[g]], by = "sample", all = FALSE),
        by = "sample",
        all = FALSE
      )
      merged <- merged[!is.na(cnv) & !is.na(expr)]
      n <- nrow(merged)
      if (n < 3) next

      ct <- suppressWarnings(cor.test(merged$cnv, merged$expr, method = method, exact = FALSE))
      r <- unname(ct$estimate[[1]])
      p <- ct$p.value

      results[[length(results) + 1]] <- data.table(
        cancer = cancer_code,
        gene = g,
        method = method,
        r = as.numeric(r),
        p = as.numeric(p),
        n = as.integer(n)
      )

      if (do_scatter) {
        scatter_dir <- file.path(outdir, "scatter")
        dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
        out_pdf <- file.path(scatter_dir, paste0("scatter_", cancer_code, "_", g, "_", method, ".pdf"))
        subtitle <- sprintf("r=%.3f, p=%.3g, n=%d", r, p, n)
        plot_scatter(merged, out_pdf, title = paste0(cancer_code, " / ", g), subtitle = subtitle)
      }
    }
  }

  res <- if (length(results) == 0) {
    data.table(cancer = character(), gene = character(), method = character(), r = numeric(), p = numeric(), n = integer())
  } else {
    rbindlist(results)
  }

  out_tsv <- file.path(outdir, paste0("tcga_cnv_expr_cor_", method, ".tsv"))
  fwrite(res, out_tsv, sep = "\t")

  message("[5/6] Writing results: ", out_tsv)

  if (length(skipped) > 0) {
    skip_dt <- unique(rbindlist(skipped, fill = TRUE))
    skip_path <- file.path(outdir, paste0("tcga_cnv_expr_cor_skipped_", method, ".tsv"))
    fwrite(skip_dt, skip_path, sep = "\t")
    message("      Skipped summary: ", skip_path)
  }

  message("[6/6] Plotting bubble plots")
  for (g in unique(res$gene)) {
    df <- res[gene == g]
    if (nrow(df) == 0) next
    df <- df[order(r)]
    out_pdf <- file.path(outdir, paste0("bubble_", g, "_", method, ".pdf"))
    if (!is.null(bubble_out) && length(unique(res$gene)) == 1) {
      out_pdf <- bubble_out
    }
    plot_bubble(
      df,
      out_pdf,
      title = paste0(g, ": CNV vs mRNA (", method, ")"),
      style = bubble_style,
      method = method,
      color_by = bubble_color
    )
  }

  message("Done.")
}

main()

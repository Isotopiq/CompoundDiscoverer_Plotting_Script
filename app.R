
suppressPackageStartupMessages({
  library(shiny)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(DT)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

sanitize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("\\.raw$", "", x)
  x <- gsub("^area:\\s*", "", x)
  x <- gsub("\\s*\\([^\\)]*\\)", "", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_feature_label <- function(df) {
  if ("Name" %in% names(df)) out <- as.character(df$Name)
  else if ("Compound Name" %in% names(df)) out <- as.character(df[["Compound Name"]])
  else if ("Compound" %in% names(df)) out <- as.character(df$Compound)
  else {
    mz_col <- names(df)[tolower(names(df)) %in% c("m/z", "mz")]
    rt_col <- names(df)[tolower(names(df)) %in% c("rt [min]", "rt")]
    if (length(mz_col) > 0 && length(rt_col) > 0) {
      out <- paste0("mz_", round(as.numeric(df[[mz_col[1]]]), 4), "_rt_", round(as.numeric(df[[rt_col[1]]]), 2))
    } else out <- paste0("feature_", seq_len(nrow(df)))
  }
  out[is.na(out) | out == ""] <- paste0("feature_", which(is.na(out) | out == ""))
  out
}

read_cd_file <- function(path, sheet = NULL) {
  sheets <- readxl::excel_sheets(path)
  target_sheet <- sheet %||% if ("Compounds" %in% sheets) "Compounds" else sheets[[1]]
  df <- readxl::read_excel(path, sheet = target_sheet)
  df <- as.data.frame(df, check.names = FALSE)

  area_cols <- grep("^Area:\\s*", names(df), value = TRUE)
  if (length(area_cols) < 2) stop("Could not find at least two sample columns beginning with 'Area: '.")
  mat <- df[, area_cols, drop = FALSE]
  for (nm in area_cols) mat[[nm]] <- suppressWarnings(as.numeric(mat[[nm]]))

  sample_names <- sanitize_name(area_cols)
  names(mat) <- sample_names
  feature_label <- safe_feature_label(df)
  feature_ids <- make.unique(feature_label)
  rownames(mat) <- feature_ids

  annotation_df <- df[, setdiff(names(df), area_cols), drop = FALSE]
  annotation_df$Feature_Label <- feature_label
  annotation_df$Feature_ID <- feature_ids

  list(
    raw = df,
    annotation = annotation_df,
    matrix = as.matrix(mat),
    sample_map = data.frame(original = area_cols, sample = sample_names, stringsAsFactors = FALSE),
    sheet = target_sheet,
    sheets = sheets
  )
}

infer_groups <- function(samples, group_a_name, group_b_name, group_a_pattern = "", group_b_pattern = "") {
  grp <- rep(NA_character_, length(samples))
  if (nzchar(group_a_pattern)) grp[stringr::str_detect(samples, regex(group_a_pattern, ignore_case = TRUE))] <- group_a_name
  if (nzchar(group_b_pattern)) grp[stringr::str_detect(samples, regex(group_b_pattern, ignore_case = TRUE))] <- group_b_name
  data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
}

normalize_for_heatmap <- function(mat, log2_transform = TRUE, pseudocount = 1) {
  x <- mat
  if (log2_transform) x <- log2(x + pseudocount)
  x <- t(scale(t(x)))
  x[is.na(x)] <- 0
  x[x > 2] <- 2
  x[x < -2] <- -2
  x
}

welch_p <- function(x, y) tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)

sig_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 1e-4) return("****")
  if (p <= 1e-3) return("***")
  if (p <= 1e-2) return("**")
  if (p <= 5e-2) return("*")
  "ns"
}

compute_stats <- function(mat, meta, group_a, group_b, pseudocount = 1) {
  meta2 <- meta %>% dplyr::filter(group %in% c(group_a, group_b))
  keep_samples <- intersect(colnames(mat), meta2$sample)
  meta2 <- meta2 %>% dplyr::filter(sample %in% keep_samples)

  a_samples <- meta2 %>% dplyr::filter(group == group_a) %>% dplyr::pull(sample)
  b_samples <- meta2 %>% dplyr::filter(group == group_b) %>% dplyr::pull(sample)
  if (length(a_samples) < 2 || length(b_samples) < 2) stop("Each comparison group needs at least 2 matched samples.")

  submat <- mat[, c(a_samples, b_samples), drop = FALSE]
  rows <- vector("list", nrow(submat))
  for (i in seq_len(nrow(submat))) {
    x <- as.numeric(submat[i, a_samples]); y <- as.numeric(submat[i, b_samples])
    mean_a <- mean(x, na.rm = TRUE); mean_b <- mean(y, na.rm = TRUE)
    rows[[i]] <- data.frame(
      mean_a = mean_a,
      mean_b = mean_b,
      log2fc = log2((mean_b + pseudocount) / (mean_a + pseudocount)),
      p_value = welch_p(x, y),
      stringsAsFactors = FALSE
    )
  }
  res <- dplyr::bind_rows(rows)
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res$Feature_ID <- rownames(submat)
  res
}

build_plot_table <- function(stats_df, annotation_df) {
  out <- dplyr::left_join(stats_df, annotation_df, by = "Feature_ID")
  label_col <- NULL
  for (nm in c("Name", "Compound Name", "Compound", "Feature_Label")) {
    if (nm %in% names(out)) { label_col <- nm; break }
  }
  if (is.null(label_col)) label_col <- "Feature_ID"
  out$Label <- as.character(out[[label_col]])
  out$Label[is.na(out$Label) | out$Label == ""] <- out$Feature_ID[is.na(out$Label) | out$Label == ""]
  out
}

plot_pca <- function(mat, meta, group_colors, point_size = 4) {
  shared <- intersect(colnames(mat), meta$sample)
  m <- t(log2(mat[, shared, drop = FALSE] + 1))
  pca <- prcomp(m, center = TRUE, scale. = TRUE)
  df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = meta$group[match(rownames(pca$x), meta$sample)], stringsAsFactors = FALSE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  ggplot(df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = point_size) +
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.08, color = NA, show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 12) +
    labs(title = "PCA", x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2]), color = "Group") +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
}

plot_volcano <- function(plot_df, group_a, group_b, group_a_color, group_b_color,
                         fc_threshold = 1, q_threshold = 0.05, label_top_n = 15,
                         identified_only_labels = TRUE) {
  df <- plot_df
  df$neglog10q <- -log10(pmax(df$q_value, 1e-300))
  df$status <- "Not sig"
  df$status[df$log2fc >= fc_threshold & df$q_value <= q_threshold] <- paste0(group_b, " higher")
  df$status[df$log2fc <= -fc_threshold & df$q_value <= q_threshold] <- paste0(group_a, " higher")
  label_candidates <- df %>% dplyr::filter(q_value <= q_threshold, abs(log2fc) >= fc_threshold)
  if (identified_only_labels) label_candidates <- label_candidates %>% dplyr::filter(!is.na(Label), Label != "", !grepl("^feature_", Label, ignore.case = TRUE))
  label_df <- label_candidates %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc))) %>% head(label_top_n)

  ggplot(df, aes(x = log2fc, y = neglog10q)) +
    geom_point(aes(color = status), alpha = 0.8, size = 2.2) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(q_threshold), linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = label_df, aes(label = Label), size = 3.2, min.segment.length = 0,
      max.overlaps = Inf, box.padding = 0.35, point.padding = 0.2, segment.alpha = 0.6, seed = 123) +
    scale_color_manual(values = c("Not sig" = "grey70", setNames(group_a_color, paste0(group_a, " higher")), setNames(group_b_color, paste0(group_b, " higher")))) +
    theme_bw(base_size = 12) +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = expression(-log[10]("q-value")), color = NULL) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
}

plot_sample_summary <- function(meta, values, title, ylab, group_colors) {
  df <- data.frame(sample = names(values), value = as.numeric(values), stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, meta, by = "sample")
  ggplot(df, aes(x = sample, y = value, fill = group)) +
    geom_col() +
    scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 12) +
    labs(title = title, x = NULL, y = ylab) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
}

save_heatmap <- function(heatmap_mat, meta, out_file, group_colors) {
  sample_order_df <- meta %>% dplyr::filter(sample %in% colnames(heatmap_mat))
  heatmap_mat <- heatmap_mat[, sample_order_df$sample, drop = FALSE]
  ann_col <- data.frame(class = sample_order_df$group, stringsAsFactors = FALSE)
  rownames(ann_col) <- sample_order_df$sample
  gp <- group_colors[names(group_colors) %in% unique(sample_order_df$group)]
  ann_colors <- list(class = gp)
  n_rows <- max(1, nrow(heatmap_mat)); n_cols <- max(1, ncol(heatmap_mat)); max_label_chars <- max(nchar(rownames(heatmap_mat)), na.rm = TRUE)
  width_in <- min(20, max(10, 4 + n_cols * 0.65 + max_label_chars * 0.09 + 1.5))
  height_in <- min(18, max(8, 2.8 + n_rows * 0.24 + 1.2))
  pheatmap::pheatmap(
    heatmap_mat,
    color = colorRampPalette(c("#6BAED6", "#F7F7F7", "#CB181D"))(100),
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = TRUE, show_colnames = TRUE,
    annotation_col = ann_col, annotation_colors = ann_colors,
    annotation_legend = TRUE, legend = TRUE,
    border_color = "grey65",
    fontsize = 10, fontsize_row = 9, fontsize_col = 10,
    angle_col = 270, treeheight_row = 55, treeheight_col = 55,
    cellwidth = 28, cellheight = 15,
    filename = out_file, width = width_in, height = height_in
  )
}

plot_individual_bar <- function(feature_name, label_name, mat, meta, group_a, group_b, group_a_color, group_b_color, out_file) {
  samples_a <- meta %>% dplyr::filter(group == group_a, sample %in% colnames(mat)) %>% dplyr::pull(sample)
  samples_b <- meta %>% dplyr::filter(group == group_b, sample %in% colnames(mat)) %>% dplyr::pull(sample)
  vals_a <- as.numeric(mat[feature_name, samples_a]); vals_b <- as.numeric(mat[feature_name, samples_b])
  df <- data.frame(value = c(vals_a, vals_b), group = c(rep(group_a, length(vals_a)), rep(group_b, length(vals_b))),
                   sample = c(samples_a, samples_b), stringsAsFactors = FALSE)
  pval <- welch_p(vals_a, vals_b); l2fc <- log2((mean(vals_b) + 1) / (mean(vals_a) + 1))
  ymax <- max(df$value, na.rm = TRUE); yline <- ymax * 1.10; ystar <- ymax * 1.16
  fill_vals <- c(group_a_color, group_b_color); names(fill_vals) <- c(group_a, group_b)

  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.7) +
    geom_jitter(aes(color = group), width = 0.08, size = 3, alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = fill_vals) +
    scale_color_manual(values = fill_vals) +
    theme_bw(base_size = 13) +
    labs(title = label_name, subtitle = sprintf("log2FC = %.3f | p = %.3g", l2fc, pval), x = NULL, y = "Peak Area") +
    theme(plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)),
          plot.subtitle = element_text(size = 11, margin = margin(b = 12)),
          panel.grid.minor = element_blank(), legend.position = "none") +
    annotate("segment", x = 1, xend = 2, y = yline, yend = yline, linewidth = 0.8) +
    annotate("segment", x = 1, xend = 1, y = yline * 0.98, yend = yline, linewidth = 0.8) +
    annotate("segment", x = 2, xend = 2, y = yline * 0.98, yend = yline, linewidth = 0.8) +
    annotate("text", x = 1.5, y = ystar, label = sig_stars(pval), size = 7, fontface = "bold") +
    expand_limits(y = ystar * 1.05)
  ggsave(out_file, p, width = 6.2, height = 6.0, dpi = 220)
}

img_panel <- function(title, src, note = NULL) {
  tags$div(class = "card-shell",
    tags$div(class = "card-head", title),
    if (!is.null(note)) tags$p(class = "card-note", note),
    tags$img(src = src, class = "plot-image")
  )
}

ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$script(src = "https://cdn.tailwindcss.com"),
    tags$style(HTML("
      body { background:#f8fafc; color:#0f172a; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }
      .app-shell { max-width: 1700px; margin: 0 auto; padding: 24px; }
      .hero { background: linear-gradient(180deg,#ffffff 0%,#f8fafc 100%); border:1px solid #e2e8f0; border-radius: 22px; padding: 24px 28px; box-shadow: 0 10px 30px rgba(15,23,42,.04); margin-bottom: 20px; }
      .hero-title { font-size: 28px; font-weight: 800; letter-spacing: -.02em; margin:0; }
      .hero-sub { margin-top: 8px; color:#475569; }
      .layout-grid { display:grid; grid-template-columns: 360px minmax(0,1fr); gap:20px; align-items:start; }
      .side-panel { position: sticky; top: 20px; background:#fff; border:1px solid #e2e8f0; border-radius: 20px; padding: 20px; box-shadow: 0 10px 30px rgba(15,23,42,.04); }
      .content-panel { min-width:0; }
      .section-title { font-size: 13px; font-weight: 700; text-transform: uppercase; letter-spacing: .08em; color:#64748b; margin-bottom: 12px; }
      .form-control, .selectize-input { border-radius: 12px !important; border:1px solid #cbd5e1 !important; box-shadow:none !important; min-height: 44px !important; padding-top: 10px !important; padding-bottom: 10px !important; }
      .btn-primary { width:100%; border:none; border-radius:14px; background:#2563eb; padding:12px 16px; font-weight:700; color:#fff; }
      .btn-default { width:100%; border-radius:14px; }
      .shiny-input-container { width:100% !important; margin-bottom: 14px !important; }
      .nav-pills { display:flex; flex-wrap:wrap; gap:10px; margin-bottom:16px; }
      .nav-pills > li { float:none !important; }
      .nav-pills > li > a { border-radius: 999px; padding: 10px 16px; border:1px solid #dbeafe; background:#eff6ff; color:#1d4ed8; font-weight:600; }
      .nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus { background:#2563eb; border-color:#2563eb; color:#fff; }
      .card-shell { background:#fff; border:1px solid #e2e8f0; border-radius:20px; padding:18px; box-shadow:0 10px 30px rgba(15,23,42,.04); margin-bottom:18px; overflow:hidden; }
      .card-head { font-size:18px; font-weight:700; margin-bottom:8px; }
      .card-note { color:#64748b; margin-bottom:12px; }
      .plot-image { width:100%; height:auto; display:block; border-radius:14px; border:1px solid #e5e7eb; background:#fff; }
      .barplot-grid { display:grid; grid-template-columns: repeat(auto-fit, minmax(420px,1fr)); gap:16px; }
      .barplot-item { background:#fff; border:1px solid #e2e8f0; border-radius:18px; padding:12px; box-shadow:0 10px 30px rgba(15,23,42,.04); }
      .barplot-item img { width:100%; height:auto; display:block; border-radius:12px; }
      .dataTables_wrapper { background:#fff; border:1px solid #e2e8f0; border-radius:20px; padding:14px; box-shadow:0 10px 30px rgba(15,23,42,.04); }
      .progress { height: 12px; border-radius: 999px; overflow:hidden; background:#e2e8f0; }
      .progress-bar { background:#2563eb; }
      @media (max-width: 1100px) { .layout-grid { grid-template-columns: 1fr; } .side-panel { position: static; } }
    "))
  ),
  div(class = "app-shell",
    div(class = "hero",
      h1(class = "hero-title", "Compound Discoverer Metabolomics Explorer"),
      p(class = "hero-sub", "Modern light dashboard wrapper with PCA, volcano plot, clustered heatmap, sample summaries, differential table, ZIP export, and per-metabolite statistical bar plots.")
    ),
    div(class = "layout-grid",
      div(class = "side-panel",
        div(class = "section-title", "Input & analysis"),
        fileInput("xlsx", "Upload Compound Discoverer Excel file", accept = c(".xlsx", ".xls")),
        uiOutput("sheet_ui"),
        fileInput("metadata", "Optional metadata CSV", accept = ".csv"),
        textInput("group_a", "Group A name", value = "Control"),
        textInput("group_b", "Group B name", value = "KO"),
        textInput("pattern_a", "Auto-assign regex for Group A", value = "yl_001|yl_002|yl_003"),
        textInput("pattern_b", "Auto-assign regex for Group B", value = "yl_004|yl_005|yl_006|yl_007|yl_008|yl_009|yl_010"),
        checkboxInput("keep_other_groups", "Keep other groups in PCA/heatmap if present", FALSE),
        checkboxInput("identified_only_heatmap", "Heatmap: identified compounds only", TRUE),
        checkboxInput("identified_only_labels", "Volcano labels: identified compounds only", TRUE),
        tags$hr(),
        div(class = "section-title", "Thresholds & display"),
        numericInput("fc_threshold", "Volcano |log2FC| threshold", value = 1, min = 0, step = 0.1),
        numericInput("q_threshold", "Volcano q-value threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("volcano_top_n", "Volcano labels top N", value = 15, min = 1, step = 1),
        numericInput("heatmap_top_n", "Heatmap top N features", value = 40, min = 5, step = 1),
        numericInput("bar_top_n", "Individual bar plots top N", value = 20, min = 1, step = 1),
        textInput("group_a_color", "Group A color", value = "#E91E63"),
        textInput("group_b_color", "Group B color", value = "#43B649"),
        actionButton("run", "Run Analysis", class = "btn btn-primary"),
        br(), br(),
        downloadButton("download_zip", "Download Results ZIP", class = "btn btn-default")
      ),
      div(class = "content-panel",
        tabsetPanel(
          type = "pills",
          tabPanel("Overview",
            div(class = "card-shell", div(class = "card-head", "Detected samples"), DTOutput("sample_table")),
            div(class = "card-shell", div(class = "card-head", "Metadata used"), DTOutput("meta_table"))
          ),
          tabPanel("PCA", uiOutput("pca_ui")),
          tabPanel("Volcano", uiOutput("volcano_ui")),
          tabPanel("Heatmap", uiOutput("heatmap_ui")),
          tabPanel("Summaries", uiOutput("tic_ui"), uiOutput("detect_ui"), uiOutput("median_ui")),
          tabPanel("Top table", DTOutput("stats_table")),
          tabPanel("Bar plots", uiOutput("barplot_gallery"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(cd=NULL, meta=NULL, stats=NULL, plot_df=NULL, outdir=NULL, resource_prefix=NULL)

  output$sheet_ui <- renderUI({
    req(input$xlsx)
    sheets <- readxl::excel_sheets(input$xlsx$datapath)
    selected_sheet <- if ("Compounds" %in% sheets) "Compounds" else sheets[[1]]
    selectInput("sheet", "Worksheet", choices = sheets, selected = selected_sheet)
  })

  parsed_data <- eventReactive(list(input$xlsx, input$sheet), {
    req(input$xlsx, input$sheet)
    read_cd_file(input$xlsx$datapath, sheet = input$sheet)
  }, ignoreNULL = FALSE)

  metadata_used <- reactive({
    cd <- parsed_data(); req(cd)
    if (!is.null(input$metadata)) {
      md <- read.csv(input$metadata$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      names(md) <- tolower(names(md))
      if (!all(c("sample", "group") %in% names(md))) stop("Metadata CSV must contain columns named 'sample' and 'group'.")
      md$sample <- sanitize_name(md$sample)
      md$group <- as.character(md$group)
      md
    } else {
      infer_groups(cd$sample_map$sample, input$group_a, input$group_b, input$pattern_a, input$pattern_b)
    }
  })

  observe({
    cd <- try(parsed_data(), silent = TRUE)
    if (!inherits(cd, "try-error")) output$sample_table <- renderDT(datatable(cd$sample_map, options = list(scrollX = TRUE, pageLength = 12)))
  })

  observe({
    md <- try(metadata_used(), silent = TRUE)
    if (!inherits(md, "try-error")) output$meta_table <- renderDT(datatable(md, options = list(scrollX = TRUE, pageLength = 12)))
  })

  analysis <- eventReactive(input$run, {
    withProgress(message = "Running analysis", value = 0, {
      cd <- parsed_data(); meta <- metadata_used()

      incProgress(0.08, detail = "Validating samples")
      meta <- meta %>% dplyr::filter(sample %in% colnames(cd$matrix))
      if (!isTRUE(input$keep_other_groups)) meta <- meta %>% dplyr::filter(group %in% c(input$group_a, input$group_b))

      incProgress(0.14, detail = "Computing statistics")
      stats_df <- compute_stats(cd$matrix, meta, input$group_a, input$group_b)
      plot_df <- build_plot_table(stats_df, cd$annotation)

      incProgress(0.12, detail = "Selecting heatmap features")
      top_for_heatmap <- plot_df %>% dplyr::filter(!is.na(q_value)) %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc)))
      if (isTRUE(input$identified_only_heatmap)) top_for_heatmap <- top_for_heatmap %>% dplyr::filter(!is.na(Label), Label != "", !grepl("^feature_", Label, ignore.case = TRUE))
      top_for_heatmap <- head(top_for_heatmap, input$heatmap_top_n)
      heatmap_features <- unique(top_for_heatmap$Feature_ID)

      hm_meta <- meta
      if (!isTRUE(input$keep_other_groups)) hm_meta <- hm_meta %>% dplyr::filter(group %in% c(input$group_a, input$group_b))
      hm_mat <- cd$matrix[heatmap_features, hm_meta$sample, drop = FALSE]
      rownames(hm_mat) <- top_for_heatmap$Label[match(rownames(hm_mat), top_for_heatmap$Feature_ID)]
      hm_mat <- normalize_for_heatmap(hm_mat)

      incProgress(0.10, detail = "Preparing output directory")
      outdir <- file.path(tempdir(), paste0("cd_modern_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(outdir, "barplots"), showWarnings = FALSE)
      resource_prefix <- paste0("results_", as.integer(Sys.time()))
      addResourcePath(resource_prefix, outdir)

      group_levels <- unique(meta$group)
      default_colors <- c(setNames(input$group_a_color, input$group_a), setNames(input$group_b_color, input$group_b), Blank = "#7b3294", QC = "#4c78a8")
      group_colors <- default_colors[names(default_colors) %in% group_levels]

      incProgress(0.12, detail = "Rendering PCA")
      ggsave(file.path(outdir, "pca.png"), plot_pca(cd$matrix, meta, group_colors), width = 7.5, height = 6.5, dpi = 220)

      incProgress(0.10, detail = "Rendering volcano")
      ggsave(file.path(outdir, "volcano.png"), plot_volcano(plot_df, input$group_a, input$group_b, input$group_a_color, input$group_b_color, input$fc_threshold, input$q_threshold, input$volcano_top_n, input$identified_only_labels), width = 9, height = 7.5, dpi = 220)

      incProgress(0.12, detail = "Rendering heatmap")
      save_heatmap(hm_mat, hm_meta, file.path(outdir, "heatmap.png"), group_colors)

      incProgress(0.10, detail = "Rendering sample summaries")
      tic_vals <- colSums(cd$matrix[, meta$sample, drop = FALSE], na.rm = TRUE)
      det_vals <- colSums(cd$matrix[, meta$sample, drop = FALSE] > 0, na.rm = TRUE)
      med_vals <- apply(cd$matrix[, meta$sample, drop = FALSE], 2, median, na.rm = TRUE)
      ggsave(file.path(outdir, "tic_like_sum.png"), plot_sample_summary(meta, tic_vals, "TIC-like Intensity Sum by Sample", "Peak Area Sum", group_colors), width = 9, height = 4.5, dpi = 220)
      ggsave(file.path(outdir, "detected_features.png"), plot_sample_summary(meta, det_vals, "Detected Features by Sample", "Detected Features", group_colors), width = 9, height = 4.5, dpi = 220)
      ggsave(file.path(outdir, "median_feature_intensity.png"), plot_sample_summary(meta, med_vals, "Median Feature Intensity by Sample", "Median Peak Area", group_colors), width = 9, height = 4.5, dpi = 220)

      incProgress(0.10, detail = "Rendering metabolite bar plots")
      top_bar <- plot_df %>% dplyr::filter(q_value <= input$q_threshold) %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc))) %>% head(input$bar_top_n)
      if (nrow(top_bar) > 0) {
        for (i in seq_len(nrow(top_bar))) {
          feature_id <- top_bar$Feature_ID[i]
          label <- top_bar$Label[i]
          safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", label)
          plot_individual_bar(feature_id, label, cd$matrix, meta, input$group_a, input$group_b, input$group_a_color, input$group_b_color, file.path(outdir, "barplots", paste0(sprintf("%03d", i), "_", safe_name, ".png")))
        }
      }

      incProgress(0.12, detail = "Saving tables")
      write.csv(plot_df, file.path(outdir, "annotated_feature_table_full.csv"), row.names = FALSE)
      write.csv(top_for_heatmap, file.path(outdir, "heatmap_features.csv"), row.names = FALSE)
      write.csv(top_bar, file.path(outdir, "individual_barplot_features.csv"), row.names = FALSE)

      rv$cd <- cd; rv$meta <- meta; rv$stats <- stats_df; rv$plot_df <- plot_df; rv$outdir <- outdir; rv$resource_prefix <- resource_prefix
      incProgress(0.02, detail = "Complete")
      list(outdir = outdir)
    })
  })

  output$pca_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("PCA", sprintf("/%s/pca.png", rv$resource_prefix)) })
  output$volcano_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Volcano Plot", sprintf("/%s/volcano.png", rv$resource_prefix)) })
  output$heatmap_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Clustered Heatmap", sprintf("/%s/heatmap.png", rv$resource_prefix), "Row z-scores with hierarchical clustering and class legend.") })
  output$tic_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("TIC-like Intensity Sum by Sample", sprintf("/%s/tic_like_sum.png", rv$resource_prefix)) })
  output$detect_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Detected Features by Sample", sprintf("/%s/detected_features.png", rv$resource_prefix)) })
  output$median_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Median Feature Intensity by Sample", sprintf("/%s/median_feature_intensity.png", rv$resource_prefix)) })

  output$stats_table <- renderDT({
    req(analysis())
    tbl <- rv$plot_df %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc))) %>% dplyr::select(any_of(c("Label", "Feature_ID", "mean_a", "mean_b", "log2fc", "p_value", "q_value", "Formula", "m/z", "RT [min]")))
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20))
  })

  output$barplot_gallery <- renderUI({
    req(analysis(), rv$resource_prefix)
    files <- list.files(file.path(rv$outdir, "barplots"), full.names = FALSE)
    if (length(files) == 0) return(tags$p("No individual bar plots generated."))
    files <- files[seq_len(min(length(files), input$bar_top_n))]
    tags$div(class = "barplot-grid",
      lapply(files, function(f) tags$div(class = "barplot-item",
        tags$img(src = sprintf("/%s/barplots/%s", rv$resource_prefix, utils::URLencode(f, reserved = TRUE)))
      ))
    )
  })

  output$download_zip <- downloadHandler(
    filename = function() paste0("cd_modern_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      req(rv$outdir)
      oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd(rv$outdir)
      files <- list.files(".", recursive = TRUE, all.files = FALSE)
      utils::zip(zipfile = file, files = files)
    }
  )
}

shinyApp(ui, server)

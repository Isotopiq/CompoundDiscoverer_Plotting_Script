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
  library(pROC)
  library(randomForest)
  library(e1071)
  library(pls)
  library(httr)
  library(jsonlite)
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
  else out <- paste0("feature_", seq_len(nrow(df)))
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
  list(raw = df, annotation = annotation_df, matrix = as.matrix(mat),
       sample_map = data.frame(original = area_cols, sample = sample_names, stringsAsFactors = FALSE),
       sheet = target_sheet, sheets = sheets)
}

infer_groups <- function(samples, group_a_name, group_b_name, group_a_pattern = "", group_b_pattern = "") {
  grp <- rep(NA_character_, length(samples))
  if (nzchar(group_a_pattern)) grp[stringr::str_detect(samples, regex(group_a_pattern, ignore_case = TRUE))] <- group_a_name
  if (nzchar(group_b_pattern)) grp[stringr::str_detect(samples, regex(group_b_pattern, ignore_case = TRUE))] <- group_b_name
  data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
}

filter_matrix <- function(mat, min_presence = 0.5) {
  keep <- rowMeans(mat > 0, na.rm = TRUE) >= min_presence
  mat[keep, , drop = FALSE]
}

impute_matrix <- function(mat, method = "halfmin") {
  x <- mat
  if (method == "zero") {
    x[is.na(x)] <- 0
  } else if (method == "median") {
    for (i in seq_len(nrow(x))) {
      med <- median(x[i, x[i, ] > 0], na.rm = TRUE)
      if (!is.finite(med)) med <- 1
      x[i, is.na(x[i, ]) | x[i, ] <= 0] <- med
    }
  } else {
    pos <- x[x > 0 & is.finite(x)]
    hm <- if (length(pos) > 0) min(pos, na.rm = TRUE) / 2 else 1
    x[is.na(x) | x <= 0] <- hm
  }
  x
}

normalize_matrix <- function(mat, method = "none") {
  x <- mat
  if (method == "sum") {
    cs <- colSums(x, na.rm = TRUE)
    cs[cs == 0] <- 1
    x <- sweep(x, 2, cs, "/") * median(cs, na.rm = TRUE)
  } else if (method == "median") {
    cm <- apply(x, 2, median, na.rm = TRUE)
    cm[cm == 0] <- 1
    x <- sweep(x, 2, cm, "/") * median(cm, na.rm = TRUE)
  }
  x
}

transform_scale_matrix <- function(mat, transform = "log2", scale_method = "auto") {
  x <- mat
  if (transform == "log2") x <- log2(x + 1)
  else if (transform == "log10") x <- log10(x + 1)
  if (scale_method == "auto") x <- t(scale(t(x)))
  else if (scale_method == "center") x <- sweep(x, 1, rowMeans(x, na.rm = TRUE), "-")
  else if (scale_method == "pareto") {
    x <- sweep(x, 1, rowMeans(x, na.rm = TRUE), "-")
    rsd <- apply(x, 1, sd, na.rm = TRUE)
    rsd[!is.finite(rsd) | rsd == 0] <- 1
    x <- sweep(x, 1, sqrt(rsd), "/")
  }
  x[is.na(x)] <- 0
  x
}

normalize_for_heatmap <- function(mat) {
  x <- t(scale(t(log2(mat + 1))))
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
    rows[[i]] <- data.frame(mean_a = mean_a, mean_b = mean_b,
                            log2fc = log2((mean_b + pseudocount) / (mean_a + pseudocount)),
                            p_value = welch_p(x, y), stringsAsFactors = FALSE)
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

plot_pca <- function(proc_mat, meta, group_colors) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  m <- t(proc_mat[, shared, drop = FALSE])
  pca <- prcomp(m, center = TRUE, scale. = FALSE)
  df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = meta$group[match(rownames(pca$x), meta$sample)], stringsAsFactors = FALSE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  ggplot(df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 4) +
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.08, color = NA, show.legend = FALSE) +
    scale_color_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 12) +
    labs(title = "PCA", x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2]), color = "Group") +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
}

plot_plsda <- function(proc_mat, meta, group_colors) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  Y <- model.matrix(~ y - 1)
  fit <- pls::plsr(Y ~ X, ncomp = 2, method = "oscorespls")
  scores <- data.frame(Comp1 = fit$scores[, 1], Comp2 = fit$scores[, 2], group = y)
  ggplot(scores, aes(x = Comp1, y = Comp2, color = group)) +
    geom_point(size = 4) +
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.08, color = NA, show.legend = FALSE) +
    scale_color_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 12) + labs(title = "PLS-DA", x = "Component 1", y = "Component 2") +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
}

plot_volcano <- function(plot_df, group_a, group_b, group_a_color, group_b_color, fc_threshold = 1, q_threshold = 0.05, label_top_n = 15, identified_only_labels = TRUE) {
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
    ggrepel::geom_text_repel(data = label_df, aes(label = Label), size = 3.2, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.2, segment.alpha = 0.6, seed = 123) +
    scale_color_manual(values = c("Not sig" = "grey70", setNames(group_a_color, paste0(group_a, " higher")), setNames(group_b_color, paste0(group_b, " higher")))) +
    theme_bw(base_size = 12) +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = expression(-log[10]("q-value")), color = NULL) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
}

plot_sample_summary <- function(meta, values, title, ylab, group_colors) {
  df <- data.frame(sample = names(values), value = as.numeric(values), stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, meta, by = "sample")
  ggplot(df, aes(x = sample, y = value, fill = group)) +
    geom_col() + scale_fill_manual(values = group_colors) + theme_bw(base_size = 12) +
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
  pheatmap::pheatmap(heatmap_mat,
    color = colorRampPalette(c("#6BAED6", "#F7F7F7", "#CB181D"))(100),
    cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
    annotation_col = ann_col, annotation_colors = ann_colors, annotation_legend = TRUE, legend = TRUE,
    border_color = "grey65", fontsize = 10, fontsize_row = 9, fontsize_col = 10,
    angle_col = 270, treeheight_row = 55, treeheight_col = 55, cellwidth = 28, cellheight = 15,
    filename = out_file, width = width_in, height = height_in)
}

plot_individual_bar <- function(feature_name, label_name, mat, meta, group_a, group_b, group_a_color, group_b_color, out_file) {
  samples_a <- meta %>% dplyr::filter(group == group_a, sample %in% colnames(mat)) %>% dplyr::pull(sample)
  samples_b <- meta %>% dplyr::filter(group == group_b, sample %in% colnames(mat)) %>% dplyr::pull(sample)
  vals_a <- as.numeric(mat[feature_name, samples_a]); vals_b <- as.numeric(mat[feature_name, samples_b])
  df <- data.frame(value = c(vals_a, vals_b), group = c(rep(group_a, length(vals_a)), rep(group_b, length(vals_b))), sample = c(samples_a, samples_b), stringsAsFactors = FALSE)
  pval <- welch_p(vals_a, vals_b); l2fc <- log2((mean(vals_b) + 1) / (mean(vals_a) + 1))
  ymax <- max(df$value, na.rm = TRUE); yline <- ymax * 1.10; ystar <- ymax * 1.16
  fill_vals <- c(group_a_color, group_b_color); names(fill_vals) <- c(group_a, group_b)
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.7) +
    geom_jitter(aes(color = group), width = 0.08, size = 3, alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = fill_vals) + scale_color_manual(values = fill_vals) +
    theme_bw(base_size = 13) +
    labs(title = label_name, subtitle = sprintf("log2FC = %.3f | p = %.3g", l2fc, pval), x = NULL, y = "Peak Area") +
    theme(plot.title = element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 8)), plot.subtitle = element_text(size = 11, margin = ggplot2::margin(b = 12)), panel.grid.minor = element_blank(), legend.position = "none") +
    annotate("segment", x = 1, xend = 2, y = yline, yend = yline, linewidth = 0.8) +
    annotate("segment", x = 1, xend = 1, y = yline * 0.98, yend = yline, linewidth = 0.8) +
    annotate("segment", x = 2, xend = 2, y = yline * 0.98, yend = yline, linewidth = 0.8) +
    annotate("text", x = 1.5, y = ystar, label = sig_stars(pval), size = 7, fontface = "bold") +
    expand_limits(y = ystar * 1.05)
  ggsave(out_file, p, width = 6.2, height = 6.0, dpi = 220)
}

compute_roc_table <- function(proc_mat, meta, group_a, group_b, top_n = 20) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  keep <- y %in% c(group_a, group_b)
  X <- X[keep, , drop = FALSE]; y <- droplevels(y[keep])
  if (nlevels(y) != 2) return(data.frame())
  rows <- lapply(colnames(X), function(f) {
    r <- tryCatch(pROC::roc(response = y, predictor = X[, f], quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
    data.frame(Feature_ID = f, AUC = if (is.null(r)) NA_real_ else as.numeric(pROC::auc(r)), stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(rows) %>% dplyr::arrange(dplyr::desc(AUC)) %>% head(top_n)
}

plot_single_roc <- function(proc_mat, meta, feature_id, group_a, group_b) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  keep <- y %in% c(group_a, group_b)
  X <- X[keep, , drop = FALSE]; y <- droplevels(y[keep])
  r <- pROC::roc(response = y, predictor = X[, feature_id], quiet = TRUE, levels = levels(y), direction = "<")
  aucv <- as.numeric(pROC::auc(r))
  df <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
  ggplot(df, aes(x = fpr, y = tpr)) + geom_path(color = "#2563eb", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
    coord_equal() + theme_bw(base_size = 12) +
    labs(title = paste0("ROC: ", feature_id), subtitle = sprintf("AUC = %.3f", aucv), x = "False Positive Rate", y = "True Positive Rate") +
    theme(plot.title = element_text(face = "bold"))
}

plot_multivariate_roc <- function(proc_mat, meta, group_a, group_b, top_features) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[top_features, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  keep <- y %in% c(group_a, group_b)
  X <- as.data.frame(X[keep, , drop = FALSE]); y <- droplevels(y[keep])
  y_num <- ifelse(y == levels(y)[2], 1, 0)
  mod <- glm(y_num ~ ., data = X, family = binomial())
  pr <- predict(mod, type = "response")
  r <- pROC::roc(response = y, predictor = pr, quiet = TRUE, levels = levels(y), direction = "<")
  aucv <- as.numeric(pROC::auc(r))
  df <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
  ggplot(df, aes(x = fpr, y = tpr)) + geom_path(color = "#7c3aed", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") + coord_equal() + theme_bw(base_size = 12) +
    labs(title = "Multivariate ROC", subtitle = sprintf("AUC = %.3f", aucv), x = "False Positive Rate", y = "True Positive Rate") +
    theme(plot.title = element_text(face = "bold"))
}

plot_rf_importance <- function(proc_mat, meta, group_a, group_b, top_n = 20) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  keep <- y %in% c(group_a, group_b)
  X <- X[keep, , drop = FALSE]; y <- droplevels(y[keep])
  rf <- randomForest::randomForest(x = X, y = y, importance = TRUE, ntree = 500)
  imp <- as.data.frame(randomForest::importance(rf)); imp$Feature_ID <- rownames(imp)
  colnm <- if ("MeanDecreaseAccuracy" %in% names(imp)) "MeanDecreaseAccuracy" else names(imp)[1]
  imp <- imp %>% dplyr::arrange(dplyr::desc(.data[[colnm]])) %>% head(top_n)
  ggplot(imp, aes(x = reorder(Feature_ID, .data[[colnm]]), y = .data[[colnm]])) + geom_col(fill = "#2563eb") + coord_flip() +
    theme_bw(base_size = 12) + labs(title = "Random Forest Feature Importance", x = NULL, y = colnm) + theme(plot.title = element_text(face = "bold"))
}

compute_svm_summary <- function(proc_mat, meta, group_a, group_b) {
  shared <- intersect(colnames(proc_mat), meta$sample)
  X <- t(proc_mat[, shared, drop = FALSE])
  y <- factor(meta$group[match(rownames(X), meta$sample)])
  keep <- y %in% c(group_a, group_b)
  X <- X[keep, , drop = FALSE]; y <- droplevels(y[keep])
  if (nrow(X) < 6) return(data.frame(metric = c("train_n", "test_n", "accuracy"), value = c(NA, NA, NA)))
  set.seed(123)
  idx <- sample(seq_len(nrow(X)), size = max(2, floor(0.7 * nrow(X))))
  mod <- e1071::svm(x = X[idx, , drop = FALSE], y = y[idx], kernel = "linear")
  pred <- predict(mod, X[-idx, , drop = FALSE])
  acc <- mean(pred == y[-idx])
  data.frame(metric = c("train_n", "test_n", "accuracy"), value = c(length(idx), nrow(X) - length(idx), round(acc, 4)))
}

simple_pathway_library <- function() {
  data.frame(
    pathway = c("Glycolysis / Gluconeogenesis", "TCA cycle", "Pentose phosphate pathway", "Amino acid metabolism", "Purine metabolism", "Pyrimidine metabolism", "Glutathione metabolism", "Fatty acid metabolism", "Methionine / one-carbon metabolism", "Lipid headgroup metabolism"),
    metabolite = c(
      "glucose;glucose-6-phosphate;fructose 6-phosphate;fructose 1,6-bisphosphate;glyceraldehyde 3-phosphate;pyruvate;lactate",
      "citrate;isocitrate;alpha-ketoglutarate;succinate;fumarate;malate;oxaloacetate",
      "ribose 5-phosphate;ribulose 5-phosphate;sedoheptulose 7-phosphate;erythrose 4-phosphate;glucose-6-phosphate",
      "glutamine;glutamate;serine;glycine;alanine;aspartate;asparagine;valine;leucine;isoleucine",
      "adenine;adenosine;amp;adp;atp;inosine;hypoxanthine;xanthine;guanine;gmp;gdp;gtp",
      "uracil;uridine;ump;udp;utp;cytidine;cmp;cdp;ctp;thymidine;tmp;tdp;ttp",
      "glutathione;gssg;cysteine;glycine;glutamate",
      "palmitate;oleate;acetyl-coa;acylcarnitine;fatty acid",
      "methionine;sam;sah;homocysteine;cystathionine;folate;serine;glycine;betaine;choline",
      "choline;phosphocholine;glycerophosphocholine;ethanolamine;phosphoethanolamine;cdp-choline"
    ), stringsAsFactors = FALSE)
}

run_simple_enrichment <- function(feature_labels, significant_labels) {
  lib <- simple_pathway_library()
  background <- unique(tolower(feature_labels))
  sig <- unique(tolower(significant_labels))
  bg_n <- length(background); sig_n <- length(sig)
  rows <- vector("list", nrow(lib))
  for (i in seq_len(nrow(lib))) {
    members <- unique(trimws(unlist(strsplit(tolower(lib$metabolite[i]), ";"))))
    members <- members[members %in% background]
    overlap <- intersect(sig, members)
    k <- length(overlap); m <- length(members); n <- bg_n - m
    pval <- if (sig_n > 0 && m > 0 && k > 0) phyper(q = k - 1, m = m, n = n, k = sig_n, lower.tail = FALSE) else 1
    rows[[i]] <- data.frame(pathway = lib$pathway[i], pathway_size = m, overlap = k,
                            hits = paste(overlap, collapse = "; "), p_value = pval,
                            enrichment_ratio = ifelse(m > 0 && sig_n > 0, (k / sig_n) / (m / bg_n), NA_real_), stringsAsFactors = FALSE)
  }
  out <- dplyr::bind_rows(rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out %>% dplyr::arrange(q_value, dplyr::desc(overlap))
}

plot_enrichment <- function(tbl, title = "Pathway / Set Enrichment") {
  d <- tbl %>% dplyr::filter(overlap > 0) %>% head(20)
  if (nrow(d) == 0) return(ggplot() + annotate("text", x = 0, y = 0, label = "No enriched pathways detected", size = 6) + theme_void())
  d$pathway <- factor(d$pathway, levels = rev(d$pathway))
  ggplot(d, aes(x = enrichment_ratio, y = pathway, size = overlap, color = -log10(pmax(q_value, 1e-300)))) +
    geom_point(alpha = 0.9) + scale_color_gradient(low = "#93c5fd", high = "#dc2626") + theme_bw(base_size = 12) +
    labs(title = title, x = "Enrichment ratio", y = NULL, color = "-log10(q)", size = "Overlap") + theme(plot.title = element_text(face = "bold"))
}


kegg_find_compound <- function(query) {
  url <- paste0("https://rest.kegg.jp/find/compound/", utils::URLencode(query, reserved = TRUE))
  res <- tryCatch(httr::GET(url, httr::timeout(20)), error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200) return(data.frame())
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  if (!nzchar(txt)) return(data.frame())
  lines <- unlist(strsplit(txt, "\n"))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(data.frame())
  parsed <- do.call(rbind, lapply(lines, function(x) {
    parts <- unlist(strsplit(x, "\t"))
    if (length(parts) < 2) return(NULL)
    data.frame(kegg_id = parts[1], kegg_name = parts[2], stringsAsFactors = FALSE)
  }))
  if (is.null(parsed)) return(data.frame())
  parsed
}

kegg_link_pathways <- function(kegg_ids) {
  kegg_ids <- unique(kegg_ids[nzchar(kegg_ids)])
  if (length(kegg_ids) == 0) return(data.frame())
  rows <- vector("list", length(kegg_ids))
  for (i in seq_along(kegg_ids)) {
    kid <- kegg_ids[i]
    url <- paste0("https://rest.kegg.jp/link/pathway/", kid)
    res <- tryCatch(httr::GET(url, httr::timeout(20)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) next
    txt <- httr::content(res, as = "text", encoding = "UTF-8")
    if (!nzchar(txt)) next
    lines <- unlist(strsplit(txt, "\n"))
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) next
    rows[[i]] <- do.call(rbind, lapply(lines, function(x) {
      parts <- unlist(strsplit(x, "\t"))
      if (length(parts) < 2) return(NULL)
      data.frame(kegg_id = parts[1], pathway_id = parts[2], stringsAsFactors = FALSE)
    }))
  }
  dplyr::bind_rows(rows)
}

kegg_get_pathway_names <- function(pathway_ids) {
  pathway_ids <- unique(pathway_ids[nzchar(pathway_ids)])
  if (length(pathway_ids) == 0) return(data.frame())
  chunks <- split(pathway_ids, ceiling(seq_along(pathway_ids) / 10))
  rows <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    ids <- paste(chunks[[i]], collapse = "+")
    url <- paste0("https://rest.kegg.jp/list/", ids)
    res <- tryCatch(httr::GET(url, httr::timeout(20)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) next
    txt <- httr::content(res, as = "text", encoding = "UTF-8")
    if (!nzchar(txt)) next
    lines <- unlist(strsplit(txt, "\n"))
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) next
    rows[[i]] <- do.call(rbind, lapply(lines, function(x) {
      parts <- unlist(strsplit(x, "\t"))
      if (length(parts) < 2) return(NULL)
      data.frame(pathway_id = parts[1], pathway_name = parts[2], stringsAsFactors = FALSE)
    }))
  }
  dplyr::bind_rows(rows)
}

hmdb_search_metabolite <- function(query) {
  url <- paste0("https://hmdb.ca/unearth/q?query=", utils::URLencode(query, reserved = TRUE), "&searcher=metabolites")
  res <- tryCatch(httr::GET(url, httr::timeout(20)), error = function(e) NULL)
  if (is.null(res)) return(data.frame())
  final_url <- tryCatch(res$url, error = function(e) NA_character_)
  data.frame(query = query, hmdb_search_url = final_url, stringsAsFactors = FALSE)
}

run_online_kegg_mapping <- function(feature_labels, significant_labels) {
  bg <- unique(feature_labels[!is.na(feature_labels) & nzchar(feature_labels)])
  sig <- unique(significant_labels[!is.na(significant_labels) & nzchar(significant_labels)])

  map_rows <- vector("list", length(bg))
  for (i in seq_along(bg)) {
    hits <- kegg_find_compound(bg[i])
    if (nrow(hits) > 0) {
      hits <- hits[1, , drop = FALSE]
      hits$query <- bg[i]
      map_rows[[i]] <- hits
    }
  }
  mapping <- dplyr::bind_rows(map_rows)
  if (nrow(mapping) == 0) {
    return(list(mapping = data.frame(), pathways = data.frame(), hmdb = data.frame()))
  }

  links <- kegg_link_pathways(mapping$kegg_id)
  names_tbl <- kegg_get_pathway_names(unique(links$pathway_id))
  joined <- mapping %>% dplyr::left_join(links, by = "kegg_id") %>% dplyr::left_join(names_tbl, by = "pathway_id")
  joined$is_significant <- joined$query %in% sig

  ptab <- joined %>%
    dplyr::filter(!is.na(pathway_id)) %>%
    dplyr::group_by(pathway_id, pathway_name) %>%
    dplyr::summarise(
      overlap = sum(is_significant, na.rm = TRUE),
      mapped = dplyr::n_distinct(query),
      sig_hits = paste(unique(query[is_significant]), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::filter(mapped > 0)

  if (nrow(ptab) > 0) {
    bg_n <- length(unique(mapping$query))
    sig_n <- length(intersect(unique(mapping$query), sig))
    ptab <- ptab %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p_value = ifelse(overlap > 0 && sig_n > 0, phyper(overlap - 1, mapped, bg_n - mapped, sig_n, lower.tail = FALSE), 1),
        enrichment_ratio = ifelse(sig_n > 0, (overlap / sig_n) / (mapped / bg_n), NA_real_)
      ) %>%
      dplyr::ungroup()
    ptab$q_value <- p.adjust(ptab$p_value, method = "BH")
    ptab <- ptab %>% dplyr::arrange(q_value, dplyr::desc(overlap))
  }

  hmdb_tbl <- dplyr::bind_rows(lapply(sig, hmdb_search_metabolite))
  list(mapping = joined, pathways = ptab, hmdb = hmdb_tbl)
}

img_panel <- function(title, src, note = NULL) {
  tags$div(class = "card-shell", tags$div(class = "card-head", title), if (!is.null(note)) tags$p(class = "card-note", note), tags$img(src = src, class = "plot-image"))
}

ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$script(src = "https://cdn.tailwindcss.com"),
    tags$style(HTML("body{background:#f8fafc;color:#0f172a;font-family:Inter,ui-sans-serif,system-ui}.app-shell{max-width:1800px;margin:0 auto;padding:24px}.hero{background:linear-gradient(180deg,#fff 0%,#f8fafc 100%);border:1px solid #e2e8f0;border-radius:22px;padding:24px 28px;box-shadow:0 10px 30px rgba(15,23,42,.04);margin-bottom:20px}.hero-title{font-size:28px;font-weight:800;letter-spacing:-.02em;margin:0}.hero-sub{margin-top:8px;color:#475569}.layout-grid{display:grid;grid-template-columns:370px minmax(0,1fr);gap:20px;align-items:start}.side-panel{position:sticky;top:20px;background:#fff;border:1px solid #e2e8f0;border-radius:20px;padding:20px;box-shadow:0 10px 30px rgba(15,23,42,.04);max-height:calc(100vh - 40px);overflow:auto}.content-panel{min-width:0}.section-title{font-size:13px;font-weight:700;text-transform:uppercase;letter-spacing:.08em;color:#64748b;margin-bottom:12px}.form-control,.selectize-input{border-radius:12px!important;border:1px solid #cbd5e1!important;box-shadow:none!important;min-height:44px!important;padding-top:10px!important;padding-bottom:10px!important}.btn-primary{width:100%;border:none;border-radius:14px;background:#2563eb;padding:12px 16px;font-weight:700;color:#fff}.btn-default{width:100%;border-radius:14px}.shiny-input-container{width:100%!important;margin-bottom:14px!important}.nav-pills{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:16px}.nav-pills>li{float:none!important}.nav-pills>li>a{border-radius:999px;padding:10px 16px;border:1px solid #dbeafe;background:#eff6ff;color:#1d4ed8;font-weight:600}.nav-pills>li.active>a,.nav-pills>li.active>a:hover,.nav-pills>li.active>a:focus{background:#2563eb;border-color:#2563eb;color:#fff}.card-shell{background:#fff;border:1px solid #e2e8f0;border-radius:20px;padding:18px;box-shadow:0 10px 30px rgba(15,23,42,.04);margin-bottom:18px;overflow:hidden}.card-head{font-size:18px;font-weight:700;margin-bottom:8px}.card-note{color:#64748b;margin-bottom:12px}.plot-image{width:100%;height:auto;display:block;border-radius:14px;border:1px solid #e5e7eb;background:#fff}.barplot-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px}.barplot-item{background:#fff;border:1px solid #e2e8f0;border-radius:18px;padding:12px;box-shadow:0 10px 30px rgba(15,23,42,.04)}.barplot-item img{width:100%;height:auto;display:block;border-radius:12px}.dataTables_wrapper{background:#fff;border:1px solid #e2e8f0;border-radius:20px;padding:14px;box-shadow:0 10px 30px rgba(15,23,42,.04)}.progress{height:12px;border-radius:999px;overflow:hidden;background:#e2e8f0}.progress-bar{background:#2563eb}@media (max-width:1150px){.layout-grid{grid-template-columns:1fr}.side-panel{position:static;max-height:none}}"))
  ),
  div(class = "app-shell",
    div(class = "hero", h1(class = "hero-title", "Compound Discoverer Metabolomics Explorer+"), p(class = "hero-sub", "Modern light dashboard with preprocessing, PCA, volcano plots, clustered heatmaps, biomarker ROC, supervised modeling, pathway analysis, enrichment, and individual metabolite statistics.")),
    div(class = "layout-grid",
      div(class = "side-panel",
        div(class = "section-title", "Input"),
        fileInput("xlsx", "Upload Compound Discoverer Excel file", accept = c(".xlsx", ".xls")),
        uiOutput("sheet_ui"),
        fileInput("metadata", "Optional metadata CSV", accept = ".csv"),
        textInput("group_a", "Group A name", value = "Control"),
        textInput("group_b", "Group B name", value = "KO"),
        textInput("pattern_a", "Auto-assign regex for Group A", value = "yl_001|yl_002|yl_003"),
        textInput("pattern_b", "Auto-assign regex for Group B", value = "yl_004|yl_005|yl_006|yl_007|yl_008|yl_009|yl_010"),
        checkboxInput("keep_other_groups", "Keep other groups in PCA/heatmap if present", FALSE),
        tags$hr(), div(class = "section-title", "Preprocessing"),
        sliderInput("min_presence", "Minimum feature presence", min = 0, max = 1, value = 0.5, step = 0.05),
        selectInput("impute_method", "Imputation", choices = c("halfmin", "zero", "median"), selected = "halfmin"),
        selectInput("norm_method", "Normalization", choices = c("none", "sum", "median"), selected = "none"),
        selectInput("transform_method", "Transformation", choices = c("log2", "log10", "none"), selected = "log2"),
        selectInput("scale_method", "Scaling", choices = c("auto", "pareto", "center", "none"), selected = "auto"),
        tags$hr(), div(class = "section-title", "Thresholds & display"),
        checkboxInput("identified_only_heatmap", "Heatmap: identified compounds only", TRUE),
        checkboxInput("identified_only_labels", "Volcano labels: identified compounds only", TRUE),
        numericInput("fc_threshold", "Volcano |log2FC| threshold", value = 1, min = 0, step = 0.1),
        numericInput("q_threshold", "Volcano q-value threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("volcano_top_n", "Volcano labels top N", value = 15, min = 1, step = 1),
        numericInput("heatmap_top_n", "Heatmap top N features", value = 40, min = 5, step = 1),
        numericInput("bar_top_n", "Individual bar plots top N", value = 20, min = 1, step = 1),
        numericInput("roc_top_n", "Top biomarker ROC features", value = 15, min = 3, step = 1),
              tags$hr(),
              div(class = "section-title", "Online identifier mapping"),
              checkboxInput("enable_online_mapping", "Enable KEGG/HMDB online lookup", FALSE),
              checkboxInput("online_kegg_only_significant", "Use significant metabolites for online focus", TRUE),
        tags$hr(), div(class = "section-title", "Colors"),
        textInput("group_a_color", "Group A color", value = "#E91E63"),
        textInput("group_b_color", "Group B color", value = "#43B649"),
        actionButton("run", "Run Analysis", class = "btn btn-primary"), br(), br(), downloadButton("download_zip", "Download Results ZIP", class = "btn btn-default")
      ),
      div(class = "content-panel",
        tabsetPanel(type = "pills",
          tabPanel("Overview", div(class = "card-shell", div(class = "card-head", "Detected samples"), DTOutput("sample_table")), div(class = "card-shell", div(class = "card-head", "Metadata used"), DTOutput("meta_table"))),
          tabPanel("PCA", uiOutput("pca_ui")),
          tabPanel("Volcano", uiOutput("volcano_ui")),
          tabPanel("Heatmap", uiOutput("heatmap_ui")),
          tabPanel("Summaries", uiOutput("tic_ui"), uiOutput("detect_ui"), uiOutput("median_ui")),
          tabPanel("ROC", uiOutput("roc_single_ui"), uiOutput("roc_multi_ui"), DTOutput("roc_table")),
          tabPanel("Modeling", uiOutput("plsda_ui"), uiOutput("rf_ui"), div(class = "card-shell", div(class = "card-head", "SVM summary"), DTOutput("svm_table"))),
          tabPanel("Pathway", uiOutput("pathway_ui"), DTOutput("pathway_table")),
          tabPanel("MSEA", uiOutput("msea_ui"), DTOutput("msea_table")),
                tabPanel("Online KEGG/HMDB", uiOutput("online_kegg_ui"), DTOutput("online_kegg_table"), DTOutput("online_hmdb_table")),
          tabPanel("Top table", DTOutput("stats_table")),
          tabPanel("Bar plots", uiOutput("barplot_gallery"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(cd = NULL, meta = NULL, outdir = NULL, resource_prefix = NULL, plot_df = NULL, proc_mat = NULL, pathway_tbl = NULL, roc_tbl = NULL, top_bar = NULL, svm_tbl = NULL, online_kegg_tbl = NULL, online_hmdb_tbl = NULL)

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

  observe({ cd <- try(parsed_data(), silent = TRUE); if (!inherits(cd, "try-error")) output$sample_table <- renderDT(datatable(cd$sample_map, options = list(scrollX = TRUE, pageLength = 12))) })
  observe({ md <- try(metadata_used(), silent = TRUE); if (!inherits(md, "try-error")) output$meta_table <- renderDT(datatable(md, options = list(scrollX = TRUE, pageLength = 12))) })

  analysis <- eventReactive(input$run, {
    withProgress(message = "Running analysis", value = 0, {
      cd <- parsed_data(); meta <- metadata_used()
      incProgress(0.06, detail = "Validating metadata")
      meta <- meta %>% dplyr::filter(sample %in% colnames(cd$matrix))
      compare_meta <- meta %>% dplyr::filter(group %in% c(input$group_a, input$group_b))
      incProgress(0.10, detail = "Preprocessing feature matrix")
      filt <- filter_matrix(cd$matrix, input$min_presence)
      imp <- impute_matrix(filt, input$impute_method)
      norm <- normalize_matrix(imp, input$norm_method)
      proc <- transform_scale_matrix(norm, input$transform_method, input$scale_method)
      incProgress(0.10, detail = "Computing differential statistics")
      stats_df <- compute_stats(norm, compare_meta, input$group_a, input$group_b)
      plot_df <- build_plot_table(stats_df, cd$annotation)
      incProgress(0.08, detail = "Selecting heatmap features")
      top_for_heatmap <- plot_df %>% dplyr::filter(!is.na(q_value)) %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc)))
      if (isTRUE(input$identified_only_heatmap)) top_for_heatmap <- top_for_heatmap %>% dplyr::filter(!is.na(Label), Label != "", !grepl("^feature_", Label, ignore.case = TRUE))
      top_for_heatmap <- head(top_for_heatmap, input$heatmap_top_n)
      heatmap_features <- unique(top_for_heatmap$Feature_ID)
      hm_meta <- if (isTRUE(input$keep_other_groups)) meta else compare_meta
      hm_norm <- norm[intersect(rownames(norm), heatmap_features), intersect(colnames(norm), hm_meta$sample), drop = FALSE]
      rownames(hm_norm) <- top_for_heatmap$Label[match(rownames(hm_norm), top_for_heatmap$Feature_ID)]
      hm_norm <- normalize_for_heatmap(hm_norm)
      incProgress(0.06, detail = "Preparing output directory")
      outdir <- file.path(tempdir(), paste0("cd_plus_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(outdir, "barplots"), showWarnings = FALSE)
      resource_prefix <- paste0("results_", as.integer(Sys.time()))
      addResourcePath(resource_prefix, outdir)
      group_levels <- unique(meta$group)
      default_colors <- c(setNames(input$group_a_color, input$group_a), setNames(input$group_b_color, input$group_b), Blank = "#7b3294", QC = "#4c78a8")
      group_colors <- default_colors[names(default_colors) %in% group_levels]
      pca_meta <- if (isTRUE(input$keep_other_groups)) meta else compare_meta
      pca_groups <- unique(pca_meta$group)
      pca_colors <- default_colors[names(default_colors) %in% pca_groups]
      cmp_colors <- default_colors[names(default_colors) %in% c(input$group_a, input$group_b)]
      incProgress(0.08, detail = "Rendering PCA / PLS-DA")
      ggsave(file.path(outdir, "pca.png"), plot_pca(proc, pca_meta, pca_colors), width = 7.5, height = 6.5, dpi = 220)
      ggsave(file.path(outdir, "plsda.png"), plot_plsda(proc, compare_meta, cmp_colors), width = 7.5, height = 6.5, dpi = 220)
      incProgress(0.06, detail = "Rendering volcano")
      ggsave(file.path(outdir, "volcano.png"), plot_volcano(plot_df, input$group_a, input$group_b, input$group_a_color, input$group_b_color, input$fc_threshold, input$q_threshold, input$volcano_top_n, input$identified_only_labels), width = 9, height = 7.5, dpi = 220)
      incProgress(0.06, detail = "Rendering heatmap")
      save_heatmap(hm_norm, hm_meta, file.path(outdir, "heatmap.png"), group_colors)
      incProgress(0.08, detail = "Rendering summaries")
      tic_vals <- colSums(norm[, meta$sample, drop = FALSE], na.rm = TRUE)
      det_vals <- colSums(norm[, meta$sample, drop = FALSE] > 0, na.rm = TRUE)
      med_vals <- apply(norm[, meta$sample, drop = FALSE], 2, median, na.rm = TRUE)
      ggsave(file.path(outdir, "tic_like_sum.png"), plot_sample_summary(meta, tic_vals, "TIC-like Intensity Sum by Sample", "Peak Area Sum", group_colors), width = 9, height = 4.5, dpi = 220)
      ggsave(file.path(outdir, "detected_features.png"), plot_sample_summary(meta, det_vals, "Detected Features by Sample", "Detected Features", group_colors), width = 9, height = 4.5, dpi = 220)
      ggsave(file.path(outdir, "median_feature_intensity.png"), plot_sample_summary(meta, med_vals, "Median Feature Intensity by Sample", "Median Peak Area", group_colors), width = 9, height = 4.5, dpi = 220)
      incProgress(0.10, detail = "Computing biomarker ROC and models")
      roc_tbl <- compute_roc_table(proc, compare_meta, input$group_a, input$group_b, input$roc_top_n)
      if (nrow(roc_tbl) > 0) ggsave(file.path(outdir, "roc_single.png"), plot_single_roc(proc, compare_meta, roc_tbl$Feature_ID[1], input$group_a, input$group_b), width = 7, height = 6, dpi = 220)
      top_roc_features <- head(roc_tbl$Feature_ID, 5)
      if (length(top_roc_features) >= 2) ggsave(file.path(outdir, "roc_multi.png"), plot_multivariate_roc(proc, compare_meta, input$group_a, input$group_b, top_roc_features), width = 7, height = 6, dpi = 220)
      ggsave(file.path(outdir, "rf_importance.png"), plot_rf_importance(proc, compare_meta, input$group_a, input$group_b, top_n = 20), width = 8, height = 7, dpi = 220)
      svm_tbl <- compute_svm_summary(proc, compare_meta, input$group_a, input$group_b)
      incProgress(0.08, detail = "Running pathway analysis and MSEA")
      sig_labels <- plot_df %>% dplyr::filter(q_value <= input$q_threshold, abs(log2fc) >= input$fc_threshold) %>% dplyr::pull(Label)
      bg_labels <- unique(plot_df$Label)
      pathway_tbl <- run_simple_enrichment(bg_labels, sig_labels)
      msea_tbl <- pathway_tbl
      ggsave(file.path(outdir, "pathway.png"), plot_enrichment(pathway_tbl, "Pathway Analysis"), width = 8.5, height = 7, dpi = 220)
      ggsave(file.path(outdir, "msea.png"), plot_enrichment(msea_tbl, "Metabolite Set Enrichment"), width = 8.5, height = 7, dpi = 220)
      
      incProgress(0.06, detail = "Optional online KEGG/HMDB mapping")
      online_kegg_tbl <- data.frame()
      online_hmdb_tbl <- data.frame()
      if (isTRUE(input$enable_online_mapping)) {
        focus_labels <- if (isTRUE(input$online_kegg_only_significant)) sig_labels else bg_labels
        online_res <- run_online_kegg_mapping(bg_labels, focus_labels)
        online_kegg_tbl <- online_res$pathways
        online_hmdb_tbl <- online_res$hmdb
        if (nrow(online_kegg_tbl) > 0) {
          ggsave(file.path(outdir, "online_kegg.png"), plot_enrichment(online_kegg_tbl, "Online KEGG Pathway Mapping"), width = 8.5, height = 7, dpi = 220)
          write.csv(online_res$mapping, file.path(outdir, "online_kegg_mapping_details.csv"), row.names = FALSE)
          write.csv(online_kegg_tbl, file.path(outdir, "online_kegg_pathways.csv"), row.names = FALSE)
        }
        if (nrow(online_hmdb_tbl) > 0) {
          write.csv(online_hmdb_tbl, file.path(outdir, "online_hmdb_links.csv"), row.names = FALSE)
        }
      }

incProgress(0.10, detail = "Rendering metabolite bar plots")
      top_bar <- plot_df %>% dplyr::filter(q_value <= input$q_threshold) %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc))) %>% head(input$bar_top_n)
      if (nrow(top_bar) > 0) for (i in seq_len(nrow(top_bar))) {
        feature_id <- top_bar$Feature_ID[i]; label <- top_bar$Label[i]
        safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", label)
        plot_individual_bar(feature_id, label, norm, compare_meta, input$group_a, input$group_b, input$group_a_color, input$group_b_color, file.path(outdir, "barplots", paste0(sprintf("%03d", i), "_", safe_name, ".png")))
      }
      incProgress(0.08, detail = "Saving tables")
      write.csv(plot_df, file.path(outdir, "annotated_feature_table_full.csv"), row.names = FALSE)
      write.csv(top_for_heatmap, file.path(outdir, "heatmap_features.csv"), row.names = FALSE)
      write.csv(top_bar, file.path(outdir, "individual_barplot_features.csv"), row.names = FALSE)
      write.csv(roc_tbl, file.path(outdir, "roc_table.csv"), row.names = FALSE)
      write.csv(pathway_tbl, file.path(outdir, "pathway_analysis.csv"), row.names = FALSE)
      write.csv(msea_tbl, file.path(outdir, "msea.csv"), row.names = FALSE)
      write.csv(svm_tbl, file.path(outdir, "svm_summary.csv"), row.names = FALSE)
      rv$cd <- cd; rv$meta <- meta; rv$outdir <- outdir; rv$resource_prefix <- resource_prefix; rv$plot_df <- plot_df; rv$proc_mat <- proc; rv$pathway_tbl <- pathway_tbl; rv$roc_tbl <- roc_tbl; rv$top_bar <- top_bar; rv$svm_tbl <- svm_tbl
      incProgress(0.04, detail = "Complete")
      list(outdir = outdir)
    })
  })

  output$pca_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("PCA", sprintf("/%s/pca.png", rv$resource_prefix)) })
  output$volcano_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Volcano Plot", sprintf("/%s/volcano.png", rv$resource_prefix)) })
  output$heatmap_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Clustered Heatmap", sprintf("/%s/heatmap.png", rv$resource_prefix), "Hierarchical clustering with class annotation.") })
  output$tic_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("TIC-like Intensity Sum by Sample", sprintf("/%s/tic_like_sum.png", rv$resource_prefix)) })
  output$detect_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Detected Features by Sample", sprintf("/%s/detected_features.png", rv$resource_prefix)) })
  output$median_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Median Feature Intensity by Sample", sprintf("/%s/median_feature_intensity.png", rv$resource_prefix)) })
  output$roc_single_ui <- renderUI({ req(analysis(), rv$resource_prefix); if (!file.exists(file.path(rv$outdir, "roc_single.png"))) return(NULL); img_panel("Top Univariate ROC", sprintf("/%s/roc_single.png", rv$resource_prefix)) })
  output$roc_multi_ui <- renderUI({ req(analysis(), rv$resource_prefix); if (!file.exists(file.path(rv$outdir, "roc_multi.png"))) return(NULL); img_panel("Multivariate ROC", sprintf("/%s/roc_multi.png", rv$resource_prefix)) })
  output$plsda_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("PLS-DA", sprintf("/%s/plsda.png", rv$resource_prefix)) })
  output$rf_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Random Forest Importance", sprintf("/%s/rf_importance.png", rv$resource_prefix)) })
  output$pathway_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Pathway Analysis", sprintf("/%s/pathway.png", rv$resource_prefix), "Over-representation against a bundled metabolite pathway library.") })
  output$msea_ui <- renderUI({ req(analysis(), rv$resource_prefix); img_panel("Metabolite Set Enrichment", sprintf("/%s/msea.png", rv$resource_prefix), "Enrichment against bundled metabolite set categories.") })
  output$stats_table <- renderDT({ req(analysis()); tbl <- rv$plot_df %>% dplyr::arrange(q_value, dplyr::desc(abs(log2fc))) %>% dplyr::select(any_of(c("Label", "Feature_ID", "mean_a", "mean_b", "log2fc", "p_value", "q_value", "Formula", "m/z", "RT [min]"))); datatable(tbl, options = list(scrollX = TRUE, pageLength = 20)) })
  output$roc_table <- renderDT({ req(analysis()); datatable(rv$roc_tbl, options = list(scrollX = TRUE, pageLength = 15)) })
  output$pathway_table <- renderDT({ req(analysis()); datatable(rv$pathway_tbl, options = list(scrollX = TRUE, pageLength = 15)) })
  output$msea_table <- renderDT({ req(analysis()); datatable(rv$pathway_tbl, options = list(scrollX = TRUE, pageLength = 15)) })
  output$svm_table <- renderDT({ req(analysis()); datatable(rv$svm_tbl, options = list(dom = 't', pageLength = 10)) })
  output$barplot_gallery <- renderUI({
    req(analysis(), rv$resource_prefix)
    files <- list.files(file.path(rv$outdir, "barplots"), full.names = FALSE)
    if (length(files) == 0) return(tags$p("No individual bar plots generated."))
    files <- files[seq_len(min(length(files), input$bar_top_n))]
    tags$div(class = "barplot-grid", lapply(files, function(f) tags$div(class = "barplot-item", tags$img(src = sprintf("/%s/barplots/%s", rv$resource_prefix, utils::URLencode(f, reserved = TRUE))))))
  })

  output$online_kegg_ui <- renderUI({
    req(analysis())
    if (is.null(rv$online_kegg_tbl) || nrow(rv$online_kegg_tbl) == 0 || !file.exists(file.path(rv$outdir, "online_kegg.png"))) {
      return(div(class = "card-shell", div(class = "card-head", "Online KEGG/HMDB"), p(class = "card-note", "Enable online lookup in the sidebar and rerun analysis. Results require internet access from the deployed app.")))
    }
    img_panel("Online KEGG Pathway Mapping", sprintf("/%s/online_kegg.png", rv$resource_prefix), "Pathways mapped from live KEGG REST lookups. HMDB links are listed in the table below.")
  })
  output$online_kegg_table <- renderDT({
    req(analysis())
    dat <- rv$online_kegg_tbl
    if (is.null(dat) || nrow(dat) == 0) dat <- data.frame(message = "No online KEGG results available")
    datatable(dat, options = list(scrollX = TRUE, pageLength = 15))
  })
  output$online_hmdb_table <- renderDT({
    req(analysis())
    dat <- rv$online_hmdb_tbl
    if (is.null(dat) || nrow(dat) == 0) dat <- data.frame(message = "No online HMDB results available")
    datatable(dat, options = list(scrollX = TRUE, pageLength = 10))
  })


  output$download_zip <- downloadHandler(
    filename = function() paste0("cd_modern_plus_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      req(rv$outdir)
      oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd(rv$outdir)
      files <- list.files(".", recursive = TRUE, all.files = FALSE)
      utils::zip(zipfile = file, files = files)
    }
  )
}

shinyApp(ui, server)

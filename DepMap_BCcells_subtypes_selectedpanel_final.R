# ==========================================================
# DepMap (CCLE) — Polyamine genes across BREAST cancer subtypes
# FINAL (clean) + readable p-values on Tukey brackets (DepMap-style)
#
# - Subtype labels (compact): LumA, LumB, HER2, B-L
# - Colors: LumA=darkkhaki, LumB=lightpink, HER2=orange, B-L=lightblue
# - Bars = mean ± SD (alpha shading + black outline)
# - Dots = individual cell lines (BLACK; size 1.8; alpha 0.75; jitter width 0.18)
# - Per gene: one-way ANOVA across subtypes (facet strip, 2 lines)
#     line 1: GENE (bold, 15 pt)
#     line 2: ANOVA p=... (regular, 11 pt)
# - Stats: Tukey HSD per gene; only significant (p<0.05) shown as brackets
# - Bracket p-value text size controlled by pval_bracket_size (set to 4.5)
# - Export: TIFF 600 dpi, LZW to Downloads
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","ggpubr","readr","purrr","tibble","stringr","ggtext")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ---------------- SETTINGS ----------------
setwd("C:/Users/u0171728/Downloads")  # adjust if needed

expr_file <- "CCLE_expression_full.csv"
meta_file <- "sample_info.csv"

# Subtype order + colors (compact labels)
subtype_levels <- c("LumA","LumB","HER2","B-L")
subtype_cols <- c(
  "LumA" = "darkkhaki",
  "LumB" = "lightpink",
  "HER2" = "orange",
  "B-L"  = "lightblue"
)

# Gene groups (facet order left->right)
transport_genes <- c("ATP13A2","ATP13A3","ATP13A4")
synth_genes     <- c("ODC1","SRM","SMS")
catab_genes     <- c("SMOX","PAOX","SAT1")

# TIFF export settings
tiff_dpi  <- 600
tiff_w_in <- 10
tiff_h_in <- 7

# ---------------- STYLE (DepMap-like) ----------------
base_size    <- 17
axis_x_size  <- 16
axis_y_size  <- 16
y_title_sz   <- 16

dot_jit_w <- 0.18
dot_size  <- 1.8
dot_alpha <- 0.75

bar_width <- 0.65
bar_alpha <- 0.75

# --- NEW: readable bracket p-values + line thickness ---
pval_bracket_size <- 4.5   # p-value label size on brackets
bracket_size      <- 0.4   # bracket line thickness
tip_len           <- 0.01

# ---------------- SAVE TIFF ----------------
save_tiff <- function(plot_obj, filename_no_ext) {
  ggplot2::ggsave(
    filename = paste0(filename_no_ext, ".tiff"),
    plot = plot_obj,
    device = "tiff",
    dpi = tiff_dpi,
    width = tiff_w_in,
    height = tiff_h_in,
    units = "in",
    compression = "lzw"
  )
}

# ---------------- LOAD ----------------
expr <- readr::read_csv(expr_file, show_col_types = FALSE)
meta <- readr::read_csv(meta_file, show_col_types = FALSE)
colnames(expr)[1] <- "DepMap_ID"

# ---------------- MAP BREAST SUBTYPES ----------------
# Requires columns in sample_info.csv: lineage, lineage_sub_subtype, DepMap_ID, cell_line_name
meta_brca <- meta %>%
  dplyr::filter(lineage == "breast") %>%
  dplyr::mutate(
    Subtype = dplyr::case_when(
      grepl("HER2pos", lineage_sub_subtype, ignore.case = TRUE) & grepl("ERneg", lineage_sub_subtype, ignore.case = TRUE) ~ "HER2",
      grepl("HER2pos", lineage_sub_subtype, ignore.case = TRUE) & grepl("ERpos", lineage_sub_subtype, ignore.case = TRUE) ~ "LumB",
      grepl("HER2neg", lineage_sub_subtype, ignore.case = TRUE) & grepl("ERpos", lineage_sub_subtype, ignore.case = TRUE) ~ "LumA",
      grepl("HER2neg", lineage_sub_subtype, ignore.case = TRUE) & grepl("ERneg", lineage_sub_subtype, ignore.case = TRUE) ~ "B-L",
      TRUE ~ NA_character_
    ),
    # ---- MANUAL OVERRIDE(S) ----
    Subtype = dplyr::case_when(
      toupper(cell_line_name) %in% c("MDA-MB-453","MDAMB453","MDA MB 453","MDA_MB_453") ~ "HER2",
      TRUE ~ Subtype
    ),
    Subtype = factor(Subtype, levels = subtype_levels)
  ) %>%
  dplyr::filter(!is.na(Subtype)) %>%
  dplyr::distinct(DepMap_ID, cell_line_name, Subtype)

# ---------------- HELPERS ----------------
find_gene_col <- function(gene, cn) {
  hit <- grep(paste0("^", gene, " \\("), cn, value = TRUE)
  if (length(hit) == 0) NA_character_ else hit[1]
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

# Tukey comparison strings can include hyphens (e.g., "B-L") so split by last hyphen.
split_tukey_comparison <- function(x) {
  group2 <- sub("^.*-", "", x)
  group1 <- sub("-[^-]+$", "", x)
  tibble::tibble(group1 = group1, group2 = group2)
}

# ---------------- BUILD LONG (SUBTYPE PANELS) ----------------
build_long_for_genes <- function(gene_set) {
  
  gene_cols <- setNames(
    vapply(gene_set, find_gene_col, character(1), cn = colnames(expr)),
    gene_set
  )
  
  missing <- names(gene_cols)[is.na(gene_cols)]
  if (length(missing) > 0) {
    warning("Missing genes (skipped): ", paste(missing, collapse = ", "))
    gene_cols <- gene_cols[!is.na(gene_cols)]
  }
  stopifnot(length(gene_cols) > 0)
  
  expr %>%
    dplyr::select(DepMap_ID, dplyr::all_of(unname(gene_cols))) %>%
    dplyr::inner_join(meta_brca, by = "DepMap_ID") %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(unname(gene_cols)),
      names_to  = "gene_col",
      values_to = "expr_value"
    ) %>%
    dplyr::mutate(
      gene    = stringr::str_trim(sub(" \\(.*", "", gene_col)),
      gene    = factor(gene, levels = gene_set),
      Subtype = factor(Subtype, levels = subtype_levels)
    ) %>%
    dplyr::filter(!is.na(expr_value), !is.na(gene), !is.na(Subtype))
}

compute_stats <- function(df_long, gene_levels) {
  
  df_long <- df_long %>% dplyr::mutate(gene = factor(gene, levels = gene_levels))
  
  # ANOVA per gene (facet strip labels)
  overall <- purrr::map_dfr(gene_levels, function(g) {
    d <- df_long %>% dplyr::filter(gene == g)
    if (dplyr::n_distinct(d$Subtype) < 2) return(tibble::tibble(gene = g, p_anova = NA_real_))
    p <- tryCatch(summary(aov(expr_value ~ Subtype, data = d))[[1]][["Pr(>F)"]][1],
                  error = function(e) NA_real_)
    tibble::tibble(gene = g, p_anova = p)
  }) %>%
    dplyr::mutate(
      gene = factor(gene, levels = gene_levels),
      gene_label = paste0(
        "<b><span style='font-size:15pt;'>", as.character(gene), "</span></b>",
        "<br><span style='font-size:11pt;'>ANOVA p=", fmt_p(p_anova), "</span>"
      )
    )
  
  # Tukey HSD per gene (only significant)
  pairwise <- purrr::map_dfr(gene_levels, function(g) {
    d <- df_long %>% dplyr::filter(gene == g)
    if (dplyr::n_distinct(d$Subtype) < 2) return(tibble::tibble())
    
    tk <- tryCatch(
      TukeyHSD(aov(expr_value ~ Subtype, data = d), "Subtype")$Subtype %>% as.data.frame(),
      error = function(e) NULL
    )
    if (is.null(tk)) return(tibble::tibble())
    
    tk$comparison <- rownames(tk)
    
    out <- tibble::as_tibble(tk) %>%
      dplyr::mutate(comparison = tk$comparison) %>%
      dplyr::bind_cols(split_tukey_comparison(tk$comparison)) %>%
      dplyr::rename(p_adj = `p adj`) %>%
      dplyr::filter(!is.na(p_adj), p_adj < 0.05)
    
    if (nrow(out) == 0) return(tibble::tibble())
    
    ymax <- max(d$expr_value, na.rm = TRUE)
    ymin <- min(d$expr_value, na.rm = TRUE)
    step <- 0.10 * (ymax - ymin + 1e-6)
    
    out %>%
      dplyr::mutate(
        gene = factor(g, levels = gene_levels),
        y.position = ymax + step * dplyr::row_number(),
        label = paste0("p=", fmt_p(p_adj))
      ) %>%
      dplyr::select(gene, group1, group2, p_adj, label, y.position)
  })
  
  list(overall = overall, pairwise = pairwise, df_long = df_long)
}

plot_panel <- function(df_long, gene_levels, ncol = 3) {
  
  st <- compute_stats(df_long, gene_levels)
  df_long <- st$df_long
  
  # mean ± SD for bars/errorbars
  sum_df <- df_long %>%
    dplyr::group_by(gene, Subtype) %>%
    dplyr::summarise(
      mean = mean(expr_value, na.rm = TRUE),
      sd   = sd(expr_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene    = factor(gene, levels = gene_levels),
      Subtype = factor(Subtype, levels = subtype_levels)
    )
  
  lab_vec <- setNames(st$overall$gene_label, as.character(st$overall$gene))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = sum_df,
      ggplot2::aes(x = Subtype, y = mean, fill = Subtype),
      width = bar_width, color = "black", alpha = bar_alpha
    ) +
    ggplot2::geom_errorbar(
      data = sum_df,
      ggplot2::aes(x = Subtype, ymin = mean - sd, ymax = mean + sd),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = df_long,
      ggplot2::aes(x = Subtype, y = expr_value),
      color = "black",
      width = dot_jit_w,
      size  = dot_size,
      alpha = dot_alpha
    ) +
    ggplot2::scale_fill_manual(values = subtype_cols, drop = FALSE) +
    ggplot2::facet_wrap(
      ~ gene,
      ncol = ncol,
      scales = "fixed",
      drop = FALSE,
      labeller = ggplot2::labeller(gene = lab_vec)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "**mRNA expression**<br><span style='font-weight:400;'>mean \u00B1 SD; DepMap</span>"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text   = ggtext::element_markdown(),
      axis.text.x  = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_y_size),
      axis.title.y = ggtext::element_markdown(size = y_title_sz),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none"
    )
  
  # --- NEW: bracket p-value text size + bracket line thickness ---
  if (nrow(st$pairwise) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data         = st$pairwise,
        label        = "label",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = tip_len,
        hide.ns      = TRUE,
        size         = pval_bracket_size,
        bracket.size = bracket_size
      )
  }
  
  p
}

# ==========================================================
# RUN: SUBTYPE PANELS
# ==========================================================
long_atp   <- build_long_for_genes(transport_genes)
long_synth <- build_long_for_genes(synth_genes)
long_catab <- build_long_for_genes(catab_genes)

p_atp   <- plot_panel(long_atp,   transport_genes, ncol = 3)
p_synth <- plot_panel(long_synth, synth_genes,     ncol = 3)
p_catab <- plot_panel(long_catab, catab_genes,     ncol = 3)

print(p_atp)
print(p_synth)
print(p_catab)

save_tiff(p_atp,   "DepMap_BRCA_transport_subtypes")
save_tiff(p_synth, "DepMap_BRCA_synthesis_subtypes")
save_tiff(p_catab, "DepMap_BRCA_catabolism_subtypes")

# ==========================================================
# CAT: CELL LINES USED (ACTUALLY PLOTTED)
# ==========================================================
cat_cell_lines_plotted <- function(df_long, label) {
  
  cat("\n\n==============================\n")
  cat("CELL LINES ACTUALLY PLOTTED\n")
  cat(label, "\n")
  cat("= unique DepMap_ID present after join + non-NA expr\n")
  cat("==============================\n")
  
  plotted <- df_long %>%
    dplyr::distinct(DepMap_ID, cell_line_name, Subtype) %>%
    dplyr::arrange(Subtype, cell_line_name)
  
  plotted %>%
    dplyr::group_by(Subtype) %>%
    dplyr::summarise(
      n_cell_lines = dplyr::n_distinct(DepMap_ID),
      cell_lines   = paste(sort(unique(cell_line_name)), collapse = ", "),
      .groups = "drop"
    ) %>%
    split(.$Subtype) %>%
    purrr::iwalk(function(d, st) {
      cat("\n------------------------------\n")
      cat("Subtype:", st, "\n")
      cat("n =", d$n_cell_lines[1], "\n")
      cat(d$cell_lines[1], "\n")
    })
  
  invisible(plotted)
}

plotted_atp   <- cat_cell_lines_plotted(long_atp,   "Transport panel (ATP13A2/ATP13A3/ATP13A4)")
plotted_synth <- cat_cell_lines_plotted(long_synth, "Synthesis panel (ODC1/SRM/SMS)")
plotted_catab <- cat_cell_lines_plotted(long_catab, "Catabolism panel (SMOX/PAOX/SAT1)")

# ==========================================================
# EXTRA: SELECTED CELL LINES PANELS (WITH SAME p-value SIZE)
# ==========================================================
sel_lines <- c("MCF7","BT474","MDAMB453","MDAMB231","SUM149PT")

sel_cols <- c(
  "MCF7"     = "darkkhaki",
  "BT474"    = "lightpink",
  "MDAMB453" = "orange",
  "MDAMB231" = "lightblue",
  "SUM149PT" = "lightblue"
)

build_long_for_lines <- function(gene_set) {
  
  gene_cols <- setNames(
    vapply(gene_set, find_gene_col, character(1), cn = colnames(expr)),
    gene_set
  )
  
  missing <- names(gene_cols)[is.na(gene_cols)]
  if (length(missing) > 0) {
    warning("Missing genes (skipped): ", paste(missing, collapse = ", "))
    gene_cols <- gene_cols[!is.na(gene_cols)]
  }
  stopifnot(length(gene_cols) > 0)
  
  expr %>%
    dplyr::select(DepMap_ID, dplyr::all_of(unname(gene_cols))) %>%
    dplyr::inner_join(meta_brca, by = "DepMap_ID") %>%
    dplyr::mutate(
      CellLine = toupper(gsub("[^A-Za-z0-9]", "", cell_line_name))
    ) %>%
    dplyr::filter(CellLine %in% sel_lines) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(unname(gene_cols)),
      names_to  = "gene_col",
      values_to = "expr_value"
    ) %>%
    dplyr::mutate(
      gene     = stringr::str_trim(sub(" \\(.*", "", gene_col)),
      gene     = factor(gene, levels = gene_set),
      CellLine = factor(CellLine, levels = sel_lines)
    ) %>%
    dplyr::filter(!is.na(expr_value), !is.na(gene), !is.na(CellLine))
}

compute_stats_grouped <- function(df_long, gene_levels, group_var = "CellLine") {
  
  df_long <- df_long %>% dplyr::mutate(gene = factor(gene, levels = gene_levels))
  
  overall <- purrr::map_dfr(gene_levels, function(g) {
    d <- df_long %>% dplyr::filter(gene == g)
    if (dplyr::n_distinct(d[[group_var]]) < 2) return(tibble::tibble(gene = g, p_anova = NA_real_))
    p <- tryCatch(summary(aov(expr_value ~ .data[[group_var]], data = d))[[1]][["Pr(>F)"]][1],
                  error = function(e) NA_real_)
    tibble::tibble(gene = g, p_anova = p)
  }) %>%
    dplyr::mutate(
      gene = factor(gene, levels = gene_levels),
      gene_label = paste0(
        "<b><span style='font-size:15pt;'>", as.character(gene), "</span></b>",
        "<br><span style='font-size:11pt;'>ANOVA p=", fmt_p(p_anova), "</span>"
      )
    )
  
  pairwise <- purrr::map_dfr(gene_levels, function(g) {
    d <- df_long %>% dplyr::filter(gene == g)
    if (dplyr::n_distinct(d[[group_var]]) < 2) return(tibble::tibble())
    
    tk <- tryCatch(
      TukeyHSD(aov(expr_value ~ .data[[group_var]], data = d), group_var)[[group_var]] %>% as.data.frame(),
      error = function(e) NULL
    )
    if (is.null(tk)) return(tibble::tibble())
    
    tk$comparison <- rownames(tk)
    
    out <- tibble::as_tibble(tk) %>%
      dplyr::mutate(comparison = tk$comparison) %>%
      dplyr::bind_cols(split_tukey_comparison(tk$comparison)) %>%
      dplyr::rename(p_adj = `p adj`) %>%
      dplyr::filter(!is.na(p_adj), p_adj < 0.05)
    
    if (nrow(out) == 0) return(tibble::tibble())
    
    ymax <- max(d$expr_value, na.rm = TRUE)
    ymin <- min(d$expr_value, na.rm = TRUE)
    step <- 0.10 * (ymax - ymin + 1e-6)
    
    out %>%
      dplyr::mutate(
        gene = factor(g, levels = gene_levels),
        y.position = ymax + step * dplyr::row_number(),
        label = paste0("p=", fmt_p(p_adj))
      ) %>%
      dplyr::select(gene, group1, group2, p_adj, label, y.position)
  })
  
  list(overall = overall, pairwise = pairwise, df_long = df_long)
}

plot_panel_lines <- function(df_long, gene_levels, ncol = 3) {
  
  st <- compute_stats_grouped(df_long, gene_levels, group_var = "CellLine")
  df_long <- st$df_long
  
  sum_df <- df_long %>%
    dplyr::group_by(gene, CellLine) %>%
    dplyr::summarise(
      mean = mean(expr_value, na.rm = TRUE),
      sd   = sd(expr_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene     = factor(gene, levels = gene_levels),
      CellLine = factor(CellLine, levels = sel_lines)
    )
  
  lab_vec <- setNames(st$overall$gene_label, as.character(st$overall$gene))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = sum_df,
      ggplot2::aes(x = CellLine, y = mean, fill = CellLine),
      width = bar_width, color = "black", alpha = bar_alpha
    ) +
    ggplot2::geom_errorbar(
      data = sum_df,
      ggplot2::aes(x = CellLine, ymin = mean - sd, ymax = mean + sd),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = df_long,
      ggplot2::aes(x = CellLine, y = expr_value),
      color = "black",
      width = dot_jit_w,
      size  = dot_size,
      alpha = dot_alpha
    ) +
    ggplot2::scale_fill_manual(values = sel_cols, drop = FALSE) +
    ggplot2::facet_wrap(
      ~ gene,
      ncol = ncol,
      scales = "fixed",
      drop = FALSE,
      labeller = ggplot2::labeller(gene = lab_vec)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "**mRNA expression**<br><span style='font-weight:400;'>mean \u00B1 SD; DepMap</span>"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text   = ggtext::element_markdown(),
      axis.text.x  = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_y_size),
      axis.title.y = ggtext::element_markdown(size = y_title_sz),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none"
    )
  
  # --- NEW: bracket p-value text size + bracket line thickness ---
  if (nrow(st$pairwise) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data         = st$pairwise,
        label        = "label",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = tip_len,
        hide.ns      = TRUE,
        size         = pval_bracket_size,
        bracket.size = bracket_size
      )
  }
  
  p
}

# ---- RUN for selected lines ----
long_atp_lines   <- build_long_for_lines(transport_genes)
long_synth_lines <- build_long_for_lines(synth_genes)
long_catab_lines <- build_long_for_lines(catab_genes)

p_atp_lines   <- plot_panel_lines(long_atp_lines,   transport_genes, ncol = 3)
p_synth_lines <- plot_panel_lines(long_synth_lines, synth_genes,     ncol = 3)
p_catab_lines <- plot_panel_lines(long_catab_lines, catab_genes,     ncol = 3)

print(p_atp_lines)
print(p_synth_lines)
print(p_catab_lines)

save_tiff(p_atp_lines,   "DepMap_BRCA_transport_selected_lines")
save_tiff(p_synth_lines, "DepMap_BRCA_synthesis_selected_lines")
save_tiff(p_catab_lines, "DepMap_BRCA_catabolism_selected_lines")



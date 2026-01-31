# ==========================================================
# DepMap (CCLE) — CAV1 only (UPDATED)
# Change requested:
# - Make the facet strip gene name + ANOVA line BIGGER (larger font)
#   (both in selected-lines plot + subtype plot)
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","ggpubr","readr","purrr","tibble","stringr","ggtext")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ---------------- SETTINGS ----------------
setwd("C:/Users/u0171728/Downloads")

expr_file <- "CCLE_expression_full.csv"
meta_file <- "sample_info.csv"

sel_lines <- c("MCF7","BT474","MDAMB453","MDAMB231","SUM149PT")

subtype_levels <- c("LumA","LumB","HER2","B-L")
subtype_cols <- c(
  "LumA" = "darkkhaki",
  "LumB" = "lightpink",
  "HER2" = "orange",
  "B-L"  = "lightblue"
)

gene_of_interest <- "CAV1"

# TIFF export settings
tiff_dpi  <- 600
tiff_w_in <- 5.8
tiff_h_in <- 4.2

save_tiff <- function(plot_obj, filename_no_ext, w = tiff_w_in, h = tiff_h_in) {
  ggplot2::ggsave(
    filename = paste0(filename_no_ext, ".tiff"),
    plot = plot_obj,
    device = "tiff",
    dpi = tiff_dpi,
    width = w,
    height = h,
    units = "in",
    compression = "lzw"
  )
}

# ---------------- STYLE ----------------
base_size    <- 17
axis_x_size  <- 16
axis_y_size  <- 16
y_title_sz   <- 16

dot_jit_w <- 0.18
dot_size  <- 1.8

dot_alpha_sel <- 0.75
dot_alpha_sub <- 0.65

bar_width <- 0.65
bar_alpha <- 0.75

pval_bracket_size <- 4.5
bracket_size      <- 0.4
tip_len           <- 0.01

# ---- NEW: strip label sizes (bigger) ----
STRIP_GENE_PT  <- 20  # was 15
STRIP_ANOVA_PT <- 14  # was 11
STRIP_TEXT_SZ  <- 16  # element_markdown size (overall scaling)

# ---------------- HELPERS ----------------
find_gene_col <- function(gene, cn) {
  gene_u <- toupper(gene)
  
  hit <- grep(paste0("^", gene_u, " \\("), toupper(cn), value = TRUE)
  if (length(hit) > 0) return(cn[match(hit[1], toupper(cn))])
  
  hit2 <- cn[toupper(cn) == gene_u]
  if (length(hit2) > 0) return(hit2[1])
  
  hit3 <- grep(paste0("^", gene_u, "\\b"), toupper(cn), value = TRUE)
  if (length(hit3) > 0) return(cn[match(hit3[1], toupper(cn))])
  
  NA_character_
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

split_tukey_comparison <- function(x) {
  group2 <- sub("^.*-", "", x)
  group1 <- sub("-[^-]+$", "", x)
  tibble::tibble(group1 = group1, group2 = group2)
}

# ---- UPDATED: bigger strip label ----
make_strip_label <- function(gene, p_anova) {
  paste0(
    "<b><span style='font-size:", STRIP_GENE_PT, "pt;'>", gene, "</span></b>",
    "<br><span style='font-size:", STRIP_ANOVA_PT, "pt;'>ANOVA p=", fmt_p(p_anova), "</span>"
  )
}

tukey_sig <- function(d, formula_str, factor_name, value_col) {
  tk <- tryCatch(
    TukeyHSD(aov(as.formula(formula_str), data = d), factor_name)[[factor_name]] %>% as.data.frame(),
    error = function(e) NULL
  )
  if (is.null(tk) || nrow(tk) == 0) return(tibble::tibble())
  
  tk$comparison <- rownames(tk)
  
  out <- tibble::as_tibble(tk) %>%
    dplyr::mutate(comparison = tk$comparison) %>%
    dplyr::bind_cols(split_tukey_comparison(tk$comparison)) %>%
    dplyr::rename(p_adj = `p adj`) %>%
    dplyr::filter(!is.na(p_adj), p_adj < 0.05)
  
  if (nrow(out) == 0) return(tibble::tibble())
  
  ymax <- max(d[[value_col]], na.rm = TRUE)
  ymin <- min(d[[value_col]], na.rm = TRUE)
  step <- 0.10 * (ymax - ymin + 1e-6)
  
  out %>%
    dplyr::mutate(
      y.position = ymax + step * dplyr::row_number(),
      label = paste0("p=", fmt_p(p_adj))
    ) %>%
    dplyr::select(group1, group2, p_adj, label, y.position)
}

# ---------------- LOAD ----------------
expr <- readr::read_csv(expr_file, show_col_types = FALSE)
meta <- readr::read_csv(meta_file, show_col_types = FALSE)
colnames(expr)[1] <- "DepMap_ID"

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
    Subtype = dplyr::case_when(
      toupper(cell_line_name) %in% c("MDA-MB-453","MDAMB453","MDA MB 453","MDA_MB_453") ~ "HER2",
      TRUE ~ Subtype
    ),
    Subtype  = factor(Subtype, levels = subtype_levels),
    CellLine = toupper(gsub("[^A-Za-z0-9]", "", cell_line_name))
  ) %>%
  dplyr::filter(!is.na(Subtype)) %>%
  dplyr::distinct(DepMap_ID, cell_line_name, CellLine, Subtype)

# CAV1 column
cav1_col <- find_gene_col(gene_of_interest, colnames(expr))
if (is.na(cav1_col)) stop("Could not find CAV1 column in CCLE_expression_full.csv")

expr_cav1 <- expr %>%
  dplyr::select(DepMap_ID, dplyr::all_of(cav1_col)) %>%
  dplyr::rename(expr_value = !!cav1_col)

df_all <- expr_cav1 %>%
  dplyr::inner_join(meta_brca, by = "DepMap_ID") %>%
  dplyr::mutate(
    gene = gene_of_interest,
    Subtype = factor(Subtype, levels = subtype_levels)
  ) %>%
  dplyr::filter(!is.na(expr_value))

df_sel <- df_all %>%
  dplyr::filter(CellLine %in% sel_lines) %>%
  dplyr::mutate(CellLine = factor(CellLine, levels = sel_lines))

# ==========================================================
# 1) Selected lines plot (FORCE ANOVA p=NA)
# ==========================================================
sum_sel <- df_sel %>%
  dplyr::group_by(gene, CellLine, Subtype) %>%
  dplyr::summarise(
    mean = mean(expr_value, na.rm = TRUE),
    sd   = sd(expr_value, na.rm = TRUE),
    .groups = "drop"
  )

strip_sel <- setNames(make_strip_label("CAV1", NA_real_), "CAV1")
br_sel <- tukey_sig(df_sel, "expr_value ~ CellLine", "CellLine", "expr_value")

p_sel <- ggplot2::ggplot() +
  ggplot2::geom_col(
    data = sum_sel,
    ggplot2::aes(x = CellLine, y = mean, fill = Subtype),
    width = bar_width, color = "black", alpha = bar_alpha
  ) +
  ggplot2::geom_errorbar(
    data = sum_sel,
    ggplot2::aes(x = CellLine, ymin = mean - sd, ymax = mean + sd),
    width = 0.2
  ) +
  ggplot2::geom_jitter(
    data = df_sel,
    ggplot2::aes(x = CellLine, y = expr_value),
    color = "black",
    width = dot_jit_w,
    size  = dot_size,
    alpha = dot_alpha_sel
  ) +
  ggplot2::scale_fill_manual(values = subtype_cols, drop = FALSE) +
  ggplot2::facet_wrap(
    ~ gene,
    ncol = 1,
    scales = "fixed",
    drop = FALSE,
    labeller = ggplot2::labeller(gene = strip_sel)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "**mRNA expression**<br><span style='font-weight:400;'>mean \u00B1 SD; DepMap</span>"
  ) +
  ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(
    strip.text   = ggtext::element_markdown(size = STRIP_TEXT_SZ, lineheight = 0.95),
    strip.background = ggplot2::element_rect(fill = "grey85", color = "black"),
    axis.text.x  = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
    axis.text.y  = ggplot2::element_text(size = axis_y_size),
    axis.title.y = ggtext::element_markdown(size = y_title_sz),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "none"
  )

if (nrow(br_sel) > 0) {
  p_sel <- p_sel +
    ggpubr::stat_pvalue_manual(
      data         = br_sel,
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

print(p_sel)
save_tiff(p_sel, "DepMap_CAV1_selected_cell_lines_ANOVA_NA_bigStrip", w = 5.8, h = 4.2)

# ==========================================================
# 2) Subtype plot (darker dots)
# ==========================================================
sum_sub <- df_all %>%
  dplyr::group_by(gene, Subtype) %>%
  dplyr::summarise(
    mean = mean(expr_value, na.rm = TRUE),
    sd   = sd(expr_value, na.rm = TRUE),
    .groups = "drop"
  )

p_anova_sub <- tryCatch(
  summary(aov(expr_value ~ Subtype, data = df_all))[[1]][["Pr(>F)"]][1],
  error = function(e) NA_real_
)

strip_sub <- setNames(make_strip_label("CAV1", p_anova_sub), "CAV1")
br_sub <- tukey_sig(df_all, "expr_value ~ Subtype", "Subtype", "expr_value")

p_sub <- ggplot2::ggplot() +
  ggplot2::geom_col(
    data = sum_sub,
    ggplot2::aes(x = Subtype, y = mean, fill = Subtype),
    width = bar_width, color = "black", alpha = bar_alpha
  ) +
  ggplot2::geom_errorbar(
    data = sum_sub,
    ggplot2::aes(x = Subtype, ymin = mean - sd, ymax = mean + sd),
    width = 0.2
  ) +
  ggplot2::geom_jitter(
    data = df_all,
    ggplot2::aes(x = Subtype, y = expr_value),
    color = "black",
    width = dot_jit_w,
    size  = dot_size,
    alpha = dot_alpha_sub
  ) +
  ggplot2::scale_fill_manual(values = subtype_cols, drop = FALSE) +
  ggplot2::facet_wrap(
    ~ gene,
    ncol = 1,
    scales = "fixed",
    drop = FALSE,
    labeller = ggplot2::labeller(gene = strip_sub)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "**mRNA expression**<br><span style='font-weight:400;'>mean \u00B1 SD; DepMap</span>"
  ) +
  ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(
    strip.text   = ggtext::element_markdown(size = STRIP_TEXT_SZ, lineheight = 0.95),
    strip.background = ggplot2::element_rect(fill = "grey85", color = "black"),
    axis.text.x  = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
    axis.text.y  = ggplot2::element_text(size = axis_y_size),
    axis.title.y = ggtext::element_markdown(size = y_title_sz),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "none"
  )

if (nrow(br_sub) > 0) {
  p_sub <- p_sub +
    ggpubr::stat_pvalue_manual(
      data         = br_sub,
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

print(p_sub)
save_tiff(p_sub, "DepMap_CAV1_bySubtype_all_lines_darkerDots_bigStrip", w = 5.8, h = 4.2)

cat("\nDone. Strip gene name + ANOVA text are larger now.\n")

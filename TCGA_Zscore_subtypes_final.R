# ==========================================================
# TCGA BRCA — Polyamine MODULE activity by PAM50 subtype (DepMap-style)
# MATCHES your TCGA gene-per-category script:
# - DepMap colors + sizes
# - Violin + box + jitter (black, lighter alpha=0.25; size 1.8; jitter width 0.18)
# - Facet strip (2 lines, markdown):
#     line 1: MODULE (15 pt, bold)
#     line 2: Kruskal–Wallis p=... (11 pt, regular)
# - Pairwise: Wilcoxon vs Normal only, BH-adjusted; only p<0.05 shown
# - Brackets ON (and preserved in TIFF using ragg::agg_tiff)
# - Save TIFF (600 dpi, LZW) to Downloads
# - Extra: cat N samples per subtype (overall; unique PATIENT_ID)
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","ggpubr","ggtext","stringr","grid","ragg")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ------------------ PATHS ------------------
DATA_DIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
OUT_DIR  <- "C:/Users/u0171728/Downloads"
setwd(DATA_DIR)

EXPR_FILE <- "data_mrna_seq_v2_rsem.txt"
CLIN_FILE <- "data_clinical_patient.txt"

# ------------------ MODULES ------------------
GENES <- c("ODC1","SMS","SRM",
           "ATP13A2","ATP13A3","ATP13A4",
           "SMOX","PAOX","SAT1")

synth_genes <- c("ODC1","SMS","SRM")
catab_genes <- c("SMOX","PAOX","SAT1")
trans_genes <- c("ATP13A2","ATP13A3","ATP13A4")

# ------------------ PAM50 ------------------
pam50_codes  <- c("BRCA_Normal","BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("Normal","LumA","LumB","HER2","Basal-like")

# DepMap-like palette (same as your CCLE/DepMap scripts)
pam50_cols <- c(
  "Normal"     = "grey80",
  "LumA"       = "darkkhaki",
  "LumB"       = "lightpink",
  "HER2"       = "orange",
  "Basal-like" = "lightblue"
)

# ------------------ STYLE (DepMap) ------------------
base_size    <- 20
strip_genept <- 20
strip_ppt    <- 20
axis_x_size  <- 20
axis_y_size  <- 20
y_title_sz   <- 20

dot_jit_w <- 0.18
dot_size  <- 1.8
dot_alpha <- 0.25

pval_text_size <- 6
bracket_size   <- 0.4
tip_len        <- 0.01

# TIFF export (brackets preserved)
tiff_dpi  <- 600
tiff_w_in <- 10
tiff_h_in <- 7

save_tiff <- function(plot_obj, filename_no_ext) {
  ggplot2::ggsave(
    filename    = file.path(OUT_DIR, paste0(filename_no_ext, ".tiff")),
    plot        = plot_obj,
    device      = ragg::agg_tiff,
    width       = tiff_w_in,
    height      = tiff_h_in,
    units       = "in",
    res         = tiff_dpi,
    compression = "lzw"
  )
}

# ------------------ HELPERS ------------------
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

stop_if_empty <- function(df, msg) {
  if (nrow(df) == 0) stop(msg, call. = FALSE)
}

# ------------------ LOAD CLINICAL ------------------
clin <- read.delim(CLIN_FILE, check.names = FALSE, comment.char = "#") %>%
  dplyr::select(PATIENT_ID, SUBTYPE, SEX) %>%
  dplyr::filter(toupper(as.character(SEX)) == "FEMALE") %>%
  dplyr::mutate(
    SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels),
    SUBTYPE = factor(as.character(SUBTYPE), levels = pam50_labels)
  ) %>%
  dplyr::filter(!is.na(SUBTYPE))

stop_if_empty(clin, "Clinical table empty after filtering to FEMALE + valid PAM50 SUBTYPE.")

# ------------------ LOAD EXPRESSION ------------------
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

present_genes <- intersect(GENES, unique(expr$Hugo_Symbol))
if (length(present_genes) == 0) stop("None of the module genes found in expression file.", call. = FALSE)

expr_long <- expr %>%
  dplyr::filter(Hugo_Symbol %in% present_genes) %>%
  tidyr::pivot_longer(cols = -Hugo_Symbol, names_to = "SAMPLE_ID", values_to = "expr_raw") %>%
  dplyr::mutate(
    expr       = log2(as.numeric(expr_raw) + 1),
    PATIENT_ID = substr(SAMPLE_ID, 1, 12),
    gene       = Hugo_Symbol
  ) %>%
  dplyr::select(PATIENT_ID, gene, expr) %>%
  dplyr::filter(!is.na(expr))

stop_if_empty(expr_long, "Expression long table empty after filtering genes/log2 transform.")

expr_with_subtype <- expr_long %>%
  dplyr::inner_join(clin %>% dplyr::select(PATIENT_ID, SUBTYPE), by = "PATIENT_ID") %>%
  dplyr::filter(!is.na(SUBTYPE))

stop_if_empty(expr_with_subtype, "No overlap between expression and clinical PATIENT_IDs.")

# ------------------ N SAMPLES PER SUBTYPE (overall) ------------------
cat("\n==============================\n")
cat("N SAMPLES PER PAM50 SUBTYPE (overall; unique PATIENT_ID)\n")
cat("==============================\n")

expr_with_subtype %>%
  distinct(PATIENT_ID, SUBTYPE) %>%
  count(SUBTYPE, name = "n_samples") %>%
  arrange(SUBTYPE) %>%
  print(n = Inf)

# ------------------ NORMAL-ANCHORED Z (per gene) ------------------
norm_params <- expr_with_subtype %>%
  dplyr::filter(SUBTYPE == "Normal") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    mu_norm = mean(expr, na.rm = TRUE),
    sd_norm = sd(expr, na.rm = TRUE),
    .groups = "drop"
  )

expr_norm <- expr_with_subtype %>%
  dplyr::inner_join(norm_params, by = "gene") %>%
  dplyr::mutate(
    z_norm = dplyr::if_else(!is.na(sd_norm) & sd_norm > 0,
                            (expr - mu_norm) / sd_norm,
                            NA_real_)
  ) %>%
  dplyr::select(PATIENT_ID, gene, SUBTYPE, z_norm)

# ------------------ MODULE SCORES (mean z_norm) ------------------
get_module_score <- function(df, genes, score_name) {
  df %>%
    dplyr::filter(gene %in% genes) %>%
    dplyr::group_by(PATIENT_ID, SUBTYPE) %>%
    dplyr::summarise(
      !!score_name := mean(z_norm, na.rm = TRUE),
      .groups = "drop"
    )
}

scores <- get_module_score(expr_norm, synth_genes, "score") %>%
  dplyr::mutate(module = "Synthesis") %>%
  dplyr::bind_rows(
    get_module_score(expr_norm, catab_genes, "score") %>% dplyr::mutate(module = "Catabolism"),
    get_module_score(expr_norm, trans_genes, "score") %>% dplyr::mutate(module = "Transport")
  ) %>%
  dplyr::filter(!is.na(score), !is.na(SUBTYPE)) %>%
  dplyr::mutate(
    module = factor(module, levels = c("Synthesis","Catabolism","Transport")),
    SUBTYPE = factor(as.character(SUBTYPE), levels = pam50_labels)
  )

stop_if_empty(scores, "Module score table is empty (all scores NA?).")

# ------------------ KRUSKAL–WALLIS per module (facet strip) ------------------
kw_df <- scores %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    p_kw = tryCatch(
      kruskal.test(score ~ SUBTYPE, data = dplyr::pick(dplyr::everything()))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    module_label = paste0(
      "<b><span style='font-size:", strip_genept, "pt;'>", as.character(module), "</span></b>",
      "<br><span style='font-size:", strip_ppt, "pt;'>Kruskal–Wallis p=", fmt_p(p_kw), "</span>"
    )
  )

scores2 <- scores %>%
  dplyr::left_join(kw_df, by = "module")

lab_vec <- setNames(kw_df$module_label, as.character(kw_df$module))

# ------------------ WILCOXON vs NORMAL (BH), only significant ------------------
pw <- ggpubr::compare_means(
  score ~ SUBTYPE,
  data            = scores2,
  method          = "wilcox.test",
  p.adjust.method = "BH",
  group.by        = "module",
  ref.group       = "Normal"
)

pw_sig <- pw %>%
  dplyr::filter(group1 == "Normal", !is.na(p.adj), p.adj < 0.05) %>%
  dplyr::mutate(
    module = factor(as.character(module), levels = levels(scores2$module)),
    label  = paste0("p=", fmt_p(p.adj))
  )

if (nrow(pw_sig) > 0) {
  yr <- scores2 %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      y_max = max(score, na.rm = TRUE),
      y_min = min(score, na.rm = TRUE),
      rng   = y_max - y_min,
      .groups = "drop"
    )
  
  pw_sig <- pw_sig %>%
    dplyr::left_join(yr, by = "module") %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(p.adj, .by_group = TRUE) %>%
    dplyr::mutate(
      step       = ifelse(rng <= 0, 0.5, 0.18 * rng),
      y.position = y_max + step * dplyr::row_number()
    ) %>%
    dplyr::ungroup()
}

# ------------------ PLOT (same style as gene script) ------------------
p <- ggplot2::ggplot(scores2, ggplot2::aes(x = SUBTYPE, y = score, fill = SUBTYPE)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.2) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, linewidth = 0.2) +
  ggplot2::geom_jitter(width = dot_jit_w, alpha = dot_alpha, size = dot_size, color = "black") +
  ggplot2::scale_fill_manual(values = pam50_cols, drop = FALSE) +
  ggplot2::facet_wrap(
    ~ module,
    nrow = 1,
    scales = "free_y",
    labeller = ggplot2::labeller(module = lab_vec)
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.25))) +
  ggplot2::labs(
    x = NULL,
    y = "**Mean z-score vs Normal**<br><span style='font-weight:400;'>module activity; mean z (Normal-anchored)</span>"
  ) +
  ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
    strip.text       = ggtext::element_markdown(),
    axis.text.x      = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
    axis.text.y      = ggplot2::element_text(size = axis_y_size),
    axis.title.y     = ggtext::element_markdown(size = y_title_sz),
    legend.position  = "none",
    panel.spacing    = grid::unit(0.7, "lines")
  )

if (nrow(pw_sig) > 0) {
  p <- p +
    ggpubr::stat_pvalue_manual(
      data         = pw_sig,
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y.position",
      label        = "label",
      tip.length   = tip_len,
      bracket.size = bracket_size,
      size         = pval_text_size
    )
}

print(p)
save_tiff(p, "TCGA_BRCA_PAM50_polyamine_module_activity_DepMapStyle")

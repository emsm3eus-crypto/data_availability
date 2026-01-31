# ==========================================================
# TCGA BRCA — CAV1 by PAM50 SUBTYPE (DepMap-style)
# - Kruskal–Wallis per gene (facet strip)
# - Wilcoxon vs Normal only, BH adjusted; only p<0.05 shown
# - CHANGE: put the p-value brackets/labels UNDER the violins (bottom),
#           because there is no space above.
# - SAVE: TIFF 600 dpi (LZW) to Downloads using ragg (keeps bracket lines)
# - EXTRA: cat N samples per subtype (overall; unique PATIENT_ID)
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

# ------------------ GENE ------------------
GENE <- "CAV1"

# ------------------ PAM50 ------------------
pam50_codes  <- c("BRCA_Normal","BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("Normal","LumA","LumB","Her2","Basal-like")

pam50_cols <- c(
  "Normal"     = "grey80",
  "LumA"       = "darkkhaki",
  "LumB"       = "lightpink",
  "Her2"       = "orange",
  "Basal-like" = "lightblue"
)

# ------------------ STYLE (DepMap) ------------------
base_size    <- 17
strip_genept <- 20
strip_ppt    <- 15
axis_x_size  <- 17
axis_y_size  <- 17
y_title_sz   <- 16

dot_jit_w <- 0.18
dot_size  <- 1.8
dot_alpha <- 0.25

pval_text_size <- 5.5
bracket_size   <- 0.4
tip_len        <- 0.01

# ------------------ TIFF EXPORT (robust brackets) ------------------
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

# ------------------ LOAD DATA ------------------
clin <- read.delim("data_clinical_patient.txt", check.names = FALSE, comment.char = "#")
expr <- read.delim("data_mrna_seq_v2_rsem.txt", check.names = FALSE)
stopifnot("Hugo_Symbol" %in% colnames(expr))

if (!(GENE %in% unique(expr$Hugo_Symbol))) {
  stop(paste0("Gene ", GENE, " not found in data_mrna_seq_v2_rsem.txt (Hugo_Symbol column)."), call. = FALSE)
}

clin2 <- clin %>%
  dplyr::select(PATIENT_ID, SUBTYPE, SEX) %>%
  dplyr::filter(toupper(as.character(SEX)) == "FEMALE") %>%
  dplyr::mutate(SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels)) %>%
  dplyr::filter(!is.na(SUBTYPE))

stop_if_empty(clin2, "Clinical table empty after filtering to FEMALE + valid PAM50 SUBTYPE.")

expr_long <- expr %>%
  dplyr::filter(Hugo_Symbol == GENE) %>%
  tidyr::pivot_longer(cols = -Hugo_Symbol, names_to = "SAMPLE_ID", values_to = "expr_raw") %>%
  dplyr::mutate(
    expr       = log2(as.numeric(expr_raw) + 1),
    PATIENT_ID = substr(SAMPLE_ID, 1, 12),
    gene       = Hugo_Symbol
  ) %>%
  dplyr::select(PATIENT_ID, gene, expr) %>%
  dplyr::filter(!is.na(expr))

stop_if_empty(expr_long, "Expression long table empty after gene filter/log2 transform.")

dat <- dplyr::inner_join(expr_long, clin2, by = "PATIENT_ID") %>%
  dplyr::filter(!is.na(expr), !is.na(SUBTYPE)) %>%
  dplyr::mutate(
    SUBTYPE = factor(as.character(SUBTYPE), levels = pam50_labels),
    gene    = factor(as.character(gene), levels = GENE)
  )

stop_if_empty(dat, "Joined expression+clinical table empty. Likely PATIENT_ID mismatch or missing SUBTYPE.")

# ==========================================================
# EXTRA: CAT N SAMPLES PER SUBTYPE (overall; unique PATIENT_ID)
# ==========================================================
cat("\n==============================\n")
cat("N SAMPLES PER PAM50 SUBTYPE (overall; unique PATIENT_ID)\n")
cat("==============================\n")

dat %>%
  distinct(PATIENT_ID, SUBTYPE) %>%
  count(SUBTYPE, name = "n_samples") %>%
  arrange(SUBTYPE) %>%
  print(n = Inf)

# ------------------ PLOT BUILDER ------------------
build_gene_plot_bottom_pvals <- function() {
  
  d_gene <- dat
  
  # y-limits + extra space below for brackets
  y_high <- max(d_gene$expr, na.rm = TRUE)
  y_low  <- min(d_gene$expr, na.rm = TRUE)
  if (!is.finite(y_high) || !is.finite(y_low)) stop("y limits not finite.", call. = FALSE)
  y_rng <- y_high - y_low
  if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1
  
  # Put brackets in the lower margin:
  #  - add ~30% extra below minimum
  y_floor <- y_low - 0.35 * y_rng
  
  # Kruskal–Wallis p -> facet strip
  kw_df <- d_gene %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      p_kw = tryCatch(
        kruskal.test(expr ~ SUBTYPE, data = dplyr::pick(dplyr::everything()))$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      facet_label = paste0(
        "<b><span style='font-size:", strip_genept, "pt;'>", as.character(gene), "</span></b>",
        "<br><span style='font-size:", strip_ppt, "pt;'>Kruskal–Wallis p=", fmt_p(p_kw), "</span>"
      )
    )
  
  d_gene2 <- d_gene %>% dplyr::left_join(kw_df, by = "gene")
  lab_vec <- setNames(kw_df$facet_label, as.character(kw_df$gene))
  
  # Wilcoxon vs Normal only (BH), only significant
  pw <- ggpubr::compare_means(
    expr ~ SUBTYPE,
    data            = d_gene2,
    method          = "wilcox.test",
    p.adjust.method = "BH",
    group.by        = "gene",
    ref.group       = "Normal"
  )
  
  # --- place p-values UNDER the data ---
  pw_sig <- pw %>%
    dplyr::filter(!is.na(p.adj), p.adj < 0.05) %>%
    dplyr::mutate(gene = factor(as.character(gene), levels = GENE)) %>%
    dplyr::group_by(gene) %>%
    dplyr::arrange(p.adj) %>%
    dplyr::mutate(
      # stack downward: first comparison closest to plot, next lower, etc.
      y.position = (y_low - 0.08*y_rng) - (0.08*y_rng)*(dplyr::row_number() - 1),
      label      = paste0("p=", fmt_p(p.adj))
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(d_gene2, ggplot2::aes(x = SUBTYPE, y = expr, fill = SUBTYPE)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, linewidth = 0.2) +
    ggplot2::geom_jitter(width = dot_jit_w, alpha = dot_alpha, size = dot_size, color = "black") +
    ggplot2::scale_fill_manual(values = pam50_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(y_floor, y_high),
      expand = ggplot2::expansion(mult = c(0.02, 0.10))
    ) +
    ggplot2::facet_wrap(
      ~ gene,
      ncol     = 1,
      scales   = "fixed",
      labeller = ggplot2::labeller(gene = lab_vec)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "**log2(RSEM + 1)**<br><span style='font-weight:400;'>TCGA-BRCA; PAM50</span>"
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
        pw_sig,
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        label        = "label",
        tip.length   = tip_len,
        bracket.size = bracket_size,
        size         = pval_text_size
      )
  }
  
  p
}

# ------------------ RUN + SAVE ------------------
p_cav1 <- build_gene_plot_bottom_pvals()
print(p_cav1)
save_tiff(p_cav1, "TCGA_BRCA_PAM50_CAV1_DepMapStyle_bottomPvals")


# ==========================================================
# qPCR (ATP13A2–4 mRNA) + WB (ODC1/SRM/SMS/SMOX protein)
# MATCH DepMap STYLE + BIGGER SIZES (like your previous scripts)
#
# - Two plots: mRNA panel + protein panel
# - Bars = mean ± SD (alpha 0.75 + black outline)
# - Dots = replicates (black; size 1.8; alpha 0.75; jitter width 0.18)
# - Facet strip (2 lines, via markdown):
#     line 1: GENE label (BIG, bold)
#     line 2: ANOVA p=... (same line)
# - Stats: Dunnett vs MCF10A only (emmeans; adjust="mvt")
# - Only significant (p<0.05) brackets shown
# - Export: TIFF 600 dpi LZW in Downloads
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","readxl","stringr","tibble",
            "ggpubr","ggtext","emmeans","purrr","stats")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

setwd("C:/Users/u0171728/Downloads")

# -----------------------------
# INPUT
# -----------------------------
df <- readxl::read_excel("PAgenes_wbqPCR_rstudio_GOOD.xlsx")

names(df) <- names(df) %>%
  stringr::str_replace_all("\\s+$", "") %>%
  stringr::str_replace_all("-", "") %>%
  stringr::str_trim()

df <- df %>%
  dplyr::rename(gene = ...1) %>%
  dplyr::mutate(gene = stringr::str_trim(as.character(gene)))

# -----------------------------
# SETTINGS
# -----------------------------
cell_levels <- c("MCF10A","MCF7","BT474","MDAMB453","MDAMB231","SUM149PT")
control_lvl <- "MCF10A"

mrna_genes    <- c("ATP13A2","ATP13A3","ATP13A4")
protein_genes <- c("ODC1","SRM","SMS","SMOX")  # ORDER FIXED

fill_cols <- c(
  "MCF10A"   = "grey80",
  "MCF7"     = "darkkhaki",
  "BT474"    = "lightpink",
  "MDAMB453" = "orange",
  "MDAMB231" = "lightblue",
  "SUM149PT" = "lightblue"
)

# TIFF export settings
tiff_dpi  <- 600
tiff_w_in <- 10
tiff_h_in <- 7

save_tiff <- function(plot_obj, filename_no_ext) {
  ggplot2::ggsave(
    filename    = file.path(getwd(), paste0(filename_no_ext, ".tiff")),
    plot        = plot_obj,
    device      = "tiff",
    dpi         = tiff_dpi,
    width       = tiff_w_in,
    height      = tiff_h_in,
    units       = "in",
    compression = "lzw"
  )
}

# -----------------------------
# STYLE = MATCH "ERVOOR" (BIGGER)
# -----------------------------
base_size     <- 17
axis_x_size   <- 16
axis_y_size   <- 16
y_title_sz    <- 16

strip_genept  <- 20   # big gene title
strip_anovapt <- 15   # ANOVA p=...

dot_size  <- 1.8
dot_alpha <- 0.75
dot_jit_w <- 0.18

bar_width <- 0.65
bar_alpha <- 0.75

pval_bracket_size <- 4.5
bracket_size      <- 0.4
tip_len           <- 0.01

# -----------------------------
# HELPERS
# -----------------------------
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

make_long <- function(df, genes, suffix_label) {
  df %>%
    dplyr::filter(gene %in% genes) %>%
    tidyr::pivot_longer(-gene, names_to = "cell_line", values_to = "expr") %>%
    dplyr::mutate(
      cell_line  = factor(stringr::str_trim(cell_line), levels = cell_levels),
      gene       = factor(gene, levels = genes),
      gene_label = factor(
        paste0(as.character(gene), " ", suffix_label),
        levels = paste0(genes, " ", suffix_label)
      )
    ) %>%
    dplyr::filter(!is.na(expr), !is.na(cell_line), !is.na(gene_label))
}

mean_sd <- function(long_df) {
  long_df %>%
    dplyr::group_by(gene_label, cell_line) %>%
    dplyr::summarise(
      mean_expr = mean(expr, na.rm = TRUE),
      sd_expr   = sd(expr,  na.rm = TRUE),
      .groups   = "drop"
    )
}

# ANOVA strip (ANOVA p=... on SAME line)
anova_strip <- function(long_df) {
  long_df %>%
    dplyr::group_by(gene_label) %>%
    dplyr::summarise(
      p_anova = {
        d <- dplyr::pick(dplyr::everything())
        tryCatch(summary(aov(expr ~ cell_line, data = d))[[1]][["Pr(>F)"]][1],
                 error = function(e) NA_real_)
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      strip = paste0(
        "<b><span style='font-size:", strip_genept, "pt;'>", as.character(gene_label), "</span></b>",
        "<br><span style='font-size:", strip_anovapt, "pt;'>ANOVA p=", fmt_p(p_anova), "</span>"
      )
    )
}

# Dunnett vs control per gene_label
dunnett_vs_control <- function(long_df, control = "MCF10A") {
  
  d0 <- long_df %>%
    dplyr::mutate(cell_line = stats::relevel(as.factor(cell_line), ref = control))
  
  gene_levels <- levels(d0$gene_label)
  
  res_list <- purrr::map(gene_levels, function(g) {
    
    d <- d0 %>% dplyr::filter(gene_label == g)
    if (dplyr::n_distinct(d$cell_line) < 2) return(tibble::tibble())
    
    fit <- aov(expr ~ cell_line, data = d)
    em  <- emmeans::emmeans(fit, ~ cell_line)
    
    ct <- emmeans::contrast(
      em,
      method = "trt.vs.ctrl",
      ref    = which(levels(d$cell_line) == control),
      adjust = "mvt"
    )
    
    out <- as.data.frame(ct) %>%
      tibble::as_tibble() %>%
      dplyr::transmute(
        gene_label = g,
        contrast   = contrast,
        p_adj      = p.value
      ) %>%
      dplyr::mutate(
        group2 = stringr::str_trim(sub(" - .*", "", contrast)),   # treatment
        group1 = stringr::str_trim(sub("^.* - ", "", contrast))   # control
      ) %>%
      dplyr::filter(group1 == control) %>%
      dplyr::filter(!is.na(p_adj), p_adj < 0.05) %>%
      dplyr::mutate(label = paste0("p=", fmt_p(p_adj)))
    
    if (nrow(out) == 0) return(tibble::tibble())
    
    yr <- d %>%
      dplyr::summarise(
        ymax = max(expr, na.rm = TRUE),
        ymin = min(expr, na.rm = TRUE),
        rng  = ymax - ymin
      )
    
    out %>%
      dplyr::arrange(p_adj) %>%
      dplyr::mutate(
        step       = ifelse(yr$rng <= 0, 0.5, 0.15 * yr$rng),
        y.position = yr$ymax + step * dplyr::row_number()
      ) %>%
      dplyr::select(gene_label, group1, group2, y.position, label)
  })
  
  dplyr::bind_rows(res_list) %>%
    dplyr::mutate(gene_label = factor(gene_label, levels = gene_levels))
}

plot_panel <- function(long_df, mean_df, strip_df, stats_df, ncol = 3, ylab_md) {
  
  lab_vec <- setNames(strip_df$strip, as.character(strip_df$gene_label))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = mean_df,
      ggplot2::aes(x = cell_line, y = mean_expr, fill = cell_line),
      width = bar_width, color = "black", alpha = bar_alpha
    ) +
    ggplot2::geom_errorbar(
      data = mean_df,
      ggplot2::aes(x = cell_line, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = long_df,
      ggplot2::aes(x = cell_line, y = expr),
      color = "black",
      width = dot_jit_w,
      size  = dot_size,
      alpha = dot_alpha
    ) +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::facet_wrap(
      ~ gene_label,
      ncol = ncol,
      scales = "free_y",
      drop = FALSE,
      labeller = ggplot2::labeller(gene_label = lab_vec)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
    ggplot2::labs(x = NULL, y = ylab_md) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text   = ggtext::element_markdown(),
      axis.text.x  = ggplot2::element_text(size = axis_x_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_y_size),
      axis.title.y = ggtext::element_markdown(size = y_title_sz),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none"
    )
  
  if (!is.null(stats_df) && nrow(stats_df) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data         = stats_df,
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        label        = "label",
        tip.length   = tip_len,
        hide.ns      = TRUE,
        size         = pval_bracket_size,
        bracket.size = bracket_size
      )
  }
  
  p
}

# ==========================================================
# PANEL 1: ATP13A2–4 mRNA
# ==========================================================
long_mrna  <- make_long(df, mrna_genes, "mRNA")
mean_mrna  <- mean_sd(long_mrna)
strip_mrna <- anova_strip(long_mrna)
stats_mrna <- dunnett_vs_control(long_mrna, control = control_lvl)

p_mrna <- plot_panel(
  long_df  = long_mrna,
  mean_df  = mean_mrna,
  strip_df = strip_mrna,
  stats_df = stats_mrna,
  ncol     = 3,
  ylab_md  = "**mRNA expression**<br><span style='font-weight:400;'>rel. to MCF10A, GAPDH-normalized; mean \u00B1 SD</span>"
)
print(p_mrna)
save_tiff(p_mrna, "qPCR_ATP13A2_A3_A4_mRNA_Dunnett_meanSD_DepMapStyle")

# ==========================================================
# PANEL 2: ODC1/SRM/SMS/SMOX protein
# ==========================================================
long_prot  <- make_long(df, protein_genes, "protein")
mean_prot  <- mean_sd(long_prot)
strip_prot <- anova_strip(long_prot)
stats_prot <- dunnett_vs_control(long_prot, control = control_lvl)

p_prot <- plot_panel(
  long_df  = long_prot,
  mean_df  = mean_prot,
  strip_df = strip_prot,
  stats_df = stats_prot,
  ncol     = 4,
  ylab_md  = "**Protein level**<br><span style='font-weight:400;'>rel. to MCF10A, GAPDH-normalized; mean \u00B1 SD</span>"
)
print(p_prot)
save_tiff(p_prot, "WB_ODC1_SRM_SMS_SMOX_protein_Dunnett_meanSD_DepMapStyle")

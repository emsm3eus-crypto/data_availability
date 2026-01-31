# ==========================================================
# BODIPY–polyamine uptake (PUT / SPD / SPM)
# facs_final.xlsx → 3-panel plot (DepMap-style)
#
# MATCHED TYPOGRAPHY to your previous DepMap scripts:
# - Larger overall fonts (base/axes/title)
# - Facet strip 2 lines:
#     line 1: polyamine title (BIGGER, bold)
#     line 2: ANOVA p=... (same line)
#
# Stats:
# - One-way ANOVA per polyamine (shown in strip)
# - Dunnett posthoc vs MCF10A only (emmeans; adjust="mvt")
# - Only significant (p<0.05) shown as brackets
#
# Export:
# - TIFF 600 dpi, LZW, saved to Downloads
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","readxl","stringr","tibble",
            "ggpubr","ggtext","emmeans","purrr","stats")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

setwd("C:/Users/u0171728/Downloads")

# -----------------------------
# SETTINGS
# -----------------------------
cell_order  <- c("MCF10A","MCF7","BT474","MDAMB453","MDAMB231","SUM149PT")
control_lvl <- "MCF10A"

poly_labels <- c(
  "PUT" = "PUT-BODIPY uptake",
  "SPD" = "SPD-BODIPY uptake",
  "SPM" = "SPM-BODIPY uptake"
)

bar_colors <- c(
  "MCF10A"   = "grey80",
  "MCF7"     = "darkkhaki",
  "BT474"    = "lightpink",
  "MDAMB453" = "orange",
  "MDAMB231" = "lightblue",
  "SUM149PT" = "lightblue"
)

# -----------------------------
# TYPOGRAPHY (MATCH "ALS ERVOOR")
# -----------------------------
base_size    <- 17
axis_x_size  <- 16
axis_y_size  <- 16
y_title_sz   <- 16

strip_genept <- 20   # big title in strip
strip_anovapt <- 15  # ANOVA p=... line

# Plot geometry
bar_width <- 0.65
bar_alpha <- 0.75

dot_size  <- 1.8
dot_alpha <- 0.75
dot_jit_w <- 0.18

# Bracket p-values
pval_bracket_size <- 4.5
bracket_size      <- 0.4
tip_len           <- 0.01

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

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

# -----------------------------
# LOAD + CLEAN
# -----------------------------
df <- readxl::read_excel("facs_final.xlsx", col_names = TRUE)

names(df) <- names(df) %>%
  stringr::str_replace_all("[\r\n]", "") %>%
  stringr::str_replace_all("-", "") %>%
  stringr::str_trim()

# First column contains PUT/SPD/SPM; in your file it was named "PUT"
df <- df %>%
  dplyr::rename(polyamine = PUT) %>%
  dplyr::mutate(polyamine = toupper(stringr::str_trim(as.character(polyamine))))

# -----------------------------
# LONG FORMAT
# -----------------------------
long <- df %>%
  tidyr::pivot_longer(
    cols      = -polyamine,
    names_to  = "cell_line",
    values_to = "uptake"
  ) %>%
  dplyr::mutate(
    cell_line  = factor(stringr::str_trim(cell_line), levels = cell_order),
    poly_label = unname(poly_labels[polyamine]),
    poly_label = factor(poly_label, levels = unname(poly_labels))
  ) %>%
  dplyr::filter(!is.na(uptake), !is.na(cell_line), !is.na(poly_label))

# -----------------------------
# MEAN ± SD
# -----------------------------
mean_df <- long %>%
  dplyr::group_by(poly_label, cell_line) %>%
  dplyr::summarise(
    mean_uptake = mean(uptake, na.rm = TRUE),
    sd_uptake   = sd(uptake,  na.rm = TRUE),
    .groups     = "drop"
  )

# -----------------------------
# ANOVA per polyamine (facet strip)
# - UPDATED: ANOVA p=... on the same line
# - UPDATED: strip font sizes match earlier scripts
# -----------------------------
anova_df <- long %>%
  dplyr::group_by(poly_label) %>%
  dplyr::summarise(
    p_anova = {
      d <- dplyr::pick(dplyr::everything())
      tryCatch(summary(aov(uptake ~ cell_line, data = d))[[1]][["Pr(>F)"]][1],
               error = function(e) NA_real_)
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    facet_label = paste0(
      "<b><span style='font-size:", strip_genept, "pt;'>", as.character(poly_label), "</span></b>",
      "<br><span style='font-size:", strip_anovapt, "pt;'>ANOVA p=", fmt_p(p_anova), "</span>"
    )
  )

lab_vec <- setNames(anova_df$facet_label, as.character(anova_df$poly_label))

# -----------------------------
# DUNNETT vs MCF10A (per poly_label)
# -----------------------------
poly_levels <- levels(long$poly_label)

dunnett_sig <- purrr::map_dfr(poly_levels, function(pl) {
  
  d <- long %>% dplyr::filter(poly_label == pl)
  if (dplyr::n_distinct(d$cell_line) < 2) return(tibble::tibble())
  
  d <- d %>% dplyr::mutate(cell_line = stats::relevel(as.factor(cell_line), ref = control_lvl))
  
  fit <- aov(uptake ~ cell_line, data = d)
  em  <- emmeans::emmeans(fit, ~ cell_line)
  
  ct <- emmeans::contrast(
    em,
    method = "trt.vs.ctrl",
    ref    = which(levels(d$cell_line) == control_lvl),
    adjust = "mvt"
  )
  
  out <- as.data.frame(ct) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      poly_label = pl,
      contrast   = contrast,
      p_adj      = p.value
    ) %>%
    dplyr::mutate(
      group2 = stringr::str_trim(sub(" - .*", "", contrast)),  # treatment
      group1 = stringr::str_trim(sub("^.* - ", "", contrast))  # control
    ) %>%
    dplyr::filter(group1 == control_lvl) %>%
    dplyr::filter(!is.na(p_adj), p_adj < 0.05) %>%
    dplyr::mutate(label = paste0("p=", fmt_p(p_adj)))
  
  if (nrow(out) == 0) return(tibble::tibble())
  
  yr <- d %>%
    dplyr::summarise(
      y_max = max(uptake, na.rm = TRUE),
      y_min = min(uptake, na.rm = TRUE),
      rng   = y_max - y_min
    )
  
  out %>%
    dplyr::arrange(p_adj) %>%
    dplyr::mutate(
      step       = ifelse(yr$rng <= 0, 0.5, 0.15 * yr$rng),
      y.position = yr$y_max + step * dplyr::row_number()
    ) %>%
    dplyr::select(poly_label, group1, group2, y.position, label)
})

# -----------------------------
# PLOT (DepMap-like)
# -----------------------------
p_bodipy <- ggplot2::ggplot() +
  ggplot2::geom_col(
    data  = mean_df,
    ggplot2::aes(x = cell_line, y = mean_uptake, fill = cell_line),
    width = bar_width,
    color = "black",
    alpha = bar_alpha
  ) +
  ggplot2::geom_errorbar(
    data = mean_df,
    ggplot2::aes(x = cell_line, ymin = mean_uptake - sd_uptake, ymax = mean_uptake + sd_uptake),
    width = 0.2,
    linewidth = 0.6
  ) +
  ggplot2::geom_jitter(
    data = long,
    ggplot2::aes(x = cell_line, y = uptake),
    color = "black",
    width = dot_jit_w,
    size  = dot_size,
    alpha = dot_alpha
  ) +
  ggplot2::scale_fill_manual(values = bar_colors, drop = FALSE) +
  ggplot2::facet_wrap(
    ~ poly_label,
    ncol = 3,
    scales = "free_y",
    drop = FALSE,
    labeller = ggplot2::labeller(poly_label = lab_vec)
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
  ggplot2::labs(
    x = NULL,
    y = "**BODIPY-polyamine**<br><span style='font-weight:400;'>rel. to MCF10A; mean \u00B1 SD</span>"
  ) +
  ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(
    strip.text   = ggtext::element_markdown(),
    axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = axis_x_size),
    axis.text.y  = ggplot2::element_text(size = axis_y_size),
    axis.title.y = ggtext::element_markdown(size = y_title_sz),
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank()
  )

if (nrow(dunnett_sig) > 0) {
  p_bodipy <- p_bodipy +
    ggpubr::stat_pvalue_manual(
      data         = dunnett_sig,
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

print(p_bodipy)

# -----------------------------
# SAVE TIFF
# -----------------------------
save_tiff(p_bodipy, "BODIPY_polyamine_uptake_Dunnett_meanSD_DepMapStyle")

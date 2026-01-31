# ==========================================================
# TCGA BRCA — Spearman correlation: CAV1 vs ATP13A3
# per PAM50 subtype (NO Normal; facets naast elkaar)
#
# - Removed Normal subtype
# - CAV1 and ATP13A3 in **bold**
# - Facet subtitle shows: rho + p
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","ggtext","stringr","ragg")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ------------------ PATHS ------------------
DATA_DIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
OUT_DIR  <- "C:/Users/u0171728/Downloads"
setwd(DATA_DIR)

EXPR_FILE <- "data_mrna_seq_v2_rsem.txt"
CLIN_FILE <- "data_clinical_patient.txt"

# ------------------ PAM50 (NO NORMAL) ------------------
pam50_codes  <- c("BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("LumA","LumB","Her2","Basal-like")

# ------------------ STYLE ------------------
base_size   <- 17
axis_x_size <- 16
axis_y_size <- 16
axis_title  <- 16
strip_pt    <- 16

pt_size  <- 1.6
pt_alpha <- 0.35

tiff_dpi  <- 600
tiff_w_in <- 12
tiff_h_in <- 3.6

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

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 4)))
}

# ------------------ LOAD CLINICAL ------------------
clin <- read.delim(CLIN_FILE, check.names = FALSE, comment.char = "#") %>%
  dplyr::select(PATIENT_ID, SUBTYPE, SEX) %>%
  dplyr::filter(toupper(as.character(SEX)) == "FEMALE") %>%
  dplyr::mutate(
    SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels),
    SUBTYPE = factor(as.character(SUBTYPE), levels = pam50_labels)
  ) %>%
  dplyr::filter(!is.na(SUBTYPE))   # <- Normal automatically excluded

# ------------------ LOAD EXPRESSION ------------------
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

genes_needed <- c("CAV1","ATP13A3")

expr_long <- expr %>%
  dplyr::filter(Hugo_Symbol %in% genes_needed) %>%
  tidyr::pivot_longer(
    cols = -Hugo_Symbol,
    names_to = "SAMPLE_ID",
    values_to = "expr_raw"
  ) %>%
  dplyr::mutate(
    expr       = log2(as.numeric(expr_raw) + 1),
    PATIENT_ID = substr(SAMPLE_ID, 1, 12),
    gene       = Hugo_Symbol
  ) %>%
  dplyr::select(PATIENT_ID, gene, expr) %>%
  dplyr::filter(!is.na(expr))

# ------------------ WIDE + MERGE ------------------
dat <- expr_long %>%
  tidyr::pivot_wider(names_from = gene, values_from = expr) %>%
  dplyr::inner_join(
    clin %>% dplyr::select(PATIENT_ID, SUBTYPE),
    by = "PATIENT_ID"
  ) %>%
  dplyr::filter(!is.na(CAV1), !is.na(ATP13A3))

if (nrow(dat) == 0) {
  stop("No data left after merging expression and clinical tables.", call. = FALSE)
}

# ------------------ SPEARMAN PER SUBTYPE (for strip) ------------------
cor_df <- dat %>%
  dplyr::group_by(SUBTYPE) %>%
  dplyr::summarise(
    rho = cor(CAV1, ATP13A3, method = "spearman", use = "complete.obs"),
    p   = cor.test(CAV1, ATP13A3, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    strip = paste0(
      "<b>", as.character(SUBTYPE), "</b>",
      "<br><span style='font-weight:400;'>rho=", sprintf("%.2f", rho),
      "; p=", fmt_p(p),
      "</span>"
    )
  )

strip_vec <- setNames(cor_df$strip, as.character(cor_df$SUBTYPE))

# ------------------ PLOT ------------------
p <- ggplot(dat, aes(x = CAV1, y = ATP13A3)) +
  geom_point(size = pt_size, alpha = pt_alpha, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black") +
  facet_wrap(
    ~ SUBTYPE,
    nrow = 1,
    labeller = labeller(SUBTYPE = strip_vec)
  ) +
  labs(
    x = "**CAV1** (log2 RSEM + 1)",
    y = "**ATP13A3** (log2 RSEM + 1)"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    strip.text   = ggtext::element_markdown(size = strip_pt),
    axis.title.x = ggtext::element_markdown(size = axis_title),
    axis.title.y = ggtext::element_markdown(size = axis_title),
    axis.text.x  = element_text(size = axis_x_size),
    axis.text.y  = element_text(size = axis_y_size),
    panel.grid.minor = element_blank()
  )

print(p)
save_tiff(p, "TCGA_BRCA_PAM50_Spearman_CAV1_vs_ATP13A3_NO_NORMAL")


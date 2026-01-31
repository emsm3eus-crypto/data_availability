# ==========================================================
# TCGA BRCA — Spearman correlation of polyamine genes
# with MKI67 and MCM2 by PAM50 subtype (DepMap-style colors)
# - NO stars
# - p-values only if p < 0.05
# - p-values placed just above (+rho) or clearly below (-rho) bars
# - p-values rotated 90°
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
})

# ------------------ SETTINGS ------------------
WORKDIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
setwd(WORKDIR)

EXPR_FILE <- "data_mrna_seq_v2_rsem.txt"
CLIN_FILE <- "data_clinical_patient.txt"

poly_genes <- c(
  "ODC1","SRM","SMS",
  "ATP13A2","ATP13A3","ATP13A4",
  "SMOX","PAOX","SAT1"
)

pam50_codes  <- c("BRCA_Normal","BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("Normal","LumA","LumB","Her2","Basal-like")

# DepMap/your TCGA scripts palette
subtype_cols <- c(
  "LumA"       = "darkkhaki",
  "LumB"       = "lightpink",
  "Her2"       = "orange",
  "Basal-like" = "lightblue"
)

subtype_levels_plot <- c("LumA","LumB","Her2","Basal-like")

# DepMap-like sizes
base_size   <- 20
axis_x_size <- 20
axis_y_size <- 20
y_title_sz  <- 20

# DepMap-like bar look
bar_alpha <- 0.75

# ------------------ HELPERS ------------------
fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 1e-4, "<1e-4", formatC(p, format = "f", digits = 3)))
}

# place p-values relative to bar end
label_y_barend <- function(rho, pad_pos = 0.05, pad_neg = 0.20) {
  ifelse(
    is.na(rho),
    NA_real_,
    ifelse(rho >= 0,
           rho + pad_pos,   # just above positive bars
           rho - pad_neg)   # clearly below negative bars
  )
}

# ------------------ LOAD DATA ------------------
clin <- read.delim(CLIN_FILE, check.names = FALSE, comment.char = "#")
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% colnames(expr))

subtype_df <- clin %>%
  select(PATIENT_ID, SUBTYPE, SEX) %>%
  filter(toupper(as.character(SEX)) == "FEMALE") %>%
  mutate(SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels)) %>%
  filter(!is.na(SUBTYPE))

genes_needed <- c(poly_genes, "MKI67", "MCM2")

expr_long <- expr %>%
  filter(Hugo_Symbol %in% genes_needed) %>%
  pivot_longer(-Hugo_Symbol, names_to = "SAMPLE_ID", values_to = "expr_raw") %>%
  mutate(
    expr       = log2(as.numeric(expr_raw) + 1),
    PATIENT_ID = substr(SAMPLE_ID, 1, 12),
    gene       = Hugo_Symbol
  ) %>%
  select(PATIENT_ID, gene, expr)

expr_joined <- expr_long %>%
  inner_join(subtype_df, by = "PATIENT_ID") %>%
  filter(SUBTYPE != "Normal") %>%
  mutate(SUBTYPE = factor(as.character(SUBTYPE), levels = subtype_levels_plot))

wide_expr <- expr_joined %>%
  select(PATIENT_ID, SUBTYPE, gene, expr) %>%
  distinct() %>%
  pivot_wider(names_from = gene, values_from = expr)

# ------------------ CORRELATIONS ------------------
compute_corr <- function(target_gene) {
  
  wide_expr %>%
    pivot_longer(
      cols = all_of(poly_genes),
      names_to = "gene",
      values_to = "expr_gene"
    ) %>%
    group_by(SUBTYPE, gene) %>%
    summarise(
      rho = suppressWarnings(
        cor(expr_gene, .data[[target_gene]],
            method = "spearman",
            use = "pairwise.complete.obs")
      ),
      p = suppressWarnings(
        tryCatch(
          cor.test(expr_gene, .data[[target_gene]], method = "spearman")$p.value,
          error = function(e) NA_real_
        )
      ),
      n = sum(!is.na(expr_gene) & !is.na(.data[[target_gene]])),
      .groups = "drop"
    ) %>%
    mutate(
      gene    = factor(gene, levels = poly_genes),
      SUBTYPE = factor(SUBTYPE, levels = subtype_levels_plot),
      p_label = ifelse(!is.na(p) & p < 0.05, paste0("p=", fmt_p(p)), ""),
      y_p     = ifelse(!is.na(p) & p < 0.05, label_y_barend(rho), NA_real_)
    )
}

# ------------------ PLOTTING ------------------
make_corr_plot <- function(df, target_label) {
  
  y_lim <- range(df$rho, na.rm = TRUE) + c(-0.25, 0.25)
  
  ggplot(df, aes(x = gene, y = rho, fill = SUBTYPE)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.6) +
    geom_col(
      position = position_dodge(width = 0.85),
      width    = 0.75,
      color    = "black",
      alpha    = bar_alpha
    ) +
    geom_text(
      aes(y = y_p, label = p_label),
      position = position_dodge(width = 0.85),
      angle    = 90,
      size     = 6,
      vjust    = 0.5,
      hjust    = 0,
      na.rm    = TRUE
    ) +
    scale_fill_manual(values = subtype_cols, drop = FALSE, name = "PAM50 subtype") +
    scale_y_continuous(limits = y_lim) +
    labs(
      x = NULL,
      y = paste0("<b>Spearman correlation</b><br><span style='font-weight:400;'>", target_label, "</span>")
    ) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = axis_x_size),
      axis.text.y  = element_text(size = axis_y_size),
      axis.title.y = element_markdown(size = y_title_sz),
      panel.grid.minor = element_blank()
    )
}

# ------------------ RUN ------------------
corr_mki67 <- compute_corr("MKI67")
corr_mcm2  <- compute_corr("MCM2")

p_corr_mki67 <- make_corr_plot(corr_mki67, "MKI67 expression")
p_corr_mcm2  <- make_corr_plot(corr_mcm2,  "MCM2 expression")

print(p_corr_mki67)
print(p_corr_mcm2)


# ==========================================================
# TCGA BRCA (PanCan Atlas 2018) — Polyamine genes vs
# mTORC1 drivers + MAPK drivers (Spearman) — ONE FIGURE
# Adds an annotation strip UNDER the x-axis with group labels:
#   PI3K/AKT/TSC/mTORC1 block  |  MAPK block
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------- Settings --------------------
WORKDIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
setwd(WORKDIR)
EXPR_FILE <- "data_mrna_seq_v2_rsem.txt"

# Polyamine genes (top → bottom in heatmap)
POLY_GENES <- c(
  "ODC1","SRM","SMS",
  "ATP13A2","ATP13A3","ATP13A4",
  "SMOX","PAOX","SAT1"
)

# ---------------- mTORC1 signalling genes (LEFT -> RIGHT) ----------------
MTOR_GENES <- c(
  "PIK3CA",
  "AKT1",
  "TSC1",
  "TSC2",
  "RHEB",
  "MTOR"
)

# ---------------- MAPK signalling drivers (LEFT -> RIGHT) ----------------
MAPK_GENES <- c(
  "KRAS",
  "NRAS",
  "HRAS",
  "BRAF",
  "MAP2K1",
  "MAP2K2",
  "MAPK1",
  "MAPK3"
)

# ---- Pretty labels for plotting (ONLY affects axis text) ----
MAPK_LABELS <- c(
  KRAS   = "KRAS",
  NRAS   = "NRAS",
  HRAS   = "HRAS",
  BRAF   = "BRAF",
  MAP2K1 = "MEK1",
  MAP2K2 = "MEK2",
  MAPK1  = "ERK2",  # MAPK1 = ERK2
  MAPK3  = "ERK1"   # MAPK3 = ERK1
)

# Optional: nicer mTOR labels (edit if you want)
MTOR_LABELS <- c(
  PIK3CA = "PI3K",
  AKT1   = "AKT",
  TSC1   = "TSC1",
  TSC2   = "TSC2",
  RHEB   = "RHEB",
  MTOR   = "mTOR"
)

# ---------------- Combined signalling axis ----------------
SIGNAL_GENES <- c(MTOR_GENES, MAPK_GENES)

# ---------------- Load expression ----------------
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

genes_needed <- unique(c(POLY_GENES, SIGNAL_GENES))

expr_long <- expr %>%
  filter(Hugo_Symbol %in% genes_needed) %>%
  pivot_longer(
    -Hugo_Symbol,
    names_to  = "SAMPLE_ID",
    values_to = "expr_raw"
  ) %>%
  mutate(
    expr       = log2(as.numeric(expr_raw) + 1),
    PATIENT_ID = substr(SAMPLE_ID, 1, 12),
    gene       = Hugo_Symbol
  ) %>%
  filter(!is.na(expr))

# Collapse multiple samples per patient (median)
expr_wide <- expr_long %>%
  group_by(PATIENT_ID, gene) %>%
  summarise(expr = median(expr), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = expr)

# ---------------- Compute Spearman correlations ----------------
cor_pairs <- tidyr::expand_grid(
  Poly_Gene   = POLY_GENES,
  Signal_Gene = SIGNAL_GENES
)

results <- cor_pairs %>%
  rowwise() %>%
  mutate({
    x <- expr_wide[[Poly_Gene]]
    y <- expr_wide[[Signal_Gene]]
    ok <- complete.cases(x, y)
    
    if (sum(ok) >= 3) {
      ct  <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      rho <- unname(ct$estimate)
      p   <- ct$p.value
    } else {
      rho <- NA_real_
      p   <- NA_real_
    }
    
    tibble(rho = rho, p = p)
  }) %>%
  unnest(c(rho, p)) %>%
  ungroup() %>%
  mutate(
    sig_star   = ifelse(!is.na(p) & p < 0.05, "*", ""),
    Poly_Gene  = factor(Poly_Gene, levels = POLY_GENES),
    Signal_Gene = factor(Signal_Gene, levels = SIGNAL_GENES)
  )

print(results, n = Inf)

# ---------------- Heatmap ----------------
# Build a named vector for x-axis labels (mTOR + MAPK)
X_LABELS <- c(MTOR_LABELS, MAPK_LABELS)

heatmap_plot <- ggplot(results, aes(x = Signal_Gene, y = Poly_Gene, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sig_star)), size = 6) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = X_LABELS) +
  scale_y_discrete(limits = rev(POLY_GENES)) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(5.5, 5.5, 0, 5.5)
  )

# ---------------- Annotation strip UNDER x-axis ----------------
# Two blocks: mTORC1 (includes PI3K/AKT/TSC/RHEB/mTOR) and MAPK
anno_df <- tibble(
  group = c("PI3K/AKT/mTOR", "MAPK"),
  xmin  = c(0.5, length(MTOR_GENES) + 0.5),
  xmax  = c(length(MTOR_GENES) + 0.5, length(SIGNAL_GENES) + 0.5),
  xmid  = c((0.5 + (length(MTOR_GENES) + 0.5))/2,
            ((length(MTOR_GENES) + 0.5) + (length(SIGNAL_GENES) + 0.5))/2)
)

anno_plot <- ggplot(anno_df) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), fill = "grey90", color = "grey60") +
  geom_text(aes(x = xmid, y = 0.5, label = group), fontface = "bold", size = 5) +
  scale_x_continuous(limits = c(0.5, length(SIGNAL_GENES) + 0.5), expand = c(0, 0)) +
  theme_void(base_size = 16) +
  theme(
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

# Combine: heatmap on top, annotation strip below
final_plot <- heatmap_plot / anno_plot + plot_layout(heights = c(12, 1.4))

print(final_plot)


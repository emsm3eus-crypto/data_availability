# ==========================================
# Correlation heatmap: Polyamine genes vs MAPK drivers
# - Same stats + style as your script
# - X-axis shows MEK1/MEK2/ERK1/ERK2 labels (instead of MAP2K1/2, MAPK1/3)
# ==========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------- Settings --------------------
WORKDIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
setwd(WORKDIR)

EXPR_FILE <- "data_mrna_seq_v2_rsem.txt"

# Polyamine genes (top → bottom order in heatmap)
POLY_GENES <- c(
  "ODC1","SRM","SMS",
  "ATP13A2","ATP13A3","ATP13A4",
  "SMOX","PAOX","SAT1"
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

# ------------- Load expression ----------------
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

genes_needed <- unique(c(POLY_GENES, MAPK_GENES))

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

# ------------- Compute Spearman correlations -------------
cor_pairs <- tidyr::expand_grid(
  Poly_Gene = POLY_GENES,
  MAPK_Gene = MAPK_GENES
)

results <- cor_pairs %>%
  rowwise() %>%
  mutate({
    x <- expr_wide[[Poly_Gene]]
    y <- expr_wide[[MAPK_Gene]]
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
    sig_star  = ifelse(!is.na(p) & p < 0.05, "*", ""),
    Poly_Gene = factor(Poly_Gene, levels = POLY_GENES),
    MAPK_Gene = factor(MAPK_Gene, levels = MAPK_GENES)
  )

print(results, n = Inf)

# ---------------- Heatmap --------------------
heatmap_plot <- ggplot(results, aes(x = MAPK_Gene, y = Poly_Gene, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = paste0(sprintf("%.2f", rho), sig_star)),
    size = 7
  ) +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    name     = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = MAPK_LABELS) +
  scale_y_discrete(limits = rev(POLY_GENES)) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  )

print(heatmap_plot)



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

# ---------------- mTORC1 signalling genes (LEFT -> RIGHT as requested) ----------------
# Requested order (concept -> gene symbols used in TCGA):
# PI3K   -> PIK3CA
# AKT    -> AKT1
# TSC    -> TSC2
# RHEB   -> RHEB
# mTOR   -> MTOR
# RAPTOR -> RPTOR
# DEPTOR -> DEPTOR
# MLST8  -> MLST8
# PRAS40 -> AKT1S1

MTOR_GENES <- c(
  "PIK3CA",
  "AKT1",
  "TSC1", 
  "TSC2",
  "RHEB",
  "MTOR"    
)

# ------------- Load expression ----------------
expr <- read.delim(EXPR_FILE, check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

genes_needed <- unique(c(POLY_GENES, MTOR_GENES))

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
  MTOR_Gene = MTOR_GENES
)

results <- cor_pairs %>%
  rowwise() %>%
  mutate({
    x <- expr_wide[[Poly_Gene]]
    y <- expr_wide[[MTOR_Gene]]
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
    MTOR_Gene = factor(MTOR_Gene, levels = MTOR_GENES)
  )

print(results, n = Inf)

# ---------------- Heatmap --------------------
heatmap_plot <- ggplot(results, aes(x = MTOR_Gene, y = Poly_Gene, fill = rho)) +
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


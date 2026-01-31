# ==========================================
# GSEA NES HEATMAP (genes on Y, pathways on X)
# - Rows = genes (Synthesis | Transport | Catabolism)
# - Cols = pathways (fixed order LEFT → RIGHT)
# - Tile fill = NES (blue → white → red, midpoint 0)
# - NES number centered
# - Star (*) in top-right of tile if FDR q < 0.25
# ==========================================

suppressPackageStartupMessages({
  pkgs <- c("readr","dplyr","stringr","tidyr","ggplot2","purrr")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# -----------------------------
# PATHS
# -----------------------------
BASE_DIR <- "C:/Users/u0171728/Downloads/GSEA_output"
setwd(BASE_DIR)

folders <- list.dirs(BASE_DIR, full.names = TRUE, recursive = FALSE)

# -----------------------------
# GROUPS / GENE ORDER (TOP → BOTTOM)
# -----------------------------
groups <- list(
  "Synthesis"  = c("ODC1", "SRM", "SMS"),
  "Transport"  = c("ATP13A2", "ATP13A3", "ATP13A4"),
  "Catabolism" = c("SMOX", "PAOX", "SAT1")
)

genes_keep <- c(
  groups$Synthesis,
  groups$Transport,
  groups$Catabolism
)

# -----------------------------
# PATHWAYS (LEFT → RIGHT)
# -----------------------------
pathways_keep <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP"
)

pathway_labels <- c(
  HALLMARK_G2M_CHECKPOINT          = "G2M_CHECKPOINT",
  HALLMARK_MYC_TARGETS_V1          = "MYC_TARGETS_V1",
  HALLMARK_E2F_TARGETS             = "E2F_TARGETS",
  HALLMARK_MTORC1_SIGNALING        = "MTORC1_SIGNALING",
  HALLMARK_PI3K_AKT_MTOR_SIGNALING = "PI3K_AKT_MTOR_SIGNALING",
  HALLMARK_KRAS_SIGNALING_UP       = "KRAS_SIGNALING_UP"
)

# -----------------------------
# STAR RULE
# -----------------------------
Q_STAR <- 0.25

# -----------------------------
# STYLE (matched to your correlation heatmap)
# -----------------------------
BASE_SIZE <- 16
NUM_SIZE  <- 7
STAR_SIZE <- 7

# Star placement (top-right of tile)
STAR_NUDGE_X <- 0.38
STAR_NUDGE_Y <- 0.38

# -----------------------------
# READ ONE GSEA FOLDER
# -----------------------------
read_gsea_folder <- function(folder_path) {
  gene <- str_replace(basename(folder_path), "_preranked_TCGA$", "")
  
  pos <- list.files(folder_path, pattern = "_pos_.*tsv$", full.names = TRUE)
  neg <- list.files(folder_path, pattern = "_neg_.*tsv$", full.names = TRUE)
  
  read_one <- function(f) {
    if (length(f) == 0) return(NULL)
    suppressMessages(read_tsv(f[1], show_col_types = FALSE))
  }
  
  df <- bind_rows(read_one(pos), read_one(neg))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    rename(pathway = NAME, nes = NES, fdr = `FDR q-val`) %>%
    mutate(gene = gene) %>%
    filter(pathway %in% pathways_keep) %>%
    select(gene, pathway, nes, fdr)
}

# -----------------------------
# LOAD ALL RESULTS
# -----------------------------
all_res <- map_dfr(folders, read_gsea_folder)
if (nrow(all_res) == 0) stop("No GSEA files read.")

all_res_full <- all_res %>%
  filter(gene %in% genes_keep) %>%
  complete(
    gene    = genes_keep,
    pathway = pathways_keep,
    fill    = list(nes = NA_real_, fdr = NA_real_)
  ) %>%
  mutate(
    gene          = factor(gene, levels = rev(genes_keep)), # top → bottom
    pathway_label = factor(pathway_labels[pathway],
                           levels = pathway_labels[pathways_keep])
  )

# -----------------------------
# PREP LABELS
# -----------------------------
hm_df <- all_res_full %>%
  mutate(
    nes_txt = ifelse(is.na(nes), NA_character_, sprintf("%.2f", nes)),
    star    = ifelse(!is.na(fdr) & fdr < Q_STAR, "*", NA_character_)
  )

NES_LIMIT <- max(abs(hm_df$nes), na.rm = TRUE)
if (!is.finite(NES_LIMIT) || NES_LIMIT == 0) NES_LIMIT <- 1

# -----------------------------
# HEATMAP
# -----------------------------
p_nes <- ggplot(
  hm_df,
  aes(x = pathway_label, y = gene, fill = nes)
) +
  geom_tile(color = "white") +
  
  # NES number (center)
  geom_text(
    aes(label = nes_txt),
    size = NUM_SIZE,
    na.rm = TRUE,
    color = "black"
  ) +
  
  # Star in top-right of tile
  geom_text(
    aes(label = star),
    size = STAR_SIZE,
    na.rm = TRUE,
    color = "black",
    hjust = 1,
    vjust = 1,
    position = position_nudge(
      x = STAR_NUDGE_X,
      y = STAR_NUDGE_Y
    )
  ) +
  
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-NES_LIMIT, NES_LIMIT),
    name     = "NES"
  ) +
  
  labs(x = NULL, y = NULL) +
  
  theme_bw(base_size = BASE_SIZE) +
  theme(
    axis.text.x  = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  )

print(p_nes)

# -----------------------------
# SAVE
# -----------------------------
out_dir <- file.path(BASE_DIR, "NES_heatmap_genesY_pathwaysX")
dir.create(out_dir, showWarnings = FALSE)

ggsave(
  file.path(out_dir, "GSEA_NES_heatmap_genesY_pathwaysX.png"),
  p_nes, width = 12.5, height = 5.5, dpi = 300
)

ggsave(
  file.path(out_dir, "GSEA_NES_heatmap_genesY_pathwaysX.pdf"),
  p_nes, width = 12.5, height = 5.5
)

cat("\nLegend: '*' shown only if FDR q < ", Q_STAR,
    " (star is top-right of tile)\n", sep = "")

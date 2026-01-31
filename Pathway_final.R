# ==========================================
# GSEA dot plots (NES + FDR q-values)
# - Grid behind
# - q-values under/over dots (auto-stagger when dots are close)
# - No dodging: dots on same line per pathway
# - Fixed grayscale mapping (gene1 black, gene2 white, gene3 grey)
# - Legend ON
#
# ORDER (top->bottom):
#   G2M, MYC V1, E2F, then the rest
#
# UPDATED:
# - Bigger q-value text
# - Bigger pathway labels
# - Bigger legend (gene names)
# ==========================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# ---- Settings ----
BASE_DIR <- "C:/Users/u0171728/Downloads/GSEA_output"
setwd(BASE_DIR)

folders <- list.dirs(BASE_DIR, full.names = TRUE, recursive = FALSE)

groups <- list(
  "Synthesis (ODC1/SRM/SMS)"    = c("ODC1", "SRM", "SMS"),
  "ATP13 (ATP13A2/3/4)"         = c("ATP13A2", "ATP13A3", "ATP13A4"),
  "Catabolism (SMOX/PAOX/SAT1)" = c("SMOX", "PAOX", "SAT1")
)

# Pathways to extract (GSEA Hallmark names) + desired order TOP->BOTTOM
pathways_keep <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP"
)

pathway_labels <- setNames(gsub("^HALLMARK_", "", pathways_keep), pathways_keep)

# If dots are closer than this in NES (within a pathway), stagger labels above/below
CLOSE_THRESH <- 0.40

# ---- TEXT SIZES (BIGGER) ----
BASE_SIZE       <- 17   # overall base size (was 12)
DOT_SIZE        <- 4.2  # dot size
DOT_STROKE      <- 1.1  # dot outline thickness
Q_TEXT_SIZE     <- 4.7  # q-value label size (was 3.4)
PATHWAY_TEXT_SZ <- 16   # pathway labels on y-axis (was 10)
XAXIS_TEXT_SZ   <- 16   # x-axis tick labels
AXIS_TITLE_SZ   <- 16   # x-axis title
LEGEND_TITLE_SZ <- 16   # legend title
LEGEND_TEXT_SZ  <- 16   # gene names in legend (was 10)

# -----------------------------
# Read one GSEA folder
# -----------------------------
read_gsea_folder <- function(folder_path) {
  gene <- str_replace(basename(folder_path), "_preranked_TCGA$", "")
  
  pos_tsv <- list.files(folder_path, pattern = "^gsea_report_for_.*_pos_.*\\.tsv$", full.names = TRUE)
  neg_tsv <- list.files(folder_path, pattern = "^gsea_report_for_.*_neg_.*\\.tsv$", full.names = TRUE)
  
  read_one <- function(f) {
    if (length(f) == 0) return(NULL)
    suppressMessages(read_tsv(f[1], show_col_types = FALSE))
  }
  
  df <- bind_rows(read_one(pos_tsv), read_one(neg_tsv))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    rename(
      pathway = NAME,
      nes     = NES,
      fdr     = `FDR q-val`
    ) %>%
    mutate(gene = gene) %>%
    filter(pathway %in% pathways_keep) %>%
    select(gene, pathway, nes, fdr)
}

# -----------------------------
# Read all folders
# -----------------------------
all_res <- map_dfr(folders, read_gsea_folder)

if (nrow(all_res) == 0) {
  stop("No GSEA report TSVs were read. Check folder structure and file names.")
}

# Ensure all gene x pathway combinations exist (so plots don't drop missing)
all_res_full <- all_res %>%
  complete(gene, pathway, fill = list(nes = NA_real_, fdr = NA_real_)) %>%
  mutate(pathway_label = pathway_labels[pathway])

# -----------------------------
# Plot function (one plot per group)
# -----------------------------
make_group_plot <- function(group_name, genes_vec, df) {
  
  # Fixed grayscale per gene position:
  # gene1 = black, gene2 = white, gene3 = grey
  fill_values <- setNames(c("black", "white", "grey60"), genes_vec)
  
  d <- df %>%
    filter(gene %in% genes_vec) %>%
    mutate(
      gene = factor(gene, levels = genes_vec),
      # critical: levels here control TOP->BOTTOM order on the plot
      pathway_label = factor(pathway_label, levels = pathway_labels[pathways_keep]),
      fdr_label = case_when(
        is.na(fdr)  ~ "",
        fdr < 0.001 ~ "q<0.001",
        TRUE        ~ sprintf("q=%.2f", fdr)
      )
    ) %>%
    group_by(pathway_label) %>%
    mutate(
      nes_rank = dense_rank(nes),
      min_nonzero_sep = {
        v <- sort(nes)
        if (length(v) < 2 || all(is.na(v))) NA_real_ else suppressWarnings(min(abs(diff(v)), na.rm = TRUE))
      },
      close = if_else(is.na(min_nonzero_sep), FALSE, min_nonzero_sep < CLOSE_THRESH),
      # if close: alternate above/below
      # else: all below
      vjust_lab = if_else(close, if_else(nes_rank %% 2 == 1, 2.0, -1.0), 2.0)
    ) %>%
    ungroup()
  
  ggplot(d, aes(x = nes, y = pathway_label, fill = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_point(
      shape  = 21,
      size   = DOT_SIZE,
      stroke = DOT_STROKE,
      color  = "black",
      na.rm  = TRUE
    ) +
    geom_text(
      aes(label = fdr_label, vjust = vjust_lab),
      size  = Q_TEXT_SIZE,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_manual(values = fill_values, name = "Gene") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = LEGEND_TITLE_SZ, face = "bold"),
      legend.text  = element_text(size = LEGEND_TEXT_SZ),
      
      axis.text.y  = element_text(size = PATHWAY_TEXT_SZ),
      axis.text.x  = element_text(size = XAXIS_TEXT_SZ),
      axis.title.x = element_text(size = AXIS_TITLE_SZ),
      
      # Grid behind
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# Make and save the 3 plots
# -----------------------------
out_dir <- file.path(BASE_DIR, "NES_group_dotplots_grid_bw_fdr_staggered_BIGTEXT")
dir.create(out_dir, showWarnings = FALSE)

plots <- imap(groups, ~make_group_plot(.y, .x, all_res_full))

iwalk(plots, ~ggsave(
  filename = file.path(out_dir, paste0(str_replace_all(.y, "[^A-Za-z0-9]+", "_"), "_NES_dot_grid_bw_fdr_BIGTEXT.png")),
  plot = .x,
  width = 9,
  height = 4.8,
  dpi = 300
))

print(plots[[1]])
print(plots[[2]])
print(plots[[3]])

cat("\nSaved plots to:\n", out_dir, "\n")



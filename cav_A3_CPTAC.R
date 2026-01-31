# ==========================================================
# CPTAC BRCA 2020 (RAW protein) — Spearman correlation: CAV1 vs ATP13A3
# per PAM50 subtype (NO Normal; facets naast elkaar) — TCGA-style layout
#
# - Removed Normal subtype
# - CAV1 and ATP13A3 in **bold**
# - Facet subtitle shows: rho + p (no n)
# - RAW protein quantification (data_protein_quantification.txt)
# - SAVE: TIFF 600 dpi (LZW) to Downloads using ragg
# - DISPLAY NAMES: Her2 and Basal-like (instead of HER2 and Basal)
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","ggplot2","ggtext","stringr","ragg")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ------------------ PATHS ------------------
DATA_DIR <- "C:/Users/u0171728/Downloads/brca_cptac_2020/brca_cptac_2020"
OUT_DIR  <- "C:/Users/u0171728/Downloads"
setwd(DATA_DIR)

PROT_FILE <- "data_protein_quantification.txt"
CLIN_FILE <- "data_clinical_sample.txt"

# ------------------ SUBTYPES (NO NORMAL) ------------------
# Keep internal labels stable for parsing/matching:
pam50_levels <- c("LumA","LumB","HER2","Basal")  # no Normal

# Display labels for facets:
pam50_display_levels <- c("LumA","LumB","Her2","Basal-like")

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

pick_first <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

# ==========================================================
# LOAD CLINICAL (sample-level)
# ==========================================================
clin_s <- read.delim(CLIN_FILE, check.names = FALSE, comment.char = "#")

sample_col_clin <- pick_first(clin_s, c("SAMPLE_ID","Sample_ID","sampleId","ID"))
subtype_col <- pick_first(clin_s, c(
  "PAM50","PAM50_SUBTYPE","SUBTYPE","SUBTYPE_PAM50",
  "PAM50.mRNA","PAM50 (mRNA)","PAM50 subtype","pam50"
))

if (is.na(sample_col_clin)) stop("Could not find SAMPLE_ID column in data_clinical_sample.txt.", call. = FALSE)
if (is.na(subtype_col))     stop("Could not find a PAM50/SUBTYPE column in data_clinical_sample.txt.", call. = FALSE)

clin <- clin_s %>%
  transmute(
    SAMPLE_ID   = as.character(.data[[sample_col_clin]]),
    SUBTYPE_RAW = as.character(.data[[subtype_col]])
  ) %>%
  filter(!is.na(SUBTYPE_RAW), SUBTYPE_RAW != "") %>%
  mutate(
    SUBTYPE = SUBTYPE_RAW,
    SUBTYPE = str_replace_all(SUBTYPE, fixed("Basal-like"), "Basal"),
    SUBTYPE = str_replace_all(SUBTYPE, fixed("HER2-enriched"), "HER2"),
    SUBTYPE = str_replace_all(SUBTYPE, fixed("Her2"), "HER2"),
    SUBTYPE = factor(SUBTYPE, levels = c("Normal", pam50_levels))  # allow parsing
  ) %>%
  filter(!is.na(SUBTYPE)) %>%
  filter(as.character(SUBTYPE) %in% pam50_levels) %>%              # DROP NORMAL
  mutate(SUBTYPE = factor(as.character(SUBTYPE), levels = pam50_levels))

# ==========================================================
# LOAD RAW PROTEIN + EXTRACT ATP13A3 & CAV1
# ==========================================================
prot <- read.delim(PROT_FILE, check.names = FALSE, comment.char = "#")

feat_col <- pick_first(prot, c("Composite.Element.REF","COMPOSITE.ELEMENT.REF","Composite Element REF","feature","ID"))
if (is.na(feat_col)) stop("Could not find Composite.Element.REF column in data_protein_quantification.txt.", call. = FALSE)

prot2 <- prot %>%
  mutate(
    feature = as.character(.data[[feat_col]]),
    gene    = stringr::str_extract(feature, "^[^|]+")
  )

genes_needed <- c("CAV1","ATP13A3")
missing_genes <- setdiff(genes_needed, unique(prot2$gene))
if (length(missing_genes) > 0) {
  stop("Missing genes in protein file: ", paste(missing_genes, collapse = ", "), call. = FALSE)
}

non_sample_cols <- unique(c(feat_col, "feature", "gene"))
sample_cols <- setdiff(colnames(prot2), non_sample_cols)

prot_long <- prot2 %>%
  filter(gene %in% genes_needed) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "SAMPLE_ID",
    values_to = "expr_raw"
  ) %>%
  transmute(
    SAMPLE_ID = as.character(SAMPLE_ID),
    gene      = as.character(gene),
    expr      = suppressWarnings(as.numeric(expr_raw))
  ) %>%
  filter(!is.na(expr))

# ==========================================================
# WIDE + MERGE
# ==========================================================
dat <- prot_long %>%
  pivot_wider(names_from = gene, values_from = expr) %>%
  inner_join(clin %>% select(SAMPLE_ID, SUBTYPE), by = "SAMPLE_ID") %>%
  filter(!is.na(CAV1), !is.na(ATP13A3))

if (nrow(dat) == 0) {
  stop("No data left after merging protein and clinical tables (and removing Normal).", call. = FALSE)
}

# ==========================================================
# DISPLAY LABELS (Her2, Basal-like) + ORDER
# ==========================================================
dat <- dat %>%
  mutate(
    SUBTYPE_LABEL = recode(
      as.character(SUBTYPE),
      "HER2"  = "Her2",
      "Basal" = "Basal-like",
      .default = as.character(SUBTYPE)
    ),
    SUBTYPE_LABEL = factor(SUBTYPE_LABEL, levels = pam50_display_levels)
  )

# ==========================================================
# SPEARMAN PER DISPLAY SUBTYPE (for strip)
# ==========================================================
cor_df <- dat %>%
  group_by(SUBTYPE_LABEL) %>%
  summarise(
    rho = tryCatch(
      cor(CAV1, ATP13A3, method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    ),
    p = tryCatch(
      cor.test(CAV1, ATP13A3, method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    strip = paste0(
      "<b>", as.character(SUBTYPE_LABEL), "</b>",
      "<br><span style='font-weight:400;'>rho=", sprintf("%.2f", rho),
      "; p=", fmt_p(p),
      "</span>"
    )
  )

strip_vec <- setNames(cor_df$strip, as.character(cor_df$SUBTYPE_LABEL))

# ==========================================================
# PLOT (TCGA-style layout)
# ==========================================================
p <- ggplot(dat, aes(x = CAV1, y = ATP13A3)) +
  geom_point(size = pt_size, alpha = pt_alpha, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black") +
  facet_wrap(
    ~ SUBTYPE_LABEL,
    nrow = 1,
    labeller = labeller(SUBTYPE_LABEL = strip_vec)
  ) +
  labs(
    x = "**CAV1** (protein)",
    y = "**ATP13A3** (protein)"
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
save_tiff(p, "CPTAC_BRCA2020_RAW_PAM50_Spearman_CAV1_vs_ATP13A3_NO_NORMAL")

cat(
  "\nSaved TIFF to:\n",
  file.path(OUT_DIR, "CPTAC_BRCA2020_RAW_PAM50_Spearman_CAV1_vs_ATP13A3_NO_NORMAL.tiff"),
  "\n"
)


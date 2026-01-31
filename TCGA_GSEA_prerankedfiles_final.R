# ==========================================================
# TCGA BRCA — GSEA Preranked (.rnk): High vs Low of ODC1/SMOX/ATP13A3
# Rank metric: limma moderated t (High vs Low)
# Cohort: TCGA-BRCA primary tumors (TCGA code "01"), FEMALE
# Output: ./GSEA_preranked/{ODC1,SMOX,ATP13A3}_High_vs_Low_overall.rnk
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","readr","limma","matrixStats")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  lapply(pkgs, library, character.only = TRUE)
})

# -------- Settings --------
WORKDIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
setwd(WORKDIR)

DRIVERS <- c(
  "ODC1", "SMS", "SRM",
  "ATP13A2", "ATP13A3", "ATP13A4",
  "SMOX", "PAOX", "SAT1"
)

KEEP_SEX  <- c("Female","FEMALE","female")

# PAM50 codes used by cBioPortal export
pam50_codes  <- c("BRCA_Normal","BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("Normal","LumA","LumB","Her2","Basal")

dir.create("GSEA_preranked", showWarnings = FALSE, recursive = TRUE)

message("[1/6] Reading inputs...")

clin <- read.delim("data_clinical_patient.txt", check.names = FALSE, comment.char = "#")
expr <- read.delim("data_mrna_seq_v2_rsem.txt", check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

# -------- Clinical subset --------
message("[2/6] Preparing clinical table...")
clin2 <- clin %>%
  transmute(
    PATIENT_ID = PATIENT_ID,
    SUBTYPE    = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels),
    SEX        = SEX
  ) %>%
  filter(!is.na(SUBTYPE), SEX %in% KEEP_SEX) %>%         # females, valid subtype
  filter(SUBTYPE != "Normal") %>%                        # drop Normal-like
  distinct(PATIENT_ID, .keep_all = TRUE)

# -------- Build clean genes × patients matrix --------
message("[3/6] Building expression matrix (genes x patients)...")

# Keep only TCGA sample columns
sample_cols <- setdiff(names(expr), "Hugo_Symbol")

# Primary tumor samples only (TCGA barcode positions 14–15 == "01")
samples_primary <- sample_cols[substr(sample_cols, 14, 15) == "01"]

# Preserve gene symbols BEFORE numeric conversion
gene_symbols <- as.character(expr$Hugo_Symbol)
X <- as.matrix(expr[, samples_primary, drop = FALSE])
rownames(X) <- gene_symbols
mode(X) <- "numeric"                       # force numeric

# Transform
X <- log2(X + 1)

# Collapse replicate samples to patient-level medians (first 12 chars)
patients <- substr(samples_primary, 1, 12)
idx_list <- split(seq_along(patients), patients)
X_pat <- sapply(idx_list, function(cols) matrixStats::rowMedians(X[, cols, drop = FALSE], na.rm = TRUE))
colnames(X_pat) <- names(idx_list)         # patient IDs as column names

# Keep patients present in clinical table
common_pat <- intersect(colnames(X_pat), clin2$PATIENT_ID)
X_pat <- X_pat[, common_pat, drop = FALSE]
clin2 <- clin2 %>% filter(PATIENT_ID %in% common_pat)

# Remove rows with empty/NA symbols and collapse duplicate symbols (median across duplicates)
keep_rows <- !is.na(rownames(X_pat)) & nchar(rownames(X_pat)) > 0
X_pat <- X_pat[keep_rows, , drop = FALSE]
X_pat <- limma::avereps(X_pat, ID = rownames(X_pat), FUN = median)  # ensures unique HGNC rownames

cat(sprintf("Matrix ready: %d genes x %d patients\n", nrow(X_pat), ncol(X_pat)))
cat("Example symbols: ", paste(head(rownames(X_pat)), collapse = ", "), "\n")

# -------- Helper: make one preranked file (High vs Low of driver) --------
make_rnk <- function(driver_symbol, X_mat, outdir = "GSEA_preranked") {
  if (!(driver_symbol %in% rownames(X_mat))) {
    warning("Driver gene ", driver_symbol, " not found in matrix; skipping.")
    return(invisible(NULL))
  }
  
  drv <- X_mat[driver_symbol, ]
  cutoff <- median(drv, na.rm = TRUE)
  group  <- factor(ifelse(drv >= cutoff, "High", "Low"), levels = c("Low","High"))
  
  # Design: High vs Low
  design <- model.matrix(~ group)  # Intercept + groupHigh
  fit <- limma::lmFit(X_mat, design)
  fit <- limma::eBayes(fit, robust = TRUE)
  
  tt <- limma::topTable(fit, coef = "groupHigh", number = Inf, sort.by = "none")
  
  # Build (GeneSymbol, t-stat) table explicitly, DO NOT rely on rownames later
  rnk <- data.frame(
    gene = rownames(tt),
    rank = tt$t,
    stringsAsFactors = FALSE
  )
  
  # Deduplicate symbols if any remain (keep max |t|)
  rnk <- rnk[order(-abs(rnk$rank)), ]
  rnk <- rnk[!duplicated(rnk$gene), ]
  rnk <- rnk[order(-rnk$rank), ]
  
  out_file <- file.path(outdir, paste0(driver_symbol, "_High_vs_Low_overall.rnk"))
  write.table(rnk, file = out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  message(sprintf("[.rnk] %s | Ngenes=%d | High/Low=%d/%d | cutoff=%.3f",
                  basename(out_file), nrow(rnk), sum(group=="High"), sum(group=="Low"), cutoff))
}

# -------- Make the three preranked files --------
message("[4/6] Generating preranked lists...")
for (drv in DRIVERS) make_rnk(drv, X_pat, outdir = "GSEA_preranked")

message("[5/6] Quick sanity peek at a file header:")
cat(readLines(file.path("GSEA_preranked","ODC1_High_vs_Low_overall.rnk"), n = 5), sep = "\n")

message("[6/6] Done. In GSEA (GSEAPreranked) use:")
message(" - Ranked list: *_High_vs_Low_overall.rnk")
message(" - Gene sets: e.g., h.all.v2025.1.Hs.symbols.gmt")
message(" - Permutations: 1000 | Collapse: No_Collapse | Chip: (blank)")
message(" - If pruning error: set Min size=5, Max size=5000 (Basic fields).")

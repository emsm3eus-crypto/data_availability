# ==========================================================
# TCGA BRCA — Continuous Cox (per 1 log2 unit)
# Genes: ODC1, SMS, SRM, ATP13A2/3/4, SMOX, PAOX, SAT1
# Output: Forest plot (HR ± 95% CI, raw p-values only)
# - Female, primary tumors only, excludes Normal subtype
# - Expression = log2(RSEM + 1); per-patient median across primary samples
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","survival","broom","matrixStats","ggplot2","forcats","limma")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ---------- Settings ----------
WORKDIR <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
setwd(WORKDIR)

GENES    <- c("ODC1","SRM","SMS","ATP13A2","ATP13A3","ATP13A4","SMOX","PAOX","SAT1")
ENDPOINT <- "OS"
KEEP_SEX <- c("Female","FEMALE","female")

pam50_codes  <- c("BRCA_Normal","BRCA_LumA","BRCA_LumB","BRCA_Her2","BRCA_Basal")
pam50_labels <- c("Normal","LumA","LumB","Her2","Basal")

# ---------- Load ----------
clin <- read.delim("data_clinical_patient.txt", check.names = FALSE, comment.char = "#")
expr <- read.delim("data_mrna_seq_v2_rsem.txt", check.names = FALSE)
stopifnot("Hugo_Symbol" %in% names(expr))

# ---------- Clinical (OS) ----------
status_col <- paste0(ENDPOINT, "_STATUS")
months_col <- paste0(ENDPOINT, "_MONTHS")
if (!all(c(status_col, months_col) %in% names(clin))) stop("Missing OS columns in clinical file.")

clin2 <- clin %>%
  transmute(
    PATIENT_ID,
    SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels),
    SEX,
    time  = suppressWarnings(as.numeric(.data[[months_col]])),
    event = case_when(
      is.na(.data[[status_col]]) ~ NA_real_,
      grepl("^1", .data[[status_col]]) |
        grepl("DECEASED|DEAD|PROGRESSION|RECURRENCE", .data[[status_col]], TRUE) ~ 1,
      grepl("^0", .data[[status_col]]) |
        grepl("LIVING|ALIVE|CENSORED|DISEASE_FREE", .data[[status_col]], TRUE) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(SUBTYPE), SUBTYPE != "Normal", SEX %in% KEEP_SEX,
         !is.na(time), !is.na(event), time >= 0)

# ---------- Expression matrix (primary tumors only) ----------
sample_cols     <- setdiff(names(expr), "Hugo_Symbol")
samples_primary <- sample_cols[substr(sample_cols, 14, 15) == "01"]

X <- as.matrix(expr[, samples_primary, drop = FALSE])
rownames(X) <- expr$Hugo_Symbol
mode(X) <- "numeric"
X <- log2(X + 1)

patients <- substr(samples_primary, 1, 12)
idx_list <- split(seq_along(patients), patients)

X_pat <- sapply(
  idx_list,
  function(cols) matrixStats::rowMedians(X[, cols, drop = FALSE], na.rm = TRUE)
)
colnames(X_pat) <- names(idx_list)

common_pat <- intersect(colnames(X_pat), clin2$PATIENT_ID)
X_pat  <- X_pat[, common_pat, drop = FALSE]
clin2  <- clin2 %>% filter(PATIENT_ID %in% common_pat)

X_pat <- X_pat[!is.na(rownames(X_pat)) & nzchar(rownames(X_pat)), , drop = FALSE]
X_pat <- limma::avereps(X_pat, ID = rownames(X_pat), FUN = median)

# ---------- Continuous Cox per gene ----------
present <- GENES[GENES %in% rownames(X_pat)]

cox_list <- lapply(present, function(g) {
  d <- tibble(
    PATIENT_ID = colnames(X_pat),
    expr = as.numeric(X_pat[g, colnames(X_pat)])
  ) %>%
    inner_join(clin2, by = "PATIENT_ID") %>%
    filter(!is.na(expr), !is.na(time), !is.na(event))
  
  if (nrow(d) < 30 || length(unique(d$event)) < 2 || sd(d$expr) == 0) {
    return(tibble(gene = g, HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_))
  }
  
  fit <- coxph(Surv(time, event) ~ expr, data = d)
  
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "expr") %>%
    transmute(gene = g, HR = estimate, LCL = conf.low, UCL = conf.high, p = p.value)
})

cox_tbl <- bind_rows(cox_list) %>%
  mutate(gene = factor(gene, levels = present))

# ---------- Forest plot ----------
plot_df <- cox_tbl %>%
  mutate(
    fill_flag = case_when(
      !is.na(p) & p < 0.05 & !is.na(HR) & HR > 1 ~ "sig_HR_gt_1",
      !is.na(p) & p < 0.05 & !is.na(HR) & HR < 1 ~ "sig_HR_lt_1",
      TRUE                                       ~ "ns"
    ),
    label = ifelse(
      is.na(HR),
      "insufficient data",
      sprintf(
        "HR %.2f (%.2f–%.2f)\n p=%s",
        HR, LCL, UCL,
        ifelse(is.na(p), "NA", formatC(p, format="e", digits=2))
      )
    )
  )

p_forest <- ggplot(plot_df, aes(x = HR, y = forcats::fct_rev(gene))) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.18, na.rm = TRUE) +
  geom_point(aes(fill = fill_flag), shape = 21, size = 3, color = "black", na.rm = TRUE) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4), labels = c("0.5", "1", "2", "4")) +
  scale_fill_manual(values = c(
    "sig_HR_gt_1" = "firebrick",
    "sig_HR_lt_1" = "forestgreen",
    "ns"          = "#4e79a7"
  )) +
  geom_text(aes(x = Inf, label = label), hjust = -0.05, size = 6) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Hazard Ratio\nlog2(RSEM+1) expression",
    y = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    plot.margin     = margin(10, 220, 10, 10)
  )

print(p_forest)


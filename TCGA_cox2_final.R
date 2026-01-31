
# ==========================================================
# TCGA BRCA — Multivariable Cox (per 1 log2 unit), adjusted
# Adjust for: Age, Stage, PAM50 subtype (as available)
# Genes: ODC1, SMS, SRM, ATP13A2/3/4, SMOX, PAOX, SAT1
# Output: forest plot of adjusted HR for expr term (raw p-values)
# ==========================================================

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tidyr","survival","broom","limma","matrixStats",
            "ggplot2","forcats","stringr","purrr")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ---------- Settings ----------
WORKDIR  <- "C:/Users/u0171728/Downloads/brca_tcga_pan_can_atlas_2018 (1)/brca_tcga_pan_can_atlas_2018"
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

# ---------- Survival columns ----------
status_col <- paste0(ENDPOINT, "_STATUS")
months_col <- paste0(ENDPOINT, "_MONTHS")
if (!all(c(status_col, months_col) %in% names(clin))) stop("Missing survival columns in clinical file.")

# ---------- Detect covariate columns ----------
age_col <- names(clin)[grepl("^AGE$|AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS|AGE_AT_DIAGNOSIS",
                             names(clin), ignore.case = TRUE)][1]

stage_col <- names(clin)[grepl("AJCC|PATHOLOGIC_STAGE|CLINICAL_STAGE",
                               names(clin), ignore.case = TRUE)][1]

norm_stage <- function(x) {
  if (is.null(x)) return(factor(NA))
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^STAGE\\s*", "", x)
  x <- gsub("\\s*[ABC]$", "", x)
  x <- gsub("\\s*\\(.*\\)$", "", x)
  x <- dplyr::case_when(
    grepl("I\\b",   x) ~ "I",
    grepl("II\\b",  x) ~ "II",
    grepl("III\\b", x) ~ "III",
    grepl("IV\\b",  x) ~ "IV",
    TRUE ~ NA_character_
  )
  factor(x, levels = c("I","II","III","IV"), ordered = TRUE)
}

# ---------- Clinical prep ----------
clin2 <- clin %>%
  transmute(
    PATIENT_ID,
    SUBTYPE = factor(SUBTYPE, levels = pam50_codes, labels = pam50_labels),
    SEX,
    time  = suppressWarnings(as.numeric(.data[[months_col]])),
    event = case_when(
      is.na(.data[[status_col]]) ~ NA_real_,
      grepl("^1", .data[[status_col]]) |
        grepl("DECEASED|DEAD|PROGRESSION|RECURRENCE", .data[[status_col]], ignore.case = TRUE) ~ 1,
      grepl("^0", .data[[status_col]]) |
        grepl("LIVING|ALIVE|CENSORED|DISEASE_FREE", .data[[status_col]], ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_
    ),
    age_raw   = if (!is.na(age_col)) suppressWarnings(as.numeric(.data[[age_col]])) else NA_real_,
    stage_raw = if (!is.na(stage_col)) as.character(.data[[stage_col]]) else NA_character_
  ) %>%
  filter(!is.na(SUBTYPE), SUBTYPE != "Normal",
         SEX %in% KEEP_SEX,
         !is.na(time), !is.na(event), time >= 0) %>%
  mutate(
    age   = age_raw,
    stage = norm_stage(stage_raw)
  )

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

# Align with clinical
common_pat <- intersect(colnames(X_pat), clin2$PATIENT_ID)
X_pat <- X_pat[, common_pat, drop = FALSE]
clin2  <- clin2 %>% filter(PATIENT_ID %in% common_pat)

# Ensure unique gene symbols
X_pat <- X_pat[!is.na(rownames(X_pat)) & nzchar(rownames(X_pat)), , drop = FALSE]
X_pat <- limma::avereps(X_pat, ID = rownames(X_pat), FUN = median)

# ---------- Analysis frame ----------
present <- GENES[GENES %in% rownames(X_pat)]

dat <- tidyr::expand_grid(PATIENT_ID = colnames(X_pat), gene = present) %>%
  mutate(expr = mapply(function(id, g) X_pat[g, id], PATIENT_ID, gene)) %>%
  inner_join(
    clin2 %>% select(PATIENT_ID, time, event, age, stage, SUBTYPE),
    by = "PATIENT_ID"
  )

# ---------- Cox formula (use covariates only if available) ----------
avail_covars <- c(
  if (any(!is.na(dat$age)))   "age",
  if (any(!is.na(dat$stage))) "stage",
  "SUBTYPE"
)

cox_formula <- as.formula(
  paste(
    "Surv(time, event) ~ expr",
    if (length(avail_covars)) paste("+", paste(avail_covars, collapse = " + ")) else ""
  )
)

# ---------- Multivariable Cox per gene ----------
fit_one <- function(df, g, form) {
  d <- df %>% filter(gene == g)
  
  needed <- all.vars(form)  # includes expr + covariates names
  needed <- setdiff(needed, c("Surv","time","event"))  # just in case
  d <- d %>% select(time, event, all_of(needed)) %>% na.omit()
  
  if (nrow(d) < 50 || length(unique(d$event)) < 2 || sd(d$expr) == 0) {
    return(tibble(gene = g, n = nrow(d),
                  HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_))
  }
  
  fit <- coxph(form, data = d)
  ti  <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  row <- ti %>% filter(term == "expr")
  
  tibble(
    gene = g,
    n    = nrow(d),
    HR   = row$estimate,
    LCL  = row$conf.low,
    UCL  = row$conf.high,
    p    = row$p.value
  )
}

cox_tbl <- purrr::map_dfr(present, ~fit_one(dat, .x, cox_formula)) %>%
  mutate(gene = factor(gene, levels = present))

# ---------- Forest plot (adjusted HR for gene term only) ----------
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
        "Adj. HR %.2f (%.2f–%.2f)\n p=%s",
        HR, LCL, UCL,
        ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))
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
    x = "Adjusted Hazard Ratio\nlog2(RSEM+1) expression",
    y = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.title.x    = element_text(size = 20),
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    plot.margin     = margin(10, 260, 10, 10)
  )

print(p_forest)

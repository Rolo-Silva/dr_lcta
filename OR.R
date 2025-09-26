# ==== Multinomial RRR tables (all references from one fit) ====================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(purrr); library(broom); library(nnet)
  library(forcats); library(stringr); library(gt)
})

# -- 1) Load data --------------------------------------------------------------
coverage_data    <- readr::read_csv("coverage_data.csv", show_col_types = FALSE)
all_named_models <- readRDS("all_named_models.rds")
demographics <- readr::read_csv("demographics.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    deprivation_index = index_standardized,
    urbanisation_classification = factor(urbanisation_classification,
                                         levels = c("Urbana","Mixta","Rural")),
    zona = factor(zona, levels = c("centro","norte","sur"))
  ) %>%
  dplyr::select(comuna, deprivation_index, urbanisation_classification, zona) %>%
  tidyr::drop_na()

# Make sure stored models point to current coverage_data
for (mn in names(all_named_models)) {
  if (!is.null(all_named_models[[mn]]$call))
    all_named_models[[mn]]$call$data <- coverage_data
}

# id -> comuna from actual modelling data
id_to_comuna <- coverage_data %>%
  dplyr::distinct(id, comuna) %>%
  dplyr::filter(!is.na(id), !is.na(comuna))

# -- 2) Human labels for your 4-class models ----------------------------------
traj_labels <- list(
  "4class_cubic_nre_dgcc_model" = c("Consistently lower","Stable medium","Stable upper medium","Highest"),
  "4class_cubic_nre_drsc_model" = c("Consistently low","Stable medium","Increasing","Highest decreasing")
)

# -- 3) Helpers ----------------------------------------------------------------
clean_term_labels <- function(df) {
  df %>%
    dplyr::mutate(term = dplyr::case_when(
      term == "deprivation_index" ~ "Deprivation Index (continuous)",
      term == "urbanisation_classificationMixta" ~ "Mixta (Ref: Urban)",
      term == "urbanisation_classificationRural" ~ "Rural (Ref: Urban)",
      term == "zonanorte" ~ "Norte (Ref: Centre)",
      term == "zonasur" ~ "Sur (Ref: Centre)",
      term == "(Intercept)" ~ "(Intercept)",
      TRUE ~ term
    ))
}

# Append explicit baseline rows ("Urban (Ref)", "Centre (Ref)") to a wide table
append_baseline_rows <- function(wide, labels, row_order = NULL) {
  ref_rows <- tibble::tibble(term = c("Urban (Ref)", "Centre (Ref)"))
  for (cc in setdiff(names(wide), "term")) ref_rows[[cc]] <- "Ref"
  out <- dplyr::bind_rows(wide, ref_rows)
  
  # default sensible ordering if not provided
  if (is.null(row_order)) {
    row_order <- c(
      "(Intercept)",
      "Deprivation Index (continuous)",
      "Urban (Ref)", "Mixta (Ref: Urban)", "Rural (Ref: Urban)",
      "Centre (Ref)", "Norte (Ref: Centre)", "Sur (Ref: Centre)"
    )
  }
  out %>%
    dplyr::mutate(term = forcats::fct_relevel(term, row_order)) %>%
    dplyr::arrange(term) %>%
    dplyr::mutate(term = as.character(term))
}

# Robustly get outcome levels used by nnet::multinom
get_multinom_levels <- function(fit) {
  levs <- fit$lev
  if (is.null(levs)) {
    mf <- stats::model.frame(fit)
    y  <- model.response(mf)
    levs <- levels(y)
  }
  levs
}

# ---- Core engine: from one fit, build tables for ALL reference classes -------
.make_all_ref_tables <- function(fit, labels, alpha = 0.95, add_baseline_rows = TRUE) {
  stopifnot(inherits(fit, "multinom"))
  levs <- get_multinom_levels(fit)                # class levels in the fit
  V    <- stats::vcov(fit)                        # (K-1)*p x (K-1)*p block covariance
  
  # Coefficient matrix B for non-baseline classes (rows), possibly a vector if K-1 = 1
  B <- coef(fit)
  if (is.vector(B)) {
    cn <- names(B)
    B  <- matrix(B, nrow = 1, dimnames = list(levs[-1], cn))
  }
  cn <- colnames(B)                               # terms
  p  <- length(cn)
  K  <- length(levs)
  
  # Full K x p beta matrix with zeros on the model baseline row
  fullB <- matrix(0, nrow = K, ncol = p, dimnames = list(levs, cn))
  fullB[levs[-1], ] <- B
  
  # helper: index block for a class in V
  block_ix <- function(lvl) {
    if (lvl == levs[1]) return(integer(0))        # model baseline has no block in V
    j <- match(lvl, levs[-1])                     # 1..K-1
    ((j - 1) * p + 1):(j * p)
  }
  
  z <- stats::qnorm(0.5 + alpha/2)
  out <- setNames(vector("list", K), labels)      # name tables by human labels
  
  # iterate desired reference r (each class)
  for (r in seq_len(K)) {
    ref_lvl <- levs[r]
    
    # estimates for k vs r: beta_k - beta_r
    est <- fullB - matrix(fullB[ref_lvl, ], nrow = K, ncol = p, byrow = TRUE)
    
    # standard errors via Var(bk - br) = Var(bk) + Var(br) - 2 Cov(bk, br)
    se <- matrix(NA_real_, nrow = K, ncol = p, dimnames = dimnames(fullB))
    idx_r <- block_ix(ref_lvl)
    for (k in seq_len(K)) {
      idx_k <- block_ix(levs[k])
      if (length(idx_k) && length(idx_r)) {
        Vk  <- V[idx_k, idx_k, drop = FALSE]
        Vr  <- V[idx_r, idx_r, drop = FALSE]
        Ckr <- V[idx_k, idx_r, drop = FALSE]
        se[k, ] <- sqrt(diag(Vk + Vr - 2 * Ckr))
      } else if (length(idx_k)) {
        se[k, ] <- sqrt(diag(V[idx_k, idx_k, drop = FALSE]))
      } else if (length(idx_r)) {
        se[k, ] <- sqrt(diag(V[idx_r, idx_r, drop = FALSE]))
      } else {
        se[k, ] <- 0
      }
    }
    
    # RRR and CI
    OR <- exp(est); LO <- exp(est - z * se); HI <- exp(est + z * se)
    
    # Build a wide table: one row per term, one column per class label
    rows <- lapply(seq_len(p), function(j) {
      vals <- vapply(seq_len(K), function(k) {
        if (k == r) "Ref" else sprintf("%.2f [%.2f, %.2f]", OR[k, j], LO[k, j], HI[k, j])
      }, FUN.VALUE = character(1))
      setNames(c(cn[j], vals), c("term", labels))
    })
    tab <- dplyr::bind_rows(rows) %>%
      clean_term_labels() %>%
      dplyr::select(term, dplyr::all_of(labels))
    
    if (add_baseline_rows) {
      tab <- append_baseline_rows(tab, labels)
    }
    out[[ labels[r] ]] <- tab
  }
  
  out
}

# -- 4) Build merged data & fit once, then all-reference tables ----------------
fit_multinom_tables_multi_ref <- function(model_name, alpha = 0.95, add_baseline_rows = TRUE) {
  model <- all_named_models[[model_name]]
  stopifnot(!is.null(model$call$data), !is.null(model$pred), !is.null(model$pprob))
  
  outcome_var <- if (grepl("drsc", model_name, ignore.case = TRUE)) "drsc" else "dgcc"
  labels      <- traj_labels[[model_name]]
  if (is.null(labels)) stop("No labels defined for ", model_name)
  
  # Merge to get one row per id with covariates + final class
  model_data <- model$call$data %>%
    dplyr::arrange(id, year) %>%
    dplyr::group_by(id) %>% dplyr::mutate(row_id = dplyr::row_number()) %>% dplyr::ungroup()
  
  pred_df <- model$pred %>%
    dplyr::group_by(id) %>% dplyr::mutate(row_id = dplyr::row_number()) %>% dplyr::ungroup()
  
  merged_df <- dplyr::left_join(model_data, pred_df, by = c("id","row_id")) %>%
    dplyr::left_join(model$pprob, by = "id") %>%
    dplyr::filter(!is.na(class))
  
  # Order latent classes low -> high by mean outcome
  class_order <- merged_df %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(mean_val = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_val) %>%
    dplyr::mutate(class_order = dplyr::row_number())
  
  df <- merged_df %>%
    dplyr::left_join(class_order, by = "class") %>%
    dplyr::mutate(class = class_order) %>%
    dplyr::distinct(id, class) %>%
    dplyr::left_join(id_to_comuna, by = "id") %>%
    dplyr::left_join(demographics, by = "comuna") %>%
    dplyr::mutate(
      class_lab = factor(class, levels = seq_along(labels), labels = labels),
      urbanisation_classification = forcats::fct_relevel(urbanisation_classification, "Urbana","Mixta","Rural"),
      zona = forcats::fct_relevel(zona, "centro","norte","sur")
    ) %>%
    tidyr::drop_na(deprivation_index, urbanisation_classification, zona)
  
  # Fit ONCE (baseline arbitrary; we’ll re-express to every ref)
  df$class <- df$class_lab
  df$class <- forcats::fct_relevel(df$class, labels[1])  # make first label the model baseline
  
  fit <- nnet::multinom(class ~ deprivation_index + urbanisation_classification + zona,
                        data = df, trace = FALSE)
  
  .make_all_ref_tables(fit, labels, alpha = alpha, add_baseline_rows = add_baseline_rows)
}

# -- 5) Pretty printers --------------------------------------------------------
as_gt <- function(tab, title_md) {
  tab %>%
    gt::gt() %>%
    gt::tab_options(
      table.font.size = 10,
      data_row.padding = gt::px(0),
      table.border.top.color = "black",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white"
    ) %>%
    gt::tab_header(title = gt::md(title_md))
}

as_gt_all <- function(tabs, title_prefix) {
  purrr::imap(tabs, ~ as_gt(.x, sprintf("**%s (Ref: %s)**", title_prefix, .y)))
}

# -- 6) Run for your two models -----------------------------------------------
# DGCC (4-class)
tabs_dgcc_all <- fit_multinom_tables_multi_ref("4class_cubic_nre_dgcc_model")
# DRSC (4-class; will include "Increasing" among refs)
tabs_drsc_all <- fit_multinom_tables_multi_ref("4class_cubic_nre_drsc_model")

# Convenience extracts to match your previous tables
tab_dgcc_highest     <- tabs_dgcc_all[["Highest"]]
tab_dgcc_sum     <- tabs_dgcc_all[["Stable upper medium"]]
tab_dgcc_sm     <- tabs_dgcc_all[["Stable medium"]]
tab_dgcc_cl     <- tabs_dgcc_all[["Consistently lower"]]

tab_drsc_highest     <- tabs_drsc_all[[ tail(traj_labels[["4class_cubic_nre_drsc_model"]], 1) ]]  # "Highest decreasing"
tab_drsc_increasing  <- tabs_drsc_all[["Increasing"]]
tab_drsc_sm  <- tabs_drsc_all[["Stable medium"]]
tab_drsc_cl  <- tabs_drsc_all[["Consistently low"]]


# Optional: render as gt objects (lists keyed by reference label)
gt_dgcc_all <- as_gt_all(tabs_dgcc_all, "Table S5. Multinomial RRR — DGCC")
gt_drsc_all <- as_gt_all(tabs_drsc_all, "Table S6–S7. Multinomial RRR — DRSC")

# Example: print the specific panels you were using before
gt_dgcc_all[["Highest"]]
gt_drsc_all[["Highest decreasing"]]
gt_drsc_all[["Increasing"]]

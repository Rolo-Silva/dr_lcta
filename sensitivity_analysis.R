# --- packages
library(lcmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(scales)
library(stringr)
library(patchwork)

coverage_data <- read.csv("coverage_data.csv")

# helper: polynomial RHS up to degree 'deg'
.poly_terms <- function(time_var = "year", deg = 3) {
  stopifnot(deg %in% 1:3)
  paste(
    c(time_var,
      if (deg >= 2) sprintf("I(%s^2)", time_var),
      if (deg >= 3) sprintf("I(%s^3)", time_var)),
    collapse = " + "
  )
}

# helper: find time var in data
.find_time_var <- function(df, prefer = c("year","ano","time","t","Time","timevar")) {
  v <- intersect(prefer, names(df))
  if (length(v)) v[1] else stop("Time variable not found in data (tried: ", paste(prefer, collapse=", "), ").")
}

run_sensitivity_analysis <- function(var_name,
                                     coverage_data,
                                     K = ifelse(var_name == "dgcc", 4L, 4L),
                                     time_var = NULL,
                                     prop = 0.8,
                                     seed = 123,
                                     poly_degree = 3) {
  # --- basic checks
  stopifnot(var_name %in% c("dgcc","drsc"))
  if (is.null(time_var)) time_var <- .find_time_var(coverage_data)
  req_cols <- c("id", var_name, time_var)
  miss <- setdiff(req_cols, names(coverage_data))
  if (length(miss)) stop("Missing required columns in coverage_data: ", paste(miss, collapse=", "))
  
  set.seed(seed)
  
  # sample ids
  ids <- unique(coverage_data$id)
  ids_sampled <- sample(ids, size = floor(prop * length(ids)))
  subset_data <- coverage_data %>% filter(id %in% ids_sampled)
  
  # build formulas
  rhs       <- .poly_terms(time_var, poly_degree)             # e.g., "year + I(year^2) + I(year^3)"
  fixed_f   <- as.formula(sprintf("%s ~ %s", var_name, rhs))  # two-sided
  mixture_f <- as.formula(paste("~", rhs))                    # <-- one-sided is REQUIRED
  
  # 1-class starter
  model_1class <- hlme(
    fixed   = fixed_f,
    random  = ~ -1,
    idiag   = FALSE,
    nwg     = FALSE,
    subject = "id",
    ng      = 1,
    data    = subset_data
  )
  
  # K-class fit with gridsearch
  model_K <- tryCatch(
    gridsearch(
      rep     = 20,
      maxiter = 1000,
      minit   = model_1class,
      hlme(
        fixed   = fixed_f,
        mixture = mixture_f,
        random  = ~ -1,
        idiag   = FALSE,
        nwg     = FALSE,
        subject = "id",
        ng      = K,
        data    = subset_data
      )
    ),
    error = function(e) e
  )
  if (inherits(model_K, "error")) {
    warning(sprintf("Gridsearch failed (seed %s, K=%s, var=%s): %s",
                    seed, K, var_name, model_K$message))
    return(list(plot = ggplot() + theme_void(),
                summary = tibble(), model = NULL, error = model_K))
  }
  
  # predictions
  pred <- model_K$pred
  if (is.null(pred)) {
    warning("No predictions returned; plotting blank panel.")
    return(list(plot = ggplot() + theme_void(),
                summary = tibble(), model = model_K))
  }
  
  # class-wise prediction columns
  pred_cols <- grep("^pred_m\\d+$", names(pred), value = TRUE)
  if (!length(pred_cols)) {
    warning("No class-specific prediction columns (pred_m*).")
    return(list(plot = ggplot() + theme_void(),
                summary = tibble(), model = model_K))
  }
  
  # align by id + row index
  id_col <- names(pred)[1]
  data_aligned <- subset_data %>%
    arrange(.data[[id_col]], .data[[time_var]]) %>%
    group_by(.data[[id_col]]) %>%
    mutate(row_id = dplyr::row_number()) %>%
    ungroup()
  
  pred_aligned <- pred %>%
    group_by(.data[[id_col]]) %>%
    mutate(row_id = dplyr::row_number()) %>%
    ungroup()
  
  combined <- left_join(
    data_aligned %>% select(all_of(id_col), row_id, all_of(time_var)),
    pred_aligned %>% select(all_of(id_col), row_id, all_of(pred_cols)),
    by = c(id_col, "row_id")
  )
  
  # long format & summarize
  pred_long <- combined %>%
    pivot_longer(all_of(pred_cols), names_to = "class_label", values_to = "pred") %>%
    mutate(
      class = readr::parse_number(class_label),
      pred  = as.numeric(pred),
      yearv = .data[[time_var]]
    )
  
  pred_summary <- pred_long %>%
    group_by(class, yearv) %>%
    summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")
  
  class_summary <- pred_long %>%
    group_by(class) %>%
    summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      n_units   = n_distinct(.data[[id_col]]),
      .groups   = "drop"
    ) %>%
    mutate(seed = seed, variable = var_name, K = K)
  
  # labels
  y_axis_label <- if (var_name == "dgcc") {
    "Predicted proportion of T2DM individuals with HbA1c < 7%"
  } else {
    "Predicted proportion of T2DM individuals with annual DR screening"
  }
  
  # x-axis labeling
  yrs <- sort(unique(pred_summary$yearv))
  x_breaks <- if (is.numeric(yrs)) pretty(yrs) else yrs
  x_labels <- x_breaks
  x_lab    <- time_var
  
  plot_title <- sprintf("Sensitivity (%.0f%%) – %s – Seed %s (K=%d)",
                        prop * 100,
                        ifelse(var_name == "dgcc", "Glycaemic control", "DR screening"),
                        seed, K)
  
  p <- ggplot(pred_summary, aes(x = yearv, y = pred, color = factor(class))) +
    geom_line(linewidth = 1) +
    ggsci::scale_color_lancet(name = "Class") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = plot_title, x = x_lab, y = y_axis_label) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text  = element_text(size = 8),
      axis.text    = element_text(size = 7),
      axis.title   = element_text(size = 8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  list(plot = p, summary = class_summary, model = model_K)
}

# -------- run: DGCC = 5 classes
sens_dgcc_123 <- run_sensitivity_analysis("dgcc", coverage_data, K = 4, seed = 123)
sens_dgcc_456 <- run_sensitivity_analysis("dgcc", coverage_data, K = 4, seed = 456)
sens_dgcc_789 <- run_sensitivity_analysis("dgcc", coverage_data, K = 4, seed = 789)

# -------- run: DRSC = 4 classes
sens_drsc_123 <- run_sensitivity_analysis("drsc", coverage_data, K = 4, seed = 123)
sens_drsc_456 <- run_sensitivity_analysis("drsc", coverage_data, K = 4, seed = 456)
sens_drsc_789 <- run_sensitivity_analysis("drsc", coverage_data, K = 4, seed = 789)

# quick panels
(dgcc_panel <- (sens_dgcc_123$plot | sens_dgcc_456$plot | sens_dgcc_789$plot) + plot_layout(guides = "collect"))
(drsc_panel <- (sens_drsc_123$plot | sens_drsc_456$plot | sens_drsc_789$plot) + plot_layout(guides = "collect"))

# save
saveRDS(sens_dgcc_123, "sens_dgcc_seed123.rds")
saveRDS(sens_dgcc_456, "sens_dgcc_seed456.rds")
saveRDS(sens_dgcc_789, "sens_dgcc_seed789.rds")
saveRDS(sens_drsc_123, "sens_drsc_seed123.rds")
saveRDS(sens_drsc_456, "sens_drsc_seed456.rds")
saveRDS(sens_drsc_789, "sens_drsc_seed789.rds")

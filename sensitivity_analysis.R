
library(lcmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(patchwork)

coverage_data <- read.csv("coverage_data.csv")

run_sensitivity_analysis <- function(var_name, original_model, original_plot, coverage_data, prop = 0.8, seed = 123) {
  set.seed(seed)
  ids_sampled <- sample(unique(coverage_data$id), size = prop * length(unique(coverage_data$id)))
  subset_data <- coverage_data %>% filter(id %in% ids_sampled)
  
  model_1class <- hlme(
    fixed = as.formula(paste0(var_name, " ~ year + I(year^2) + I(year^3)")),
    random = ~ -1,
    idiag = FALSE,
    nwg = FALSE,
    subject = "id",
    ng = 1,
    data = subset_data
  )
  
  model_5class <- gridsearch(
    rep = 20,
    maxiter = 1000,
    minit = model_1class,
    hlme(
      fixed = as.formula(paste0(var_name, " ~ year + I(year^2) + I(year^3)")),
      mixture = ~ year + I(year^2) + I(year^3),
      random = ~ -1,
      idiag = FALSE,
      nwg = FALSE,
      subject = "id",
      ng = 5,
      data = subset_data
    )
  )
  
  pred <- model_5class$pred
  
  # Añadir row_id para hacer el join
  data <- subset_data %>%
    arrange(id, year) %>%
    group_by(id) %>%
    mutate(row_id = row_number()) %>%
    ungroup()
  
  
  pred <- pred %>%
    group_by(id) %>%
    mutate(row_id = row_number()) %>%
    ungroup()
  
  combined <- left_join(
    data %>% select(id, row_id, year),
    pred %>% select(id, row_id, starts_with("pred_m")),
    by = c("id", "row_id")
  )
  
  cat("🔍 year summary:\n")
  print(summary(combined$year))
  cat("📋 unique years:\n")
  print(unique(combined$year))
  cat("🔬 year typeof:\n")
  print(typeof(combined$year))
  
  
  # Pivotear y limpiar completamente
  pred_long <- combined %>%
    pivot_longer(cols = starts_with("pred_m"), names_to = "class_label", values_to = "pred") %>%
    mutate(class = as.integer(gsub("pred_m", "", class_label)))
  
  # 💥 Reparación definitiva de estructuras para ggplot
  pred_long <- as.data.frame(pred_long)
  pred_long$class <- as.integer(unlist(pred_long$class))
  pred_long$year <- as.integer(unlist(pred_long$year))
  pred_long$pred <- as.numeric(unlist(pred_long$pred))
  
  # Tabla resumen por clase y año
  pred_summary <- pred_long %>%
    group_by(class, year) %>%
    dplyr::summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")
  
  print(str(pred_summary))
  print(sapply(pred_summary, typeof))
  
  
  # Tabla resumen por clase total
  class_summary <- pred_long %>%
    group_by(class) %>%
    dplyr::summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      n_units = n_distinct(id),
      .groups = "drop"
    ) %>%
    mutate(seed = seed, variable = var_name)
  
  # Etiquetas del gráfico
  y_axis_label <- if (var_name == "dgcc") {
    "Predicted proportion of T2DM individuals with HbA1C < 7%"
  } else {
    "Predicted proportion of T2DM individuals with annual DR screening"
  }
  
  plot_title <- paste0("Sensitivity (80%) – ",
                       ifelse(var_name == "dgcc", "Glycaemic control", "Diabetic retinopathy screening"),
                       " – Seed ", seed)
  
  sensitivity_plot <- ggplot(pred_summary, aes(x = year, y = pred, color = factor(class))) +
    geom_line(size = 1) +
    scale_color_lancet(name = "Class") +
    scale_x_continuous(breaks = 0:12, labels = 2011:2023) +
    scale_y_continuous(labels = percent) +
    labs(
      title = plot_title,
      x = "Year",
      y = y_axis_label
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  return(list(plot = sensitivity_plot, summary = class_summary))
}

# Ejecutar sensibilidad para DGCC
sens1 <- run_sensitivity_analysis("dgcc", NULL, NULL, coverage_data, seed = 123)
sens2 <- run_sensitivity_analysis("dgcc", NULL, NULL, coverage_data, seed = 456)
sens3 <- run_sensitivity_analysis("dgcc", NULL, NULL, coverage_data, seed = 789)

# Ejecutar sensibilidad para DRSC
sens4 <- run_sensitivity_analysis("drsc", NULL, NULL, coverage_data, seed = 123)
sens5 <- run_sensitivity_analysis("drsc", NULL, NULL, coverage_data, seed = 456)
sens6 <- run_sensitivity_analysis("drsc", NULL, NULL, coverage_data, seed = 789)


saveRDS(sens1, file = "sens_dgcc_seed123.rds")
saveRDS(sens2, file = "sens_dgcc_seed456.rds")
saveRDS(sens1, file = "sens_dgcc_seed789.rds")

saveRDS(sens4, file = "sens_drsc_seed123.rds")
saveRDS(sens5, file = "sens_drsc_seed456.rds")
saveRDS(sens6, file = "sens_drsc_seed789.rds")

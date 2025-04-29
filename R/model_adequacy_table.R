

# Load required libraries
library(tidyverse) #dplyr, tibble, ggplot2, readr, tidyr, stringr, forcats, lubridate, purr
library(knitr)
library(kableExtra)
library(ggtext)
library(Cairo)
library(extrafont)
library(hrbrthemes)
library(directlabels)
library(ggrepel)
library(readxl)
library(extrafont)
library(scales)
library(ggsci)
library(lcmm)
library(MixAll)
require(kml)
require(traj)
require(lmerTest)
require(plyr)
require(psych)
require(fpc)
require(mclust)
require(rcompanion)
library(gridExtra)
library(tidyLPA)
library(MASS)
library(broom)
library(skimr)
library(gtExtras)
library(pander)
library(BayesFactor)
library(modelsummary)
library(gt)
library(gtsummary)
library(survival)
library(xtable)
library(skimr)
library(htmltools)
library(future)
library(future.apply)
library(LCTMtools)
#devtools::install_github("hlennon/LCTMtools")


# New version

#all_models_by_period <- readRDS("/Users/rolo/Documents/dr_lcmm/all_models_by_period.rds")# Cuando termina de correr el functions script, corro el de model adquacy


# New version
all_models_by_period <- readRDS("all_models_by_period.rds")

# Check the structure of the loaded data
str(all_models_by_period)
# 1. Renombrar modelos con sufijo del periodo
create_named_models_by_period <- function(all_models_by_period, period_suffixes = c("2011_2023", "2011_2019")) {
  names(all_models_by_period) <- period_suffixes
  named_models <- list()
  
  for (i in seq_along(all_models_by_period)) {
    models <- all_models_by_period[[i]]
    suffix <- period_suffixes[i]
    renamed_models <- setNames(models, paste0(names(models), "_", suffix))
    named_models[[i]] <- renamed_models
  }
  
  all_named_models <- do.call(c, named_models)
  return(all_named_models)
}


# Función robusta y corregida para crear summary_table
build_summary_table <- function(model_names, model_list) {
  required_metrics <- c("G", "loglik", "conv", "npm", "AIC", "BIC", 
                        "SABIC", "entropy", "ICL1", "ICL2", "%class")
  
  summary_list <- lapply(model_names, function(model_name) {
    model_obj <- model_list[[model_name]]
    
    if (is.null(model_obj)) {
      warning(paste("Model", model_name, "is NULL — skipped"))
      return(tibble::tibble(Model = model_name))
    }
    
    tryCatch({
      # Extrae todo lo posible
      summ <- summarytable(model_obj, which = required_metrics)
      summ_df <- as.data.frame(summ)
      
      # Tomar la primera fila con métricas principales
      first_row <- summ_df[1, , drop = FALSE]
      
      
      
      # Completar %class1 si no viene
      if (!any(grepl("^%class1$", names(first_row))) && model_obj$ng == 1) {
        first_row$`%class1` <- 100
      }
      
      # Completar %class2 a %class7 si no existen
      for (i in 1:7) {
        col <- paste0("%class", i)
        if (!(col %in% names(first_row))) {
          first_row[[col]] <- NA_real_
        }
      }
      
      # Completar cualquier otro campo ausente
      for (m in c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "ICL1", "ICL2")) {
        if (!(m %in% names(first_row))) {
          first_row[[m]] <- NA_real_
        }
      }
      
      # Agregar nombre del modelo
      first_row$Model <- model_name
      
      # Reordenar
      ordered_cols <- c("Model", "G", "loglik", "conv", "npm", "AIC", "BIC", 
                        "SABIC", "entropy", "ICL1", "ICL2", paste0("%class", 1:7))
      first_row <- first_row[, intersect(ordered_cols, names(first_row))]
      return(first_row)
    }, error = function(e) {
      warning(paste("Error in model", model_name, ":", e$message))
      return(tibble::tibble(Model = model_name))
    })
  })
  
  summary_table <- bind_rows(summary_list)
  return(summary_table)
}


# 3. Métricas de adecuación
extract_postprob <- function(model) {
  tryCatch({
    postprob_values <- postprob(model)
    smallest_class_size_perc <- min(postprob_values[[1]][2, ], na.rm = TRUE)
    smallest_class_count <- min(postprob_values[[1]][1, ], na.rm = TRUE)
    list(smallest_class_size_perc = smallest_class_size_perc, smallest_class_count = smallest_class_count)
  }, error = function(e) {
    list(smallest_class_size_perc = NA, smallest_class_count = NA)
  })
}

extract_occ_appa_mismatch <- function(model) {
  tryCatch({
    toolkit <- LCTMtoolkit(model)
    lower_occ <- if (!is.null(toolkit$occ) && any(!is.na(toolkit$occ))) min(as.numeric(toolkit$occ[1, ]), na.rm = TRUE) else NA
    lower_appa <- if (!is.null(toolkit$appa) && any(!is.na(toolkit$appa))) min(as.numeric(toolkit$appa[1, ]), na.rm = TRUE) else NA
    highest_mismatch <- if (!is.null(toolkit$mismatch) && any(!is.na(toolkit$mismatch))) max(as.numeric(toolkit$mismatch[1, ]), na.rm = TRUE) else NA
    list(lower_occ = lower_occ, lower_appa = lower_appa, highest_mismatch = highest_mismatch)
  }, error = function(e) {
    list(lower_occ = NA, lower_appa = NA, highest_mismatch = NA)
  })
}

extract_vllrt <- function(prev_model, curr_model) {
  tryCatch({
    if (is.null(prev_model) || is.null(curr_model)) return(NA_real_)
    if (!inherits(prev_model, c("hlme", "lcmm")) || !inherits(curr_model, c("hlme", "lcmm"))) return(NA_real_)
    prev_npm <- as.numeric(summarytable(prev_model)[1, "npm"])
    curr_npm <- as.numeric(summarytable(curr_model)[1, "npm"])
    prev_loglik <- as.numeric(prev_model$loglik)
    curr_loglik <- as.numeric(curr_model$loglik)
    LRT_stat <- 2 * (curr_loglik - prev_loglik)
    df_diff <- curr_npm - prev_npm
    p_val <- pchisq(LRT_stat, df = df_diff, lower.tail = FALSE)
    return(p_val)
  }, error = function(e) NA_real_)
}

# 4. Proceso individual por modelo
process_model <- function(model_name, model, prev_model = NULL) {
  tryCatch({
    if (!is.list(model) || is.null(model$ng)) {
      return(tibble(Model = model_name, Error = "Invalid or missing structure"))
    }
    postprob_results <- if (model$ng > 1) extract_postprob(model) else list(smallest_class_size_perc = NA, smallest_class_count = NA)
    occ_appa_mismatch_results <- extract_occ_appa_mismatch(model)
    vllrt_p_value <- if (!is.null(prev_model)) extract_vllrt(prev_model, model) else NA_real_
    
    tibble(
      Model = model_name,
      Smallest_Class_Size_Percentage = postprob_results$smallest_class_size_perc,
      Smallest_Class_Count = postprob_results$smallest_class_count,
      Lowest_OCC = occ_appa_mismatch_results$lower_occ,
      Lowest_APPA = occ_appa_mismatch_results$lower_appa,
      Highest_Mismatch = occ_appa_mismatch_results$highest_mismatch,
      VLMRLRT_P_Value = vllrt_p_value
    )
  }, error = function(e) {
    tibble(Model = model_name, Error = as.character(e))
  })
}

# 5. Procesar todos los modelos
process_all_models <- function(models_list) {
  results <- vector("list", length(models_list))
  prev_model <- NULL
  for (i in seq_along(models_list)) {
    model_name <- names(models_list)[i]
    model <- models_list[[i]]
    results[[i]] <- process_model(model_name, model, prev_model)
    prev_model <- model
  }
  bind_rows(results)
}

# 6. Tabla de adecuación final
create_model_adequacy_table <- function(summary_table, all_named_models) {
  consolidated_summary <- process_all_models(all_named_models)
  summary_table$Model <- as.character(summary_table$Model)
  consolidated_summary$Model <- as.character(consolidated_summary$Model)
  model_adequacy_table <- summary_table %>%
    left_join(consolidated_summary, by = "Model")
  return(model_adequacy_table)
}


# --- EXTRAER COEFICIENTES POR CLASE (ROBUSTO) ---
extract_coeffs_by_class <- function(model) {
  tryCatch({
    if (!inherits(model, c("hlme", "lcmm"))) return(NULL)
    if (model$ng <= 1) return(NULL)
    
    beta <- model$best
    beta_names <- names(beta)
    
    # Detectar coeficientes relacionados con interceptos y pendientes
    intercept_idx <- grep("interc|Intercept|beta0", beta_names, ignore.case = TRUE)
    slope_idx     <- grep("slope|beta1", beta_names, ignore.case = TRUE)
    
    # Si los nombres no ayudan, usar la estructura habitual: primeros ng coef = intercepts, siguientes ng = slopes
    if (length(intercept_idx) < model$ng || length(slope_idx) < model$ng) {
      intercepts <- beta[1:model$ng]
      slopes <- beta[(model$ng + 1):(2 * model$ng)]
    } else {
      intercepts <- beta[intercept_idx][1:model$ng]
      slopes <- beta[slope_idx][1:model$ng]
    }
    
    coeff_df <- data.frame(intercept = as.numeric(intercepts), slope = as.numeric(slopes))
    return(coeff_df)
  }, error = function(e) {
    return(NULL)
  })
}

# --- VERSIÓN CORREGIDA ---
calculate_Mahalanobis_DoS <- function(coeff_df, cov_matrix = NULL, use_identity_if_missing = TRUE) {
  # ⚠️ Validación por si coeff_df es NULL
  if (is.null(coeff_df) || !is.data.frame(coeff_df) || nrow(coeff_df) <= 1) {
    return(NA_real_)
  }
  
  ng <- nrow(coeff_df)
  if (any(is.na(coeff_df$intercept)) || any(is.na(coeff_df$slope))) return(NA_real_)
  
  coef_mat <- as.matrix(coeff_df[, c("intercept", "slope")])
  
  if (is.null(cov_matrix)) {
    if (use_identity_if_missing) {
      cov_matrix <- diag(ncol(coef_mat))
    } else {
      cov_matrix <- cov(coef_mat)
    }
  }
  
  inv_cov <- tryCatch(solve(cov_matrix), error = function(e) NULL)
  if (is.null(inv_cov)) return(NA_real_)
  
  distances <- c()
  for (i in 1:(ng - 1)) {
    for (j in (i + 1):ng) {
      diff_vec <- coef_mat[i, ] - coef_mat[j, ]
      d <- t(diff_vec) %*% inv_cov %*% diff_vec
      distances <- c(distances, sqrt(as.numeric(d)))
    }
  }
  
  return(mean(distances, na.rm = TRUE))
}

# --- AÑADIR MAHALANOBIS DoS A LA TABLA FINAL ---
add_Mahalanobis_DoS_to_model_table <- function(model_list, model_adequacy_table) {
  DoS_df <- lapply(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    coeff_df <- extract_coeffs_by_class(model)
    DoS_value <- calculate_Mahalanobis_DoS(coeff_df)
    data.frame(Model = model_name, DoS_Mahalanobis = round(DoS_value, 4))
  }) %>% do.call(rbind, .)
  
  model_adequacy_table <- model_adequacy_table %>%
    left_join(DoS_df, by = "Model")
  
  return(model_adequacy_table)
}

# --- AÑADIR Y NORMALIZAR ---

# Paso 1: crear modelos renombrados
all_named_models <- create_named_models_by_period(all_models_by_period)



# Paso 2: crear tabla resumen
summary_table <- build_summary_table(names(all_named_models), all_named_models) 

# Paso 3: tabla de adecuación final
model_adequacy_table <- create_model_adequacy_table(summary_table, all_named_models)

model_adequacy_table <- add_Mahalanobis_DoS_to_model_table(all_named_models, model_adequacy_table)

model_adequacy_table <- model_adequacy_table %>% 
  mutate(
    structure = case_when(
      str_detect(Model, "linear_nre_homocedastic") ~ "A",
      str_detect(Model, "linear_nre_heterocedastic") ~ "B",
      str_detect(Model, "quadratic_nre") ~ "C",
      str_detect(Model, "cubic_nre") ~ "D",
      str_detect(Model, "linear_random_intercept_slope") ~ "F",  # ⚠️ must come before "linear_random_intercept"
      str_detect(Model, "linear_random_intercept") ~ "E",
      str_detect(Model, "quadratic_random_effects_prop") ~ "H",  # ⚠️ must come before "quadratic_random_effects"
      str_detect(Model, "quadratic_random_effects") ~ "G",
      str_detect(Model, "cubic_random_effects_prop") ~ "J",      # ⚠️ must come before "cubic_random_effects"
      str_detect(Model, "cubic_random_effects") ~ "I",
      TRUE ~ "Other"
    )) 

write_csv(model_adequacy_table, "model_adequacy_table.csv")

options(scipen = 999)  # very high penalty for scientific notation

View(model_adequacy_table)



# Best Models -------------------------------------------------------------


model_adequacy_table %>% 
  filter(str_detect(Model, "dgcc_model_2011_2023") |
           str_detect(Model, "dgcc_model_2011_2019") |
           str_detect(Model, "drsc_model_2011_2023") |
           str_detect(Model, "drsc_model_2011_2019")) %>%
  
  mutate(
    VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value)),
    coverage = case_when(
      str_detect(Model, "dgcc") ~ "dgcc",
      str_detect(Model, "drsc") ~ "drsc"
    ),
    period = case_when(
      str_detect(Model, "2011_2023") ~ "2011_2023",
      str_detect(Model, "2011_2019") ~ "2011_2019"
    )
  ) %>%
  
  filter(
    Lowest_APPA > 0.70,
    Lowest_OCC > 5,
    entropy > 0.6,
    Smallest_Class_Size_Percentage > 2
  ) %>%
  
  arrange(BIC, -entropy, -Lowest_APPA, -Lowest_OCC, -DoS_Mahalanobis) %>%
  
  group_by(coverage, period) %>%
  #slice(1) %>%
  ungroup() %>% filter(coverage=="dgcc",
                       period=="2011_2023") %>% data.frame() 




# Some testing and filtering in model adeqacy ------------------------------


#Step 1
model_adequacy_table %>% 
  filter(str_detect(Model, pattern = "4class_linear_nre_homocedastic_dgcc_model_2011_2023")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC) %>% 
  slice(1:3)


#Step 3
model_adequacy_table %>% 
  filter(str_detect(Model, "4class") & str_detect(Model, "dgcc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC, -Lowest_APPA, -Lowest_OCC, -entropy, -DoS_Mahalanobis)  %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 2) %>% View()


model_adequacy_table %>% 
  filter(str_detect(Model, pattern = "dgcc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC, -Lowest_APPA, -Lowest_OCC, -entropy, -DoS_Mahalanobis) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 2) 
#slice(1:2) %>% 
#arrange(Model)
model_adequacy_table %>% 
  filter(str_detect(Model, pattern = "dgcc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC, -Lowest_APPA, -Lowest_OCC, -entropy, -DoS_Mahalanobis)  %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 2) %>% View()

model_adequacy_table %>% 
  filter(str_detect(Model, pattern = "class_cubic_random_effects_dgcc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC) 

#Step 4
model_adequacy_table %>% 
  filter(str_detect(Model, "4class") & str_detect(Model, "drsc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC, -Lowest_APPA, -Lowest_OCC, -entropy, -DoS_Mahalanobis) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 2) %>% View()



model_adequacy_table %>% 
  filter(str_detect(Model, "dgcc_model_2011_2023")|
           str_detect(Model, "dgcc_model_2011_2019")|
           str_detect(Model, "drsc_model_2011_2023")|
           str_detect(Model, "drsc_model_2011_2019")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC, -Lowest_APPA, -Lowest_OCC, -entropy, -DoS_Mahalanobis) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 2) %>%
  filter(Model %in% c( "4class_cubic_nre_dgcc_model_2011_2023",
                       "4class_quadratic_nre_dgcc_model_2011_2019",
                       "4class_cubic_nre_drsc_model_2011_2023",
                       "4class_quadratic_random_effects_prop_drsc_model_2011_2019") ) %>% 
  select(Model, G, BIC, entropy, Smallest_Class_Size_Percentage, Lowest_APPA, Highest_Mismatch, DoS_Mahalanobis)



selected_models <- c(
  "4class_cubic_nre_dgcc_model_2011_2023",
  "4class_quadratic_nre_dgcc_model_2011_2019",
  "4class_cubic_nre_drsc_model_2011_2023",
  "4class_quadratic_random_effects_prop_drsc_model_2011_2019"
)

model_adequacy_table %>% 
  filter(Model %in% selected_models) %>% 
  mutate(Model = factor(Model, levels = selected_models)) %>%
  arrange(Model) %>% 
  mutate(Structure= c("cubic_nre",
                      "quadratic_nre", 
                      "cubic_nre",
                      "quadratic_random_effects_prop"),
         Period= c("2011-2023",
                   "2011-2019",
                   "2011-2023",
                   "2011-2019")) %>% 
  select(Period, -Model, Structure, G, BIC, entropy, Smallest_Class_Size_Percentage, Lowest_APPA, Lowest_OCC, Highest_Mismatch, DoS_Mahalanobis) %>% 
  dplyr::rename(
    "Nº of classes" = G, 
    "Relative entropy"= entropy,
    "Smallest class size (%)" = Smallest_Class_Size_Percentage,
    "Lowest APPA" = Lowest_APPA,
    "Lowest OCC" = Lowest_OCC,
    "Highest MMV" = Highest_Mismatch,
    "Mahalanobis distance"= DoS_Mahalanobis) %>% 
  gt() %>%
  tab_options(
    table.font.size = 10,
    data_row.padding = px(1),
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
  tab_row_group(
    group = "Diabetic retinopathy screening coverage",
    rows = 3:4
  ) %>% 
  tab_row_group(
    group = "Diabetic glycemic control coverage",
    rows =1:2
  ) %>% 
  tab_source_note(source_note = md(" **Abbreviations:** LCMM - Latent class mixture model; BIC - Bayesian information criteria; SCS - Smallest class size; APPA: Average class posterior probability; MMV: Mistmatch value; OCC: Odds of correct classification"))



model_adequacy_table %>% 
  filter(Model %in% c( "4class_cubic_nre_dgcc_model_2011_2023",
                       "4class_quadratic_nre_dgcc_model_2011_2019",
                       "4class_cubic_nre_drsc_model_2011_2023",
                       "4class_quadratic_random_effects_prop_drsc_model_2011_2019") ) %>% 
  select(Model, G, BIC, entropy, Smallest_Class_Size_Percentage, Lowest_APPA, Highest_Mismatch, DoS_Mahalanobis) %>% 
  mutate(structure= c("cubic_nre",
                      "cubic_nre", 
                      "quadratic_nre",
                      "quadratic_random_effects_prop"))


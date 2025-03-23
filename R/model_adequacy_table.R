

# Load required libraries
library(dplyr)
library(lcmm)
library(tibble)
library(LCTMtools)
library(stringr)

# Load required libraries
library(dplyr)
library(lcmm)
library(tibble)
library(LCTMtools)
library(stringr)

# 1. Renombrar modelos con sufijo del periodo
create_named_models_by_period <- function(all_models_by_period, period_suffixes = c("2011_2023", "2011_2019", "2020_2023")) {
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


# Paso 1: crear modelos renombrados
all_named_models <- create_named_models_by_period(all_models_by_period)

# Paso 2: crear tabla resumen
summary_table <- build_summary_table(names(all_named_models), all_named_models) 

# Paso 3: tabla de adecuación final
model_adequacy_table <- create_model_adequacy_table(summary_table, all_named_models)

View(model_adequacy_table)




               
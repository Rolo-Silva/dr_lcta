# Install required packages (only run this once)
install.packages(c("lcmm", "dplyr", "future", "future.apply"))

# Load required libraries
library(lcmm)
library(dplyr)
library(future)
library(future.apply)

start_time <- Sys.time()

# Read the CSV file and clean the data
coverage <- read.csv("coverage_2011_2023_noq1.csv") %>%
  arrange(comuna2, ano) %>%
  mutate(
    year = ano - min(ano) + 1,
    id = comuna2,
    drsc = drs_coverage,
    dgcc = dm_coverage
  )

# Define all model structures
model_structures <- list(
  linear_nre_homocedastic = list(fixed = "1 + year", random = "~ -1", nwg = FALSE, idiag = FALSE),
  linear_nre_heterocedastic = list(fixed = "1 + year", random = "~ -1", nwg = TRUE, idiag = FALSE),
  quadratic_nre = list(fixed = "1 + year + I(year^2)", random = "~ -1", nwg = FALSE, idiag = FALSE),
  cubic_nre = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ -1", nwg = FALSE, idiag = FALSE),
  linear_random_intercept = list(fixed = "1 + year", random = "~ 1", nwg = FALSE, idiag = FALSE),
  linear_random_intercept_slope = list(fixed = "1 + year", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  quadratic_random_effects = list(fixed = "1 + year + I(year^2)", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  quadratic_random_effects_prop = list(fixed = "1 + year + I(year^2)", random = "~ 1 + year", nwg = TRUE, idiag = FALSE),
  cubic_random_effects = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  cubic_random_effects_prop = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ 1 + year", nwg = TRUE, idiag = FALSE)
)

# Define function to build models for one period
build_models <- function(data, subject, variable, structure_params, max_classes = 7) {
  results <- list()
  
  
  # Build the 1-class model explicitly
  one_class_model <- tryCatch({
    hlme(
      fixed = as.formula(paste0(variable, " ~ ", structure_params$fixed)),
      random = as.formula(structure_params$random),
      ng = 1,
      nwg = FALSE,
      idiag = structure_params$idiag,
      data = data,
      subject = subject
    )
  }, error = function(e) {
    message("❌ Error building 1class model: ", e$message)
    return(NULL)
  })
  
  if (!is.null(one_class_model)) {
    results[[paste0("1class_", variable, "_model")]] <- one_class_model
    message("✅ Successfully built 1class model")
  }
  
  # Build multi-class models using gridsearch
  for (ng in 2:max_classes) {
    model_name <- paste0(ng, "class_", variable, "_model")
    tryCatch({
      message("Building model: ", model_name)
      multi_class_model <- hlme(
          fixed = as.formula(paste0(variable, " ~ ", structure_params$fixed)),
          mixture = as.formula(paste0(variable, " ~ ", structure_params$fixed)),
          random = as.formula(structure_params$random),
          ng = ng,
          nwg = structure_params$nwg,
          idiag = structure_params$idiag,
          data = data,
          subject = subject
        )
      results[[model_name]] <- multi_class_model
      message("✅ Successfully built: ", model_name)
    }, error = function(e) {
      message("❌ Error in model: ", model_name, "\n", e$message)
    })
  }
  
  return(results)
}

# Function to run models for multiple periods
run_models_for_periods <- function(data, periods, subject, variables, structures, max_classes) {
  all_period_results <- list()
  
  for (period_name in names(periods)) {
    message("Running models for period: ", period_name)
    period_data <- periods[[period_name]]
    
    period_results <- build_all_models(
      data = period_data,
      subject = subject,
      variables = variables,
      structures = structures,
      max_classes = max_classes
    )
    
    all_period_results[[period_name]] <- period_results
  }
  
  return(all_period_results)
}

# Define periods
periods <- list(
  "2011_2023" = coverage,
  "2011_2019" = coverage %>% filter(year >= 1 & year <= 9),
  "2020_2023" = coverage %>% filter(year >= 10 & year <= 13)
)

# Run models for all periods
results_all_periods <- run_models_for_periods(
  data = coverage,
  periods = periods,
  subject = "id",
  variables = c("drsc", "dgcc"),
  structures = model_structures,
  max_classes = 7
)

# Save the results
saveRDS(results_all_periods, file = "results_all_periods.rds")

# End timing
end_time <- Sys.time()

# Calculate execution time
execution_time <- end_time - start_time

# Print the execution time
print(paste("Execution time:", execution_time))












# Load required libraries
library(dplyr)
library(lcmm)
library(targets)
library(tarchetypes)
library(future)

# Load additional required packages
library(tidyverse)
library(lcmm)
library(gridExtra)
library(LCTMtools)
library(tidyLPA)

plan(multisession, workers = availableCores() - 1)

periods <- list(
  "2011_2023" = 1:13,
  "2011_2019" = 1:9,
  "2020_2023" = 10:13
)

# Function to load and preprocess data
load_and_preprocess_data <- function(filepath) {
  read.csv(filepath) %>%
    arrange(comuna2, ano) %>%
    mutate(
      year = ano - min(ano) + 1,
      id = comuna2,
      drsc = drs_coverage,
      dgcc = dm_coverage
    )
}

# Function to subset data for a specific period
subset_data <- function(data, period) {
  data %>% filter(year %in% period)
}

# Function to build models for all periods
build_models_for_all_periods <- function(data, periods, subject, max_classes = 7, variables = c("drsc", "dgcc")) {
  all_models <- list()
  for (period_name in names(periods)) {
    period_data <- subset_data(data, periods[[period_name]])
    all_models[[period_name]] <- build_all_models(
      data = period_data,
      subject = subject,
      max_classes = max_classes,
      variables = variables
    )
  }
  return(all_models)
}

# Function to build all models for a dataset
build_all_models <- function(data, subject, max_classes = 7, variables = c("drsc", "dgcc")) {
  model_structures <- list(
    linear_nre_homocedastic = list(fixed = "1 + year", random = "~ -1", nwg = FALSE, idiag = FALSE),
    linear_nre_heterocedastic = list(fixed = "1 + year", random = "~ -1", nwg = TRUE, idiag = FALSE),
    quadratic_nre = list(fixed = "1 + year + I(year^2)", random = "~ -1", nwg = FALSE, idiag = FALSE),
    cubic_nre = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ -1", nwg = FALSE, idiag = FALSE),
    linear_random_intercept = list(fixed = "1 + year", random = "~ 1", nwg = FALSE, idiag = FALSE),
    linear_random_intercept_slope = list(fixed = "1 + year", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
    quadratic_random_effects = list(fixed = "1 + year + I(year^2)", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
    quadratic_random_effects_prop = list(fixed = "1 + year + I(year^2)", random = "~ 1 + year", nwg = TRUE, idiag = FALSE),
    cubic_random_effects = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
    cubic_random_effects_prop = list(fixed = "1 + year + I(year^2) + I(year^3)", random = "~ 1 + year", nwg = TRUE, idiag = FALSE)
  )
  
  all_models <- list()
  for (structure_name in names(model_structures)) {
    structure_params <- model_structures[[structure_name]]
    for (var in variables) {
      one_class_model <- NULL
      for (ng in 1:max_classes) {
        model_name <- paste0(ng, "class_", structure_name, "_", var, "_model")
        nwg_value <- if (ng == 1) FALSE else structure_params$nwg
        model_args <- list(
          fixed = as.formula(paste0(var, " ~ ", structure_params$fixed)),
          random = as.formula(structure_params$random),
          ng = ng,
          nwg = nwg_value,
          idiag = structure_params$idiag,
          data = data,
          subject = subject
        )
        if (ng > 1) {
          model_args$mixture <- as.formula(paste0(var, " ~ ", structure_params$fixed))
          model_args$B <- one_class_model
        }
        tryCatch({
          fitted_model <- do.call(lcmm::hlme, model_args)
          all_models[[model_name]] <- fitted_model
          if (ng == 1) {
            one_class_model <- fitted_model
          }
        }, error = function(e) {
          message("Error in model fitting for: ", model_name, "\n", e)
          all_models[[model_name]] <- NULL
        })
      }
    }
  }
  return(all_models)
}

coverage_data <- load_and_preprocess_data("coverage_2011_2023_noq1.csv")

all_models_by_period <- build_models_for_all_periods(coverage_data, 
                             periods, 
                             subject = "id",
                             max_classes = 7,
                             variables = c("drsc", "dgcc"))


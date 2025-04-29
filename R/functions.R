# Install required packages (only run this once)
install.packages(c("lcmm", "dplyr", "future", "future.apply"))

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

## All series --------------------------------------------------------------


coverage_2011_2023 <- read.csv("all_series_diabetes.csv", sep=",")

coverage_2011_2023 <- coverage_2011_2023 %>% 
  mutate(zona = ifelse(region %in% c("De Arica y Parinacota", "De Tarapacá", "De Antofagasta", "De Atacama", "De Coquimbo"), 1, 
                       ifelse(region %in% c("De Valparaíso", "Metropolitana de Santiago", "Del Libertador B. O'Higgins", "Del Maule"), 2, 3)))


## coverage ---------------------------------------------------------------

coverage_2011_2023 <- coverage_2011_2023 %>% 
  filter(codigo_prestacion %in% c("P4150602","P4190950","P4190400", "P4180300")) %>% 
  dplyr::group_by(ano, comuna, id_servicio, id_region, zona, codigo_prestacion) %>% 
  dplyr::summarise(cantidad = round(sum(col01))) %>% 
  tidyr::spread(codigo_prestacion, cantidad) %>% 
  dplyr::rename(dm = P4150602,
                dm_fo = P4190950,
                dm_fo_2 =P4190400,
                dm_hg_menor7 = P4180300) %>% 
  dplyr::mutate(dm_fo= coalesce(dm_fo, dm_fo_2)) %>% 
  dplyr::select(-dm_fo_2) %>% 
  dplyr::mutate(drs_coverage = dm_fo/dm,
                dm_coverage = dm_hg_menor7/dm)



# fix duplicated comunas and unify in one only ------------------------------------------------


pac <- coverage_2011_2023 %>% 
  ungroup() %>% 
  filter(str_detect(comuna, "Cerda")) %>% 
  select(-id_region, -id_servicio, -drs_coverage, -dm_coverage) %>% 
  group_by(ano, comuna) %>%
  summarise_all(~ if(is.numeric(.)) sum(., na.rm = TRUE) else first(.)) %>% 
  mutate(drs_coverage = dm_fo/dm,
         dm_coverage = dm_hg_menor7/dm,
         id_servicio = 13,
         id_region =13) %>% 
  select(ano,comuna,id_servicio,id_region,dm,dm_hg_menor7,dm_fo,drs_coverage, dm_coverage )


santiago <- coverage_2011_2023 %>% 
  ungroup() %>% 
  filter(comuna== "Santiago") %>% 
  select(-id_region, -id_servicio, -drs_coverage, -dm_coverage) %>% 
  group_by(ano, comuna) %>%
  summarise_all(~ if(is.numeric(.)) sum(., na.rm = TRUE) else first(.)) %>%
  mutate(drs_coverage = dm_fo/dm,
         dm_coverage = dm_hg_menor7/dm,
         id_servicio = 11,
         id_region =13) %>% 
  select(ano,comuna,id_servicio,id_region,dm,dm_hg_menor7,dm_fo,drs_coverage, dm_coverage )

la_granja <- coverage_2011_2023 %>% 
  ungroup() %>% 
  filter(comuna== "La Granja") %>% 
  select(-id_region, -id_servicio, -drs_coverage, -dm_coverage) %>% 
  group_by(ano, comuna) %>%
  summarise_all(~ if(is.numeric(.)) sum(., na.rm = TRUE) else first(.)) %>% 
  mutate(drs_coverage = dm_fo/dm,
         dm_coverage = dm_hg_menor7/dm,
         id_servicio = 14,
         id_region =13)%>% 
  select(ano,comuna,id_servicio,id_region,dm,dm_hg_menor7,dm_fo,drs_coverage, dm_coverage )


# Remove duplicated comunes -----------------------------------------------


coverage_2011_2023 <- coverage_2011_2023 %>% 
  filter(!comuna %in% c("Santiago", "La Granja"),  
         !str_detect(comuna, "Cerda")) 

#Bind fixed comunes with the main dataset

coverage_2011_2023 <- bind_rows(coverage_2011_2023, pac, santiago, la_granja)



coverage_2011_2023 %>% 
  filter(grepl("La Granja", comuna))



coverage_2011_2023 <- coverage_2011_2023 %>% 
  mutate(zona = ifelse(id_region %in% 13, 2, zona)) 



#Categorise dm quintiles -------------------------------------------------

coverage_2011_2023_noq1 <- coverage_2011_2023 %>%
  ungroup() %>% 
  mutate(quintil_dm_category = cut(dm,
                                   breaks = quantile(dm, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
                                   labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                   include.lowest = TRUE))

  
  
  coverage_2011_2023 %>%
  ungroup() %>% 
  mutate(quintil_dm_category = cut(dm,
                                   breaks = quantile(dm, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
                                   labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                                   include.lowest = TRUE)) %>% 
  dplyr::filter(quintil_dm_category!= "Q1") %>%
  dplyr::mutate(comuna2 = as.integer(factor(comuna, levels = unique(comuna)))) %>% 
  dplyr::mutate(ano2 = ano - 2011) %>% 
  arrange(-drs_coverage) %>% 
  mutate(drs_coverage=replace(drs_coverage, drs_coverage>1, 1)) %>% 
  arrange(ano, comuna) 



## Delete Q1 category ------------------------------------------------------

coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% 
  dplyr::filter(quintil_dm_category!= "Q1")

# Add a simple indentifier to comuna --------------------------------------

coverage_2011_2023_noq1$comuna2 <- match(coverage_2011_2023_noq1$comuna, unique(coverage_2011_2023_noq1$comuna))  ##creating a new comuna variable as a numeric variable (comuna2)

coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% 
  dplyr::mutate(ano2 = ano - 2011)

coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% 
  arrange(-drs_coverage) %>% 
  mutate(drs_coverage=replace(drs_coverage, drs_coverage>1, 1)) %>% 
  arrange(ano, comuna)



isde <- read_excel("SocioEconominoSaludComunas.xlsx")

isde <- isde[,c(2,4)]

isde <- isde %>% 
  dplyr::rename(comuna = ...2)

coverage_isde <- coverage_2011_2023_noq1 %>% 
  mutate(ano2=match(ano, unique(ano)),
         id_region = ifelse(id_region==16, 8, id_region), #Tratar a Ñuble como si hubiera siempre pertencido a una misma region  
         id_region2 = match(id_region, unique(id_region)),
         zona = factor(zona),
         id_servicio2 = match(id_servicio, unique(id_servicio))) %>% 
  dplyr::group_by(comuna, comuna2,id_servicio, id_region, zona) %>% 
  dplyr::summarise(ano=ano,
                  ano2=ano2,
                  dm=dm,
                   dm_hg_menor7=dm_hg_menor7,
                   dm_fo=dm_fo,
                   mean_drs_coverage = mean(drs_coverage, na.rm=T),
                   mean_dm_coverage = mean(dm_coverage, na.rm=T)) %>% 
  distinct(comuna, .keep_all = TRUE) %>% 
  mutate(comuna2 = cur_group_id()) 

coverage_isde[!(coverage_isde$comuna %in% isde$comuna), ] 

isde$comuna[isde$comuna == "Aisen"] <- "Aisén"
isde$comuna[isde$comuna =="Padre las Casas"] <- "Padre Las Casas"
coverage_isde$comuna[174] = isde$comuna[192] ## esto arregla PAC
isde$comuna[isde$comuna =="San Juan de La Costa"] <- "San Juan de la Costa" 


coverage_isde[!(coverage_isde$comuna %in% isde$comuna), ] 




coverage_2011_2023_noq1 <- left_join(coverage_isde, isde) 


coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% 
  ungroup() %>% 
  dplyr::mutate(index_standardized = scale(index)) %>% 
  arrange(comuna)

summary(coverage_2011_2023_noq1)

coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% 
mutate(zona = ifelse(zona == 1, 'norte', 
                     ifelse(zona == 2, 'centro', 'sur')))

rurality <- read_excel("Clasificacion-comunas-PNDR.xlsx")
rurality <- janitor::clean_names(rurality) %>% 
  select(comuna, tipo_com, clasificacion) %>% 
  mutate(urbanisation_level = tipo_com,
         comuna = str_to_title(comuna),
         urbanisation_classification= clasificacion) %>% 
  select(-tipo_com)

rurality$comuna <- gsub("\\bDe\\b", "de", rurality$comuna, ignore.case = TRUE)
rurality$comuna <- gsub("\\bDel\\b", "del", rurality$comuna, ignore.case = TRUE)
rurality$comuna <- gsub("\\b La\\b", " la", rurality$comuna, ignore.case = TRUE)
rurality$comuna[rurality$comuna == "Aysén"] <- "Aisén"
rurality$comuna[rurality$comuna == "Coyhaique"] <- "Coihaique"
rurality$comuna[rurality$comuna == "Alto Biobío"] <- "Alto Bío-Bío"

rurality[!(rurality$comuna %in% coverage_2011_2023_noq1$comuna), ] 

coverage_2011_2023_noq1[!(coverage_2011_2023_noq1$comuna %in% rurality$comuna), ] 

coverage_2011_2023_noq1 <- left_join(coverage_2011_2023_noq1, rurality)

#coverage_2011_2023_noq1 <- coverage_2011_2023_noq1 %>% drop_na(mean_drs_coverage, mean_dm_coverage)


## Save coverage.csv -------------------------------------------------------


write.csv(coverage_2011_2023_noq1, "coverage_2011_2023_noq1.csv")

# Load the .rds file



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





linear_nre_homocedastic 
linear_nre_heterocedastic 
quadratic_nre 
cubic_nre 
linear_random_intercept 
linear_random_intercept_slope 
quadratic_random_effects 
quadratic_random_effects_prop 
cubic_random_effects 
cubic_random_effects_prop


 two outcomes: "drsc", "dgcc"

 two periods: 2011_2023 and 2011_2019
 
 
 
 class_linear_nre_homocedastic_drsc_model_2011_2023
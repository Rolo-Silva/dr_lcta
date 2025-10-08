library(dplyr)
library(stringr)
library(stringi)
# install.packages("fuzzyjoin") if needed
library(fuzzyjoin)
library(readxl)

coverage_data    <- readr::read_csv("coverage_data.csv", show_col_types = FALSE)

idc <- read_excel("indice_desarrollo_comunal_pdf_excel.xlsx")

idc <- idc |> janitor::clean_names()

idc |> skimr::skim()
# 1) Normaliser for comuna names
norm_comuna <- function(x) {
  x |>
    str_trim() |>
    str_squish() |>
    str_to_upper() |>
    # remove accents / diacritics (Á->A, É->E, Ñ->N, etc.)
    stringi::stri_trans_general("Latin-ASCII")
}

# 2) Create canonical keys in both data frames
coverage_clean <- coverage_data|> 
  mutate(comuna_key = norm_comuna(comuna))

idc_clean <- idc|> 
  mutate(comuna_key = norm_comuna(comuna))

# 3) Exact join using the key
joined_exact <- coverage_clean|> 
  left_join(idc_clean|>  select(comuna_key, everything()),
            by = "comuna_key",
            suffix = c("_cov","_idc"))

# 4) Check what didn’t match (if any)
unmatched <- joined_exact|> 
  filter(is.na(bienestar))|>                # replace 'comuna_idc' with a column from idc to test match
  distinct(comuna_cov = comuna_cov, comuna_key)

# You already have: norm_comuna(), coverage_clean, idc_clean
# Build the set of normalized keys that exist in `idc`
idc_keys <- idc_clean|>  distinct(comuna_key)|>  pull()

# Alias map (choose targets that actually exist in `idc`)
alias_map <- tibble::tibble(
  from = c("AISEN", "COIHAIQUE", "CALERA", "NATALES"),
  to   = c(
    if ("AYSEN" %in% idc_keys) "AYSEN" else "AISEN",        # Aisén → Aysén
    if ("COYHAIQUE" %in% idc_keys) "COYHAIQUE" else "COIHAIQUE", # Coihaique → Coyhaique
    if ("LA CALERA" %in% idc_keys) "LA CALERA" else "CALERA",    # Calera (Valpo) → La Calera
    if ("PUERTO NATALES" %in% idc_keys) "PUERTO NATALES" else "NATALES" # Natales → Puerto Natales
  )
)

# Apply aliases to coverage
coverage_fix <- coverage_clean|> 
  left_join(alias_map, by = c("comuna_key" = "from"))|> 
  mutate(comuna_key_fix = coalesce(to, comuna_key))|> 
  select(-to)

# Re-join
coverage_data_idc <- coverage_fix|> 
  left_join(
    idc_clean|>  select(comuna_key, region:rangos, comuna_idc = comuna),
    by = c("comuna_key_fix" = "comuna_key"),
    suffix = c("_cov","_idc")
  )


comunas_pndr <- read_excel("comunas_pndr.xlsx") |> janitor::clean_names()

# --- put this once, after your existing code that builds coverage_data_idc ---
library(dplyr)
library(stringr)
library(stringi)

# normaliser (reuse yours)
norm_comuna <- function(x) x |> str_trim() |> str_squish() |> str_to_upper() |> stringi::stri_trans_general("Latin-ASCII")

# 1) Main join key (use fix if present)
main <- coverage_data_idc %>%
  mutate(join_key = coalesce(comuna_key_fix, comuna_key))

target_keys <- unique(main$join_key)

# 2) Prep PNDR: normalise + only redirect variants TOWARD keys that exist in `main`
pndr_clean <- comunas_pndr %>%
  dplyr::mutate(join_key = norm_comuna(comuna)) %>%
  dplyr::mutate(
    join_key = case_when(
      join_key == "CALERA"         & "LA CALERA"       %in% target_keys ~ "LA CALERA",
      join_key == "NATALES"        & "PUERTO NATALES"  %in% target_keys ~ "PUERTO NATALES",
      join_key == "COIHAIQUE"      & "COYHAIQUE"       %in% target_keys ~ "COYHAIQUE",
      TRUE ~ join_key
    )
  ) %>%
  dplyr::distinct(join_key, .keep_all = TRUE) %>%           # de-dup by key
  dplyr::rename(comuna_pndr = comuna)                        # avoid .x/.y later

# 3) Join: pick the PNDR cols you actually want; rename any that might clash
coverage_data_idc_pndr <- main %>%
  dplyr::left_join(
    pndr_clean %>% select(
      join_key, comuna_pndr,
      cod_reg, region_pndr = region, cod_com, n_habitantes, km_2, densidad, tipo_com, clasificacion
    ),
    by = "join_key"
  )

coverage_data_idc_pndr <- coverage_data_idc_pndr |> 
  dplyr::mutate(zone = ifelse(zona==1, "Northern", zona), #|> View()
         zone = ifelse(zona ==2, "Centre", zone),
         zone = ifelse(zona ==3, "Southern", zone),
         urbanicity = ifelse(clasificacion== "Urbana", "Urban", clasificacion),
         urbanicity = ifelse(clasificacion== "Mixta", "Mixed", urbanicity),
         urbanicity = ifelse(clasificacion== "Rural", "Rural", urbanicity)
         ) 

coverage_data_clean <- coverage_data_idc_pndr

coverage_data_clean |> names()
library(dplyr)

coverage_summary <- coverage_data_clean %>%
  dplyr::filter(comuna != "Melipeuco") %>%
  dplyr::group_by(comuna) %>%
  dplyr::arrange(ano, .by_group = TRUE) %>%   # ensure chronological within comuna
  dplyr::summarise(
    # counts of valid (non-NA) observations
    n_years         = n_distinct(ano),
    n_dm            = sum(!is.na(dm)),
    n_drsc          = sum(!is.na(drsc)),
    n_dgcc          = sum(!is.na(dgcc)),
    n_all3          = sum(!is.na(dm) & !is.na(drsc) & !is.na(dgcc)),
    
    # (optional) proportions of years observed
    prop_dm         = n_dm   / n_years,
    prop_drsc       = n_drsc / n_years,
    prop_dgcc       = n_dgcc / n_years,
    
    # summaries
    mean_dm         = mean(dm,   na.rm = TRUE),
    sd_dm           = sd(dm,     na.rm = TRUE),
    
    median_drsc     = median(drsc, na.rm = TRUE),
    mean_drsc       = mean(drsc,   na.rm = TRUE),
    sd_drsc         = sd(drsc,     na.rm = TRUE),
    
    median_dgcc     = median(dgcc, na.rm = TRUE),
    mean_dgcc       = mean(dgcc,   na.rm = TRUE),
    sd_dgcc         = sd(dgcc,     na.rm = TRUE),
    
    # first/last NON-MISSING (chronological)
    first_drsc      = dplyr::first(na.omit(drsc)),
    last_drsc       = dplyr::last(na.omit(drsc)),
    change_drsc     = last_drsc - first_drsc,
    rel_change_drsc = if_else(is.na(first_drsc) | first_drsc == 0, NA_real_,
                              100 * change_drsc / first_drsc),
    
    first_dgcc      = dplyr::first(na.omit(dgcc)),
    last_dgcc       = dplyr::last(na.omit(dgcc)),
    change_dgcc     = last_dgcc - first_dgcc,
    rel_change_dgcc = if_else(is.na(first_dgcc) | first_dgcc == 0, NA_real_,
                              100 * change_dgcc / first_dgcc),
    .groups = "drop"
  )

additional_data <- coverage_data_clean |> 
  dplyr::select(comuna,
                id,
                id_comuna, 
                id_region, 
                region_cov,
                zone,
                urbanicity,
                bienestar, 
                economia, 
                educacion, 
                idc, 
                ranking,
                n_habitantes, 
                km_2, 
                densidad
                ) |> #, id_region, region_cov, zone, bienestar, economia, educacion, idc, ranking, n_habitantes, km_2, densidad, zone, urbanicity) |> 
  dplyr::distinct() |> 
  dplyr::mutate(scaled_idc = scale(idc)) |> 
  dplyr::filter(comuna!= "Melipeuco")

coverage_summary_demographics <- coverage_summary |> 
  dplyr::left_join(additional_data, by ="comuna")



make_coverage_columns_ms <- function(coverage_data_clean, coverage) {
  cov_name <- rlang::as_name(rlang::ensym(coverage))
  out <- coverage_data_clean %>%
    dplyr::filter(comuna != "Melipeuco") %>%
    dplyr::group_by(comuna, year) %>%
    dplyr::summarise(
      mean = mean({{ coverage }}, na.rm = TRUE)) |> 
    tidyr::pivot_wider(
      names_from  = year,
      values_from = c(mean),
      names_glue  = "{.value}_{year}"
    ) %>%
    dplyr::rename_with(~ paste0(cov_name, "_", .x), -comuna)
  
  out
}

# Example:
coverage_columns_ms_dgcc <- make_coverage_columns_ms(coverage_data_clean, dgcc)
coverage_columns_ms_drsc <- make_coverage_columns_ms(coverage_data_clean, drsc)
coverage_columns_ms_dm <- make_coverage_columns_ms(coverage_data_clean, dm)



coverage_columns <- coverage_columns_ms_dgcc |> 
  dplyr::left_join(coverage_columns_ms_drsc) |> 
  dplyr::left_join(coverage_columns_ms_dm)


coverage_columns_summary_demographics <- coverage_summary_demographics |> 
  dplyr::left_join(coverage_columns, by ="comuna")


all_named_models[["4class_cubic_nre_dgcc_model"]]$pred

dgcc_class <- all_named_models[["4class_cubic_nre_dgcc_model"]]$pprob
dgcc_class <- dgcc_class |> 
              dplyr::rename(class_dgcc = class,
                            prob1_dgcc = prob1,
                            prob2_dgcc = prob2,
                            prob3_dgcc = prob3,
                            prob4_dgcc = prob4)

drsc_class <- all_named_models[["4class_cubic_nre_drsc_model"]]$pprob 
drsc_class <- drsc_class |> 
  dplyr::rename(class_drsc = class,
                prob1_drsc = prob1,
                prob2_drsc = prob2,
                prob3_drsc = prob3,
                prob4_drsc = prob4)


coverage_classes <- dgcc_class |> 
  dplyr::left_join(drsc_class)


coverage_demographics_class <- coverage_columns_summary_demographics |> 
  dplyr::left_join(coverage_classes, by ="id") #ESTE ES OUTPUT TARGETS

  
coverage_demographics_class |> View()


# Summary Table coverage-demographics -------------------------------------



zone_table <- coverage_demographics_class |> 
  dplyr::group_by(zone) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   #mean_dm_population = mean(mean_dm, na.rm=T),
                   #sd_dm = sd(mean_dm, na.rm=T),
                   #T2DM population: median (IQR) across municipalities
                   dm_med   = median(mean_dm, na.rm = TRUE),
                   dm_q1    = quantile(mean_dm, 0.25, na.rm = TRUE),
                   dm_q3    = quantile(mean_dm, 0.75, na.rm = TRUE),
                   #Coverage mean
                   mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                   sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                   mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   sd_drsc = sd(mean_drsc*100, na.rm=T),
                   #CDI mean
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
  gt::gt() |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(dm_med, dm_q1, dm_q3)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
  tab_spanner(label = "DGCC (%)",            columns = c(dm_med, sd_dgcc)) %>%
  #tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency_n      = "N",
    frequency_percent = "%",
    #mean_dm_population        = "Mean",
    dm_med ="Median",
    dm_q1 ="p25",
    dm_q3="p75",
    #sd_dm     = "SD",
    mean_dgcc_class  = "Mean",
    sd_dgcc    = "SD",
    mean_drsc_class  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent, mean_dgcc_class, sd_dgcc
                         ,mean_drsc_class, sd_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(dm_med, dm_q1, dm_q3), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))




urbanicity_table <- coverage_demographics_class |> 
  dplyr::group_by(urbanicity) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   #mean_dm_population = mean(mean_dm, na.rm=T),
                   #sd_dm = sd(mean_dm, na.rm=T),
                   #T2DM population: median (IQR) across municipalities
                   dm_med   = median(mean_dm, na.rm = TRUE),
                   dm_q1    = quantile(mean_dm, 0.25, na.rm = TRUE),
                   dm_q3    = quantile(mean_dm, 0.75, na.rm = TRUE),
                   #Coverage mean
                   mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                   sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                   mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   sd_drsc = sd(mean_drsc*100, na.rm=T),
                   #CDI mean
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
  gt::gt() |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(dm_med, dm_q1, dm_q3)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
  tab_spanner(label = "DGCC (%)",            columns = c(dm_med, sd_dgcc)) %>%
  tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency_n      = "N",
    frequency_percent = "%",
    #mean_dm_population        = "Mean",
    dm_med ="Median",
    dm_q1 ="p25",
    dm_q3="p75",
    #sd_dm     = "SD",
    mean_dgcc_class  = "Mean",
    sd_dgcc    = "SD",
    mean_drsc_class  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent, mean_dgcc_class, sd_dgcc
                         ,mean_drsc_class, sd_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(dm_med, dm_q1, dm_q3), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))


zone_table 
urbanicity_table


# Summary Table Trajectories ----------------------------------------------

dgcc_trajectory_table <- coverage_demographics_class |> 
  dplyr::group_by(class_dgcc) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   #mean_dm_population = mean(mean_dm, na.rm=T),
                   #sd_dm = sd(mean_dm, na.rm=T),
                   #T2DM population: median (IQR) across municipalities
                   dm_med   = median(mean_dm, na.rm = TRUE),
                   dm_q1    = quantile(mean_dm, 0.25, na.rm = TRUE),
                   dm_q3    = quantile(mean_dm, 0.75, na.rm = TRUE),
                   #Coverage mean
                   mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                   sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                   #mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   #sd_drsc = sd(mean_drsc*100, na.rm=T),
                   #CDI mean
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
  gt::gt() |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(dm_med, dm_q1, dm_q3)) %>%
  tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(dm_med, sd_dgcc)) %>%
  #tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency_n      = "N",
    frequency_percent = "%",
    #mean_dm_population        = "Mean",
    dm_med ="Median",
    dm_q1 ="p25",
    dm_q3="p75",
    #sd_dm     = "SD",
    mean_dgcc_class  = "Mean",
    sd_dgcc    = "SD",
    #mean_drsc_class  = "Mean",
    #sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent, mean_dgcc_class, sd_dgcc
                         #,mean_drsc_class, sd_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(dm_med, dm_q1, dm_q3), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))


drsc_trajectory_table <- coverage_demographics_class |> 
  dplyr::group_by(class_drsc) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   #mean_dm_population = mean(mean_dm, na.rm=T),
                   #sd_dm = sd(mean_dm, na.rm=T),
                   #T2DM population: median (IQR) across municipalities
                   dm_med   = median(mean_dm, na.rm = TRUE),
                   dm_q1    = quantile(mean_dm, 0.25, na.rm = TRUE),
                   dm_q3    = quantile(mean_dm, 0.75, na.rm = TRUE),
                   #Coverage mean
                   #mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                   #sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                   mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   sd_drsc = sd(mean_drsc*100, na.rm=T),
                   #CDI mean
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
  gt::gt() |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(dm_med, dm_q1, dm_q3)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(dm_med, sd_dgcc)) %>%
  tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency_n      = "N",
    frequency_percent = "%",
    #mean_dm_population        = "Mean",
    dm_med ="Median",
    dm_q1 ="p25",
    dm_q3="p75",
    #sd_dm     = "SD",
    #mean_dgcc_class  = "Mean",
    #sd_dgcc    = "SD",
    mean_drsc_class  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent #, mean_dgcc_class, sd_dgcc
                         ,mean_drsc_class, sd_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(dm_med, dm_q1, dm_q3), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))

dgcc_trajectory_table


drsc_trajectory_table







make_data_table <- function(demographics_table, group_vars) {
  
  demographics_table <- coverage_demographics_class %>%
    ungroup() |> 
    #group_by(urbanicity) |> 
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      frequency        = dplyr::n(),
      frequency_percent= (frequency /299)*100,
      
      mean_drsc_points = mean(n_drsc),
      mean_dgcc_points = mean(n_dgcc),
      
      # T2DM Population mean
      mean_dm        = mean(mean_dm, na.rm = TRUE),
      sd_mean_dm     = mean(sd_dm,   na.rm = TRUE),
      sd_dm          = sd(mean_dm,   na.rm = TRUE),
      
      
      
      # DGCC (coverage of glycaemic control)
      mean_dgcc        = mean(mean_dgcc*100, na.rm = TRUE),
      sd_mean_dgcc     = mean(sd_dgcc*100,   na.rm = TRUE),  # average within-municipality SD over years
      sd_dgcc          = sd(mean_dgcc*100,   na.rm = TRUE),
      # DRSC (screening coverage)
      mean_drsc        = mean(mean_drsc*100, na.rm = TRUE),
      sd_mean_drsc     = mean(sd_drsc*100,   na.rm = TRUE),
      sd_drsc          = sd(mean_drsc*100,   na.rm = TRUE),
      # between-municipality SD of means
        # average within-municipality SD over years
      
      # IDC
      mean_idc         = mean(idc, na.rm = TRUE),
      sd_idc           = sd(idc,   na.rm = TRUE),
      
      mean_change_dgcc      = mean(change_dgcc, na.rm=TRUE),
      relative_mean_change_dgcc = mean(rel_change_dgcc, na.rm=TRUE),
      
      mean_change_drsc      = mean(change_drsc, na.rm=TRUE),
      relative_mean_change_drsc = mean(rel_change_drsc, na.rm=TRUE),

      .groups = "drop"
    )
  return(demographics_table)
}


zone_table <- make_data_table(demographics_table, "zone")
urbanicity_table <- make_data_table(demographics_table, "urbanicity")
zone_urbanicity_table <- make_data_table(demographics_table, c("zone","urbanicity"))
dgcc_class_trajectory <- make_data_table(demographics_table, "class_dgcc")
drsc_class_trajectory <- make_data_table(demographics_table, "class_drsc")

zone_table_clean <- zone_table |> 
  dplyr::rename(Demographics = zone) 

urbanicity_table_clean <- urbanicity_table |> 
  dplyr::rename(Demographics = urbanicity)


summary_table <- bind_rows(zone_table_clean, urbanicity_table_clean)

summary_table |> names()

gt_tbl <- summary_table |> 
dplyr::select(
  Demographics,
  frequency,
  frequency_percent,
  mean_dm, 
  sd_dm,
  mean_dgcc, 
  sd_dgcc,
  mean_drsc, 
  sd_drsc,
  mean_idc, 
  sd_idc) |> 
gt::gt(rowname_col = "Demographic")  |> 
gt::tab_header(title = "Summary statistics of diabetes coverage and sociodemographics in Chilean municipalities")


coverage_demographics_table <- gt_tbl |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(mean_dm, sd_dm)) %>%
  tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc, sd_dgcc)) %>%
  tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) %>%
  tab_row_group(label= "Geographical Zone", row = 1:3) |> 
  tab_row_group(label= "Urbanicity", row = 4:6) |> 
  # column labels
  cols_label(
    frequency      = "N",
    frequency_percent = "%",
    mean_dm        = "Mean",
    sd_dm     = "SD",
    mean_dgcc  = "Mean",
    sd_dgcc    = "SD",
    mean_drsc  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
fmt_number(columns = c(frequency_percent, mean_dgcc, sd_dgcc,mean_drsc, sd_drsc),
           decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(mean_dm, sd_dm), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))


















dgcc_class_trajectory <- make_data_table(demographics_table, "class_dgcc")
drsc_class_trajectory <- make_data_table(demographics_table, "class_drsc")


# Trajectory table drsc --------------------------------------------------------



gt_tbl_dgcc_trajectories <- dgcc_class_trajectory |>  
  dplyr::select(
    class_dgcc,
    frequency,
    frequency_percent,
    mean_dm, 
    sd_dm,
    mean_dgcc, 
    sd_dgcc,
    mean_idc, 
    sd_idc) |> 
  gt::gt(rowname_col = "Demographic")  |> 
  gt::tab_header(title = "Summary statistics of diabetes coverage and sociodemographics in Chilean municipalities")


coverage_dgcc_trajectory_table <- gt_tbl_dgcc_trajectories |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(mean_dm, sd_dm)) %>%
  tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc, sd_dgcc)) %>%
  #tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc, sd_mean_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency      = "N",
    frequency_percent = "%",
    mean_dm        = "Mean",
    sd_dm     = "SD",
    mean_dgcc  = "Mean",
    sd_dgcc    = "SD",
    #mean_drsc  = "Mean",
    #sd_mean_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent, mean_dgcc, sd_mean_dgcc
                         #,mean_drsc, sd_mean_drsc
                         ),
             decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(mean_dm, sd_dm), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))




gt_tbl_drsc_trajectories <- drsc_class_trajectory |>  
  dplyr::select(
    class_drsc,
    frequency,
    frequency_percent,
    mean_dm, 
    sd_dm,
    mean_drsc, 
    sd_drsc,
    mean_idc, 
    sd_idc) |> 
  gt::gt(rowname_col = "Demographic")  |> 
  gt::tab_header(title = "Summary statistics of diabetes coverage and sociodemographics in Chilean municipalities")


coverage_drsc_trajectory_table <- gt_tbl_drsc_trajectories |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(mean_dm, sd_dm)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc, sd_mean_dgcc)) %>%
  tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc, sd_drsc)) %>%
  tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency      = "N",
    frequency_percent = "%",
    mean_dm        = "Mean",
    sd_dm     = "SD",
    mean_drsc  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent, mean_drsc, sd_mean_drsc
                         #,mean_drsc, sd_mean_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(mean_dm, sd_dm), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))

coverage_dgcc_trajectory_table
coverage_drsc_trajectory_table


coverage_demographics_class |> 
  dplyr::filter(class_dgcc==2) |> 
  select(mean_dm) |> View() |> skimr::skim()


  summarise(mean_dm = mean(mean_dm, na.rm=T),
    sd_dm = sd(mean_dm, na.rm=T))
  
  
  drsc_trajectory_table <- coverage_demographics_class |> 
  dplyr::group_by(class_drsc) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   mean_dm_population = mean(mean_dm, na.rm=T),
                   #median_dm_population = median(mean_dm, na.rm=T),
                   sd_dm = sd(mean_dm, na.rm=T),
                   mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   sd_drsc = sd(mean_drsc *100, na.rm=T),
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
    gt::gt() |> 
    # column groups
    tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
    tab_spanner(label = "T2DM population",     columns = c(mean_dm_population, sd_dm)) %>%
    #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
    tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
    tab_spanner(label = "CDI",   columns = c(mean_idc, sd_idc)) |> 
    # column labels
    cols_label(
      frequency_n      = "N",
      frequency_percent = "%",
      mean_dm_population        = "Mean",
      sd_dm     = "SD",
      mean_drsc_class  = "Mean",
      sd_drsc    = "SD",
      #mean_drsc  = "Mean",
      #sd_mean_drsc    = "SD",
      mean_idc       = "Mean",
      sd_idc         = "SD"
    ) |> 
    fmt_number(columns = c(frequency_percent, mean_drsc_class, sd_drsc
                           #,mean_drsc, sd_mean_drsc
    ),
    decimals = 1) |> 
    fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
    fmt_number(columns = c(mean_dm_population, sd_dm), decimals = 0, use_seps = TRUE) |> 
    cols_align(align = "right", columns = everything()) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
    tab_options(table.font.names = c("Helvetica", "Arial", "sans"))
  
  
  dgcc_trajectory_table <- coverage_demographics_class |> 
    dplyr::group_by(class_dgcc) |> 
    dplyr::summarise(frequency_n = n(),
                     frequency_percent= (frequency_n /299)*100,
                     mean_dm_population = mean(mean_dm, na.rm=T),
                     #median_dm_population = median(mean_dm, na.rm=T),
                     sd_dm = sd(mean_dm, na.rm=T),
                     mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                     sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                     mean_idc = mean(idc, na.rm=T),
                     sd_idc = sd(idc, na.rm=T)) |> 
    gt::gt() |> 
    # column groups
    tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
    tab_spanner(label = "T2DM population",     columns = c(mean_dm_population, sd_dm)) %>%
    tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
    #tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc, sd_mean_drsc)) %>%
    tab_spanner(label = "ICD",   columns = c(mean_idc, sd_idc)) |> 
    # column labels
    cols_label(
      frequency_n      = "N",
      frequency_percent = "%",
      mean_dm_population        = "Mean",
      sd_dm     = "SD",
      mean_dgcc_class  = "Mean",
      sd_dgcc    = "SD",
      #mean_drsc  = "Mean",
      #sd_mean_drsc    = "SD",
      mean_idc       = "Mean",
      sd_idc         = "SD"
    ) |> 
    fmt_number(columns = c(frequency_percent, mean_dgcc_class, sd_dgcc
                           #,mean_drsc, sd_mean_drsc
    ),
    decimals = 1) |> 
    fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
    fmt_number(columns = c(mean_dm_population, sd_dm), decimals = 0, use_seps = TRUE) |> 
    cols_align(align = "right", columns = everything()) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
    tab_options(table.font.names = c("Helvetica", "Arial", "sans"))
 
  drsc_trajectory_table
  dgcc_trajectory_table 

  
multi_trajectory_table_dgcc <- coverage_demographics_class |> 
    dplyr::group_by(class_dgcc, zone, urbanicity) |> 
    dplyr::summarise(frequency_n = n(),
                     frequency_percent= (frequency_n /299)*100,
                     mean_dm_population = mean(mean_dm, na.rm=T),
                     #median_dm_population = median(mean_dm, na.rm=T),
                     sd_dm = sd(mean_dm, na.rm=T),
                     mean_dgcc_class = mean(mean_dgcc*100, na.rm=T),
                     sd_dgcc = sd(mean_dgcc*100, na.rm=T),
                     mean_idc = mean(idc, na.rm=T),
                     sd_idc = sd(idc, na.rm=T)) |> 
    gt::gt() |> 
    # column groups
    tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
    tab_spanner(label = "T2DM population",     columns = c(mean_dm_population, sd_dm)) %>%
    tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
    #tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc, sd_mean_drsc)) %>%
    tab_spanner(label = "ICD",   columns = c(mean_idc, sd_idc)) |> 
    # column labels
    cols_label(
      frequency_n      = "N",
      frequency_percent = "%",
      mean_dm_population        = "Mean",
      sd_dm     = "SD",
      #mean_dgcc_class  = "Mean",
      #sd_dgcc    = "SD",
      mean_drsc_class  = "Mean",
      sd_drsc    = "SD",
      mean_idc       = "Mean",
      sd_idc         = "SD"
    ) |> 
    fmt_number(columns = c(frequency_percent#, mean_dgcc_class, sd_dgcc
                           ,mean_drsc, sd_mean_drsc
    ),
    decimals = 1) |> 
    fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
    fmt_number(columns = c(mean_dm_population, sd_dm), decimals = 0, use_seps = TRUE) |> 
    cols_align(align = "right", columns = everything()) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
    tab_options(table.font.names = c("Helvetica", "Arial", "sans"))
  
multi_trajectory_table_drsc <- coverage_demographics_class |> 
  dplyr::group_by(class_drsc, zone, urbanicity) |> 
  dplyr::summarise(frequency_n = n(),
                   frequency_percent= (frequency_n /299)*100,
                   mean_dm_population = mean(mean_dm, na.rm=T),
                   #median_dm_population = median(mean_dm, na.rm=T),
                   sd_dm = sd(mean_dm, na.rm=T),
                   mean_drsc_class = mean(mean_drsc*100, na.rm=T),
                   sd_drsc = sd(mean_drsc*100, na.rm=T),
                   mean_idc = mean(idc, na.rm=T),
                   sd_idc = sd(idc, na.rm=T)) |> 
  gt::gt() |> 
  # column groups
  tab_spanner(label = "Frequency",           columns = c(frequency_n, frequency_percent)) %>%
  tab_spanner(label = "T2DM population",     columns = c(mean_dm_population, sd_dm)) %>%
  #tab_spanner(label = "DGCC (%)",            columns = c(mean_dgcc_class, sd_dgcc)) %>%
  tab_spanner(label = "DRSC (%)",            columns = c(mean_drsc_class, sd_drsc)) %>%
  tab_spanner(label = "ICD",   columns = c(mean_idc, sd_idc)) |> 
  # column labels
  cols_label(
    frequency_n      = "N",
    frequency_percent = "%",
    mean_dm_population        = "Mean",
    sd_dm     = "SD",
    #mean_dgcc_class  = "Mean",
    #sd_dgcc    = "SD",
    mean_drsc_class  = "Mean",
    sd_drsc    = "SD",
    mean_idc       = "Mean",
    sd_idc         = "SD"
  ) |> 
  fmt_number(columns = c(frequency_percent#, mean_dgcc_class, sd_dgcc
                         ,mean_drsc_class, sd_drsc
  ),
  decimals = 1) |> 
  fmt_number(columns = c(mean_idc, sd_idc), decimals = 2) |> 
  fmt_number(columns = c(mean_dm_population, sd_dm), decimals = 0, use_seps = TRUE) |> 
  cols_align(align = "right", columns = everything()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_title("title")) |> 
  tab_options(table.font.names = c("Helvetica", "Arial", "sans"))

multi_trajectory_table_dgcc
multi_trajectory_table_drsc
  
  

# plots -------------------------------------------------------------------

coverage_data_clean_model <- coverage_data_clean |>
  dplyr::filter(comuna!= "Melipeuco") 
# Load necessary libraries (assuming ggplot2 and scales are used)
# library(dplyr)
# library(ggplot2)
# library(scales) # for percent()

# 1. Data Preparation Function (Revised with Tidy Evaluation)
make_data_coverage <- function(data, outcome){
  # Note: The 'plot_data' argument in the original function is redundant 
  # as it is immediately overwritten. I've renamed 'data' for clarity if you pass 
  # a different initial dataset later, but left the body as is for now.
  plot_data <-coverage_data_clean_model |> 
    dplyr::left_join(coverage_classes, by="id") |>
    dplyr::group_by(year) %>%
    dplyr::summarise(
      mean_coverage = mean({{ outcome }}, na.rm = TRUE),
      sd_coverage = sd({{ outcome }}, na.rm = TRUE),
      n_coverage = n(),
      se_mean = sd_coverage / sqrt(n_coverage),
      margin_error = 1.96 * se_mean,
      lower_ci = mean_coverage - margin_error,
      upper_ci = mean_coverage + margin_error,
      .groups = "drop"
    )
  
  return(plot_data)
}

make_data_coverage(plot_data, drsc)
make_data_coverage(plot_data, dgcc)

# 2. Plotting Function (make_coverage_mean_trends)
# This function takes the summary data and the outcome name for the title.

make_coverage_mean_trends <- function(summary_data, outcome_name){
  
  # Ensure the percent function is available (from the 'scales' package)
  if (!exists("percent")) {
    percent <- scales::percent # Requires scales package
  }
  
  # Create the plot using the provided summary_data
  plot_out <- ggplot(summary_data, aes(x = year, y = mean_coverage)) +
    geom_line(color = "red", linewidth = 0.3) +
    geom_point(size = 1.5, color = "red") +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                  width = 0.25, color = "red") +
    # The target of 0.8 is likely domain-specific (e.g., 80% coverage)
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    scale_x_continuous(
      breaks = seq(0, 12, by = 2), # Start breaks from 0 as years are 0-12
      labels = seq(2010, 2022, by = 2) # Adjust labels to match start year 
    ) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      x = "Year", 
      y = "Average Glycaemic Coverage",
      title = paste("Trend in Average Glycaemic Coverage (", toupper(outcome_name), ")")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold") # Center title
    )
  
  return(plot_out)
}

# 3. Generating the Data and Plots
# Note: Since the plotting function takes the data as its first argument, 
# you can't use the 'plot_data' placeholder anymore for the call.

# 3.1. Generate data for DRSC
drsc_data_summary <- make_data_coverage(data, drsc) # Assuming 'drsc' is a column name

# 3.2. Generate data for DGCC
dgcc_data_summary <- make_data_coverage(data, dgcc) # Assuming 'dgcc' is a column name

# 3.3. Generate and display the DRSC plot
drsc_plot <- make_coverage_mean_trends(
  summary_data = drsc_data_summary, 
  outcome_name = "drsc"
)

# print(drsc_plot) # Uncomment to display the plot

# 3.4. Generate and display the DGCC plot
dgcc_plot <- make_coverage_mean_trends(
  summary_data = dgcc_data_summary, 
  outcome_name = "dgcc"
)



# Plots para trajectories -------------------------------------------------


# PLot trrajectories ------------------------------------------------------

# Original function (renamed for clarity: NO CLASS GROUPING)
make_overall_data_coverage <- function(data, outcome){
  plot_data <- coverage_data_clean_model |> 
    dplyr::left_join(coverage_classes, by="id") |>
    dplyr::group_by(year) %>%
    dplyr::summarise(
      mean_coverage = mean({{ outcome }}, na.rm = TRUE),
      sd_coverage = sd({{ outcome }}, na.rm = TRUE),
      n_coverage = n(),
      se_mean = sd_coverage / sqrt(n_coverage),
      margin_error = 1.96 * se_mean,
      lower_ci = mean_coverage - margin_error,
      upper_ci = mean_coverage + margin_error,
      .groups = "drop"
    )
  return(plot_data)
}

# NEW function (WITH CLASS GROUPING)
# It takes an extra argument 'class_var' for the class column name (e.g., class_drsc)
make_class_data_coverage <- function(data, outcome, class_var){
  plot_data <- coverage_data_clean_model |> 
    dplyr::left_join(coverage_classes, by="id") |>
    # Group by BOTH year and the dynamic class variable
    dplyr::group_by(year, {{ class_var }}) %>%
    dplyr::summarise(
      mean_coverage = mean({{ outcome }}, na.rm = TRUE),
      sd_coverage = sd({{ outcome }}, na.rm = TRUE),
      n_coverage = n(),
      se_mean = sd_coverage / sqrt(n_coverage),
      margin_error = 1.96 * se_mean,
      lower_ci = mean_coverage - margin_error,
      upper_ci = mean_coverage + margin_error,
      .groups = "drop"
    ) %>%
    # Rename the class column to 'class' for consistent plotting later
    dplyr::rename(class = {{ class_var }})
  
  return(plot_data)
}

# NEW Plotting Function for Class-Specific Trajectories
make_class_trajectory_plot <- function(summary_data, outcome_name, class_labels, y_axis_label){
  
  # Set up the class factor for plotting, using the provided labels
  # Assuming the 'class' column in summary_data is 1, 2, 3, 4 after processing
  n_classes <- length(class_labels)
  
  plot_out <- summary_data |>
    dplyr::mutate(class = factor(class, levels = 1:n_classes, labels = class_labels)) |>
    
    ggplot2::ggplot(ggplot2::aes(x = year, y = mean_coverage, colour = class)) +
    
    # Lines, Points, and Error bars (colored by class)
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), width = 0.25) +
    
    # Target line (can be static or dynamic based on class)
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
    
    # Axes and Labels
    ggplot2::scale_x_continuous(
      breaks = seq(0, 12, by = 2), 
      labels = seq(2010, 2022, by = 2) 
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.2), 
      labels = scales::percent
    ) +
    
    ggplot2::xlab("Year") + 
    ggplot2::ylab(y_axis_label) +
    ggplot2::labs(
      title = paste("Class Trajectories for", toupper(outcome_name)),
      color = "Latent Class" # Legend title
    ) +
    
    # Theme and Colors
    ggplot2::theme_bw() +
    # Use a color scale suitable for distinct categories
    ggsci::scale_color_lancet() + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(plot_out)
}

# Pre-defined Labels (as provided in your script)
trajectory_labels <- list(
  "dgcc" = c("Consistently lower", "Stable medium", "Stable upper medium", "Highest"),
  "drsc" = c("Consistently lower", "Stable medium", "Increasing", "Highest")
)

y_axis_labels <- list(
  "dgcc" = "Proportion of T2DM individuals with HbA1c < 7%",
  "drsc" = "Proportion of T2DM individuals with annual DR screening"
)

# --- Example of generating data and plotting for DRSC ---

# 1. Generate Class-Specific Summary Data
drsc_class_summary <- make_class_data_coverage(
  data = coverage_data_clean_model,       # Your main dataset
  outcome = drsc,             # The coverage variable
  class_var = class_drsc      # The class grouping variable
)

View(drsc_class_summary)
# 2. Generate the Class Trajectory Plot
drsc_trajectory_plot <- make_class_trajectory_plot(
  summary_data = drsc_class_summary,
  outcome_name = "drsc",
  class_labels = trajectory_labels$drsc,
  y_axis_label = y_axis_labels$drsc
)

# print(drsc_trajectory_plot) # Uncomment to display the plot

# --- Example for DGCC ---

# 1. Generate Class-Specific Summary Data
dgcc_class_summary <- make_class_data_coverage(
  data = coverage_data_clean_model,
  outcome = dgcc,
  class_var = class_dgcc
)

# 2. Generate the Class Trajectory Plot
dgcc_trajectory_plot <- make_class_trajectory_plot(
  summary_data = dgcc_class_summary,
  outcome_name = "dgcc",
  class_labels = trajectory_labels$dgcc,
  y_axis_label = y_axis_labels$dgcc
)

# print(dgcc_trajectory_plot) # Uncomment to display the plot

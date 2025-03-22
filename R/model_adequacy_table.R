
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




model_adequacy_table <- read_csv("model_adequacy_table.csv")
View(model_adequacy_table)



#Step 1
model_adequacy_table %>% 
  data.frame() %>% 
  filter(str_detect(Model, pattern = "class_linear_nre_homocedastic_drsc_model_2011_2023")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC) %>% data.frame()


#Step 2
model_adequacy_table %>% 
  data.frame() %>% 
  filter(str_detect(Model, pattern = "class_linear_nre_homocedastic_drsc_model_2011_2023")) %>% 
  mutate(VLMRLRT_P_Value = sprintf("%.7f", as.numeric(VLMRLRT_P_Value))) %>% 
  arrange(BIC) %>% 
  head(1)

#Step 3
model_adequacy_table %>% 
  data.frame() %>% 
  filter(str_detect(Model, "3class") & str_detect(Model, "drsc_model_2011_2023")) %>% 
  arrange(BIC)

#Step 4
model_adequacy_table %>% 
  data.frame() %>% 
  filter(str_detect(Model, "2class") & str_detect(Model, "drsc_model_2011_2023")) %>% 
  arrange(BIC) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 5) 





model_adequacy_table %>% 
  data.frame() %>% 
  filter(str_detect(Model, "drsc_model_2011_2023")) %>% 
  arrange(BIC, -Lowest_APPA, -entropy) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         entropy > 0.6,
         Smallest_Class_Size_Percentage > 5) %>% 
  arrange(BIC)





#Summary all steps
model_adequacy_table %>% 
  data.frame() %>% 
  filter(!str_detect(Model, "nre") & str_detect(Model, "drsc_model_2011_2023")) %>% 
  arrange(BIC, -Lowest_APPA, -entropy) %>% 
  filter(Lowest_APPA > 0.70,
         Lowest_OCC> 5,
         Smallest_Class_Size_Percentage > 5) %>% 
  arrange(BIC)












  all_models_by_period[["2011_2023"]][["3class_linear_nre_homocedastic_drsc_model"]]
  
  model_adequacy_table %>% 
    filter(str_detect(Model, "drsc_model_2011_2023") & str_detect(Model, "6class")) %>% 
    arrange(BIC)
  
  model_adequacy_table %>% 
    filter(str_detect(Model, pattern = "cubic_random_effects_drsc_model_2011_2023"))

# Reliability t# Reliability t# Reliability test --------------------------------------------------------


model_adequacy_table[308,]
summary(all_models_by_period[["2011_2023"]][["6class_cubic_random_effects_dgcc_model"]])
LCTMtoolkit(all_models_by_period[["2011_2023"]][["6class_cubic_random_effects_dgcc_model"]])

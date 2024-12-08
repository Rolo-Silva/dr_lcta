
# Eight step framework for Latent class trajectory modelling --------------
library(survival)
library(tidyr)
library(tidyverse)
# library(devtools)
# devtools::install_github("hlennon/LCTMtools")
library(LCTMtools)
library(lcmm)
library(ggplot2)
#install.packages("DT")
library(DT)
summary(oneclass_linear_drsc_model)
summary(oneclass_quadratic_drsc_model)
summary(twoclass_linear_drsc_model)
summary(twoclass_quadratic_drsc_model)
library(dplyr)
library(LCTMtools)
#> 
#> Attaching package: 'LCTMtools'
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     gg_color_hue
#data(bmi_long, package="LCTMtools")



coverage_long_2011_2023 <- coverage_2011_2023_noq1 %>%
  select(comuna2, comuna, ano, drs_coverage, dm_coverage) %>%
  arrange(comuna2, ano) %>% # Sort ascending by year within each comuna2
  mutate(
    year = ano - min(ano) + 1,   # Calculate year as difference from 2011
    id = comuna2,
    drsc = drs_coverage,
    dgcc = dm_coverage
  )            


# Second period: 2011-2019 (actual calendar years 2011 to 2019)
coverage_long_2011_2019 <- coverage_long_2011_2023 %>%
  filter(year %in% 1:9)  # This corresponds to years 1 to 9 (2011-2019)

# Third period: 2020-2023 (actual calendar years 2020 to 2023)
coverage_long_2020_2023 <- coverage_long_2011_2023 %>%
  filter(year %in% 10:13)  # This corresponds to years 10 to 23 (2020-2023)


# Step 1 - models 2011-2023 --------------------------------------------------


## Linear non random effects model (homocedastic) -----------------------------------------


oneclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                              #mixture = ~1+year+I(year^2),
                              random = ~-1,
                              ng = 1,
                              nwg = FALSE, 
                              data=coverage_long_2011_2023,
                              subject = "id")

twoclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                         mixture = drsc~1+year,
                                         random = ~-1,
                                         ng = 2,
                                         nwg = FALSE, 
                                         data=coverage_long_2011_2023,
                                         subject = "id", 
                                         B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)

threeclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 3,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)

fourclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)

fiveclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)

sixclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)

sevenclass_linear_nre_homocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 7,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_homocedastic_drsc_model_2011_2023)



# Save specific models to a file
save(oneclass_linear_nre_homocedastic_drsc_model_2011_2023,
     twoclass_linear_nre_homocedastic_drsc_model_2011_2023,
     threeclass_linear_nre_homocedastic_drsc_model_2011_2023,
     fourclass_linear_nre_homocedastic_drsc_model_2011_2023,
     fiveclass_linear_nre_homocedastic_drsc_model_2011_2023,
     sixclass_linear_nre_homocedastic_drsc_model_2011_2023,
     sevenclass_linear_nre_homocedastic_drsc_model_2011_2023, file = "linear_nre_homocedastic_drsc_model_2011_2023.RData")



## Linear non random heterocedastic -----------------------------------------


oneclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   #mixture = ~1+year+I(year^2),
                                                   random = ~-1,
                                                   ng = 1,
                                                   #nwg = TRUE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id")

twoclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 2,
                                                   nwg = TRUE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)

threeclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 3,
                                                     nwg = TRUE, 
                                                     data=coverage_long_2011_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)

fourclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 4,
                                                    nwg = TRUE, 
                                                    data=coverage_long_2011_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)

fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = TRUE, 
                                                    data=coverage_long_2011_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)

sixclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = TRUE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)

sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 7,
                                                     nwg = TRUE, 
                                                     data=coverage_long_2011_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2023)



save(oneclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     twoclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     threeclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     fourclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     sixclass_linear_nre_heterocedastic_drsc_model_2011_2023,
     sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023, file = "linear_nre_heterocedastic_drsc_model_2011_2023.RData")

## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                            #mixture = ~1+year+I(year^2),
                                            random = ~-1,
                                            ng = 1,
                                            nwg = FALSE, 
                                            data=coverage_long_2011_2023,
                                            subject = "id")

twoclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_nre_drsc_model_2011_2023)


threeclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 3,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_nre_drsc_model_2011_2023)


fourclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 4,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_nre_drsc_model_2011_2023)

fiveclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_nre_drsc_model_2011_2023)

sixclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 6,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_nre_drsc_model_2011_2023)

sevenclass_quadratic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 7,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_nre_drsc_model_2011_2023)


save(oneclass_quadratic_nre_drsc_model_2011_2023,
      twoclass_quadratic_nre_drsc_model_2011_2023,
      threeclass_quadratic_nre_drsc_model_2011_2023,
      fourclass_quadratic_nre_drsc_model_2011_2023,
      fiveclass_quadratic_nre_drsc_model_2011_2023,
      sixclass_quadratic_nre_drsc_model_2011_2023,
      sevenclass_quadratic_nre_drsc_model_2011_2023, file = "quadratic_nre_drsc_model_2011_2023.RData")
## Cubic non random effects model ---------------------------------------


oneclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id")


twoclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_nre_drsc_model_2011_2023)


threeclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 3,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_nre_drsc_model_2011_2023)

fourclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_nre_drsc_model_2011_2023)


fiveclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_nre_drsc_model_2011_2023)

sixclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 6,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_nre_drsc_model_2011_2023)


sevenclass_cubic_nre_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 7,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_nre_drsc_model_2011_2023)

save(oneclass_cubic_nre_drsc_model_2011_2023,
     twoclass_cubic_nre_drsc_model_2011_2023,
     threeclass_cubic_nre_drsc_model_2011_2023,
     fourclass_cubic_nre_drsc_model_2011_2023,
     fiveclass_cubic_nre_drsc_model_2011_2023,
     sixclass_cubic_nre_drsc_model_2011_2023,
     sevenclass_cubic_nre_drsc_model_2011_2023, file = "cubic_nre_drsc_model_2011_2023.RData")

### Random effects intercept ------------------------------------------------


oneclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    #mixture = ~1+year+I(year^2),
                                                                    random = ~1,
                                                                    ng = 1,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2023,
                                                                    subject = "id")

twoclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 2,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2023)

threeclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 3,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2011_2023)

fourclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~1,
                                                                     ng = 4,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2011_2023,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_drsc_model_random_intercept_2011_2023)

fiveclass_linear_drsc_model_random_intercept_2011_2023<- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 5,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2023)

sixclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 6,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2023)

sevenclass_linear_drsc_model_random_intercept_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 7,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2011_2023)

save(oneclass_linear_drsc_model_random_intercept_2011_2023,
     twoclass_linear_drsc_model_random_intercept_2011_2023,
     threeclass_linear_drsc_model_random_intercept_2011_2023,
     fourclass_linear_drsc_model_random_intercept_2011_2023,
     fiveclass_linear_drsc_model_random_intercept_2011_2023,
     sixclass_linear_drsc_model_random_intercept_2011_2023,
     sevenclass_linear_drsc_model_random_intercept_2011_2023, file= "linear_drsc_model_random_intercept_2011_2023.RData")

### Random effects slope ------------------------------------------


oneclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          #mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 1,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2023,
                                                                          subject = "id")

twoclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 2,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2023,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)

threeclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 3,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2011_2023,
                                                                            subject = "id",
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)


fourclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 4,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2011_2023,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)

fiveclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 5,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2011_2023,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)

sixclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 6,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2023,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)

sevenclass_linear_drsc_model_random_intercept_slope_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 7,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2011_2023,
                                                                            subject = "id",  
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2011_2023)


save(oneclass_linear_drsc_model_random_intercept_slope_2011_2023,
     twoclass_linear_drsc_model_random_intercept_slope_2011_2023,
     threeclass_linear_drsc_model_random_intercept_slope_2011_2023,
     fourclass_linear_drsc_model_random_intercept_slope_2011_2023,
     fiveclass_linear_drsc_model_random_intercept_slope_2011_2023,
     sixclass_linear_drsc_model_random_intercept_slope_2011_2023,
     sevenclass_linear_drsc_model_random_intercept_slope_2011_2023, file= "linear_drsc_model_random_intercept_slope_2011_2023.RData")

## Quadratic random effects model Equal ---------------------------------------


oneclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~ 1 + year, 
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data=coverage_long_2011_2023,
                                                      subject = "id")


twoclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 2,                                 
                                                      nwg = FALSE,   
                                                      idiag=FALSE,
                                                      data = coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2023
)


threeclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 3,                                 
                                                        nwg = FALSE,
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2023
)



fourclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 4,                                 
                                                       nwg = FALSE, 
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2023
)


fiveclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 5,                                 
                                                       nwg = FALSE,  
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2023
)

sixclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 6,                                 
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data = coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2023
)

sevenclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 7,                                 
                                                        nwg = FALSE,    
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2023
)


save(oneclass_quadratic_drsc_model_2011_2023,
     twoclass_quadratic_drsc_model_2011_2023,
     threeclass_quadratic_drsc_model_2011_2023,
     fourclass_quadratic_drsc_model_2011_2023,
     fiveclass_quadratic_drsc_model_2011_2023,
     sixclass_quadratic_drsc_model_2011_2023,
     sevenclass_quadratic_drsc_model_2011_2023, file= "quadratic_drsc_model_2011_2023.RData")



## Quadratic random effects model Proportional ---------------------------------------


oneclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           #mixture = ~1+year+I(year^2),
                                                           random = ~ 1 + year, 
                                                           ng = 1,
                                                           #nwg = TRUE, 
                                                           idiag=FALSE,
                                                           data=coverage_long_2011_2023,
                                                           subject = "id")


twoclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 2,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2011_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2011_2023
)



threeclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 3,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2011_2023,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2011_2023
)

fourclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 4,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2011_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2011_2023
)



fiveclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 5,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2011_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2011_2023
)

sixclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 6,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2011_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2011_2023
)



sevenclass_quadratic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 7,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2011_2023,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2011_2023
)


save(oneclass_quadratic_prop_drsc_model_2011_2023,
     twoclass_quadratic_prop_drsc_model_2011_2023,
     threeclass_quadratic_prop_drsc_model_2011_2023,
     fourclass_quadratic_prop_drsc_model_2011_2023,
     fiveclass_quadratic_prop_drsc_model_2011_2023,
     sixclass_quadratic_prop_drsc_model_2011_2023,
     sevenclass_quadratic_prop_drsc_model_2011_2023, file= "quadratic_prop_drsc_model_2011_2023.RData")

## Cubic random effects model Equal ---------------------------------------


oneclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                  #mixture = ~1+year+I(year^2)+I(year^3),
                                                  random = ~ 1 + year, 
                                                  ng = 1,
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data=coverage_long_2011_2023,
                                                  subject = "id")


twoclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 2,                                 
                                                  nwg = FALSE,   
                                                  idiag=FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2011_2023
)


threeclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 3,                                 
                                                    nwg = FALSE,
                                                    idiag=FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2011_2023
)



fourclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 4,                                 
                                                   nwg = FALSE, 
                                                   idiag=FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2011_2023
)


fiveclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 5,                                 
                                                   nwg = FALSE,  
                                                   idiag=FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2011_2023
)

sixclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 6,                                 
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2011_2023
)

sevenclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 7,                                 
                                                    nwg = FALSE,    
                                                    idiag=FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2011_2023
)



save(oneclass_cubic_drsc_model_2011_2023,
     twoclass_cubic_drsc_model_2011_2023,
     threeclass_cubic_drsc_model_2011_2023,
     fourclass_cubic_drsc_model_2011_2023,
     fiveclass_cubic_drsc_model_2011_2023,
     sixclass_cubic_drsc_model_2011_2023,
     sevenclass_cubic_drsc_model_2011_2023, file= "cubic_drsc_model_2011_2023.RData")


## Cubic random effects model Proportional ---------------------------------------


oneclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                       #mixture = ~1+year+I(year^2)+I(year^3),
                                                       random = ~ 1 + year, 
                                                       ng = 1,
                                                       #nwg = TRUE, 
                                                       idiag=FALSE,
                                                       data=coverage_long_2011_2023,
                                                       subject = "id")


twoclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 2,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2011_2023
)



threeclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 3,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2011_2023,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2011_2023
)

fourclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 4,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2011_2023
)



fiveclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 5,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2011_2023
)

sixclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 6,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2011_2023
)



sevenclass_cubic_prop_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 7,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2011_2023,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2011_2023
)

save(oneclass_cubic_prop_drsc_model_2011_2023,
     twoclass_cubic_prop_drsc_model_2011_2023,
     threeclass_cubic_prop_drsc_model_2011_2023,
     fourclass_cubic_prop_drsc_model_2011_2023,
     fiveclass_cubic_prop_drsc_model_2011_2023,
     sixclass_cubic_prop_drsc_model_2011_2023,
     sevenclass_cubic_prop_drsc_model_2011_2023, file= "cubic_prop_drsc_model_2011_2023.RData")

# Step 1 - models 2011-2019 --------------------------------------------------


## Linear non random effects model (homocedastic) -----------------------------------------


oneclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                        #mixture = ~1+year+I(year^2),
                                                                        random = ~-1,
                                                                        ng = 1,
                                                                        nwg = FALSE, 
                                                                        data=coverage_long_2011_2019,
                                                                        subject = "id")

twoclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~-1,
                                                                    ng = 2,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

threeclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 3,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

fourclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~-1,
                                                                     ng = 4,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2011_2019,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

fiveclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~-1,
                                                                     ng = 5,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2011_2019,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

sixclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~-1,
                                                                    ng = 6,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

sevenclass_linear_nre_homocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 7,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_homocedastic_drsc_model_2011_2019)

save(oneclass_linear_nre_homocedastic_drsc_model_2011_2019,
     twoclass_linear_nre_homocedastic_drsc_model_2011_2019,
     threeclass_linear_nre_homocedastic_drsc_model_2011_2019,
     fourclass_linear_nre_homocedastic_drsc_model_2011_2019,
     fiveclass_linear_nre_homocedastic_drsc_model_2011_2019,
     sixclass_linear_nre_homocedastic_drsc_model_2011_2019,
     sevenclass_linear_nre_homocedastic_drsc_model_2011_2019, file = "linear_nre_homocedastic_drsc_model_2011_2019.RData")

## Linear non random heterocedastic -----------------------------------------


oneclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      #mixture = ~1+year+I(year^2),
                                                                      random = ~-1,
                                                                      ng = 1,
                                                                      #nwg = TRUE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id")

twoclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 2,
                                                                      nwg = TRUE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)

threeclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                        mixture = drsc~1+year,
                                                                        random = ~-1,
                                                                        ng = 3,
                                                                        nwg = TRUE, 
                                                                        data=coverage_long_2011_2019,
                                                                        subject = "id", 
                                                                        B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)

fourclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                       mixture = drsc~1+year,
                                                                       random = ~-1,
                                                                       ng = 4,
                                                                       nwg = TRUE, 
                                                                       data=coverage_long_2011_2019,
                                                                       subject = "id", 
                                                                       B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)

fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                       mixture = drsc~1+year,
                                                                       random = ~-1,
                                                                       ng = 5,
                                                                       nwg = TRUE, 
                                                                       data=coverage_long_2011_2019,
                                                                       subject = "id", 
                                                                       B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)

sixclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 6,
                                                                      nwg = TRUE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)

sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                        mixture = drsc~1+year,
                                                                        random = ~-1,
                                                                        ng = 7,
                                                                        nwg = TRUE, 
                                                                        data=coverage_long_2011_2019,
                                                                        subject = "id", 
                                                                        B= oneclass_linear_nre_heterocedastic_drsc_model_2011_2019)


save(oneclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     twoclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     threeclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     fourclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     sixclass_linear_nre_heterocedastic_drsc_model_2011_2019,
     sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019, file = "linear_nre_heterocedastic_drsc_model_2011_2019.RData")

## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          #mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 1,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2011_2019,
                                                          subject = "id")

twoclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 2,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2011_2019,
                                                          subject = "id",
                                                          B=oneclass_quadratic_nre_drsc_model_2011_2019)


threeclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                            mixture = ~1+year+I(year^2),
                                                            random = ~-1,
                                                            ng = 3,
                                                            nwg = FALSE, 
                                                            data=coverage_long_2011_2019,
                                                            subject = "id",
                                                            B=oneclass_quadratic_nre_drsc_model_2011_2019)


fourclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           mixture = ~1+year+I(year^2),
                                                           random = ~-1,
                                                           ng = 4,
                                                           nwg = FALSE, 
                                                           data=coverage_long_2011_2019,
                                                           subject = "id",
                                                           B=oneclass_quadratic_nre_drsc_model_2011_2019)

fiveclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           mixture = ~1+year+I(year^2),
                                                           random = ~-1,
                                                           ng = 5,
                                                           nwg = FALSE, 
                                                           data=coverage_long_2011_2019,
                                                           subject = "id",
                                                           B=oneclass_quadratic_nre_drsc_model_2011_2019)

sixclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 6,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2011_2019,
                                                          subject = "id",
                                                          B=oneclass_quadratic_nre_drsc_model_2011_2019)

sevenclass_quadratic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                            mixture = ~1+year+I(year^2),
                                                            random = ~-1,
                                                            ng = 7,
                                                            nwg = FALSE, 
                                                            data=coverage_long_2011_2019,
                                                            subject = "id",
                                                            B=oneclass_quadratic_nre_drsc_model_2011_2019)

save(oneclass_quadratic_nre_drsc_model_2011_2019,
     twoclass_quadratic_nre_drsc_model_2011_2019,
     threeclass_quadratic_nre_drsc_model_2011_2019,
     fourclass_quadratic_nre_drsc_model_2011_2019,
     fiveclass_quadratic_nre_drsc_model_2011_2019,
     sixclass_quadratic_nre_drsc_model_2011_2019,
     sevenclass_quadratic_nre_drsc_model_2011_2019, file = "quadratic_nre_drsc_model_2011_2019.RData")
## Cubic non random effects model ---------------------------------------


oneclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE,
                                                      data = coverage_long_2011_2019,
                                                      subject = "id")


twoclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE,
                                                      data = coverage_long_2011_2019,
                                                      subject = "id",
                                                      B= oneclass_cubic_nre_drsc_model_2011_2019)


threeclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B= oneclass_cubic_nre_drsc_model_2011_2019)

fourclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B= oneclass_cubic_nre_drsc_model_2011_2019)


fiveclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B= oneclass_cubic_nre_drsc_model_2011_2019)

sixclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE,
                                                      data = coverage_long_2011_2019,
                                                      subject = "id",
                                                      B= oneclass_cubic_nre_drsc_model_2011_2019)


sevenclass_cubic_nre_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        random = ~-1,
                                                        ng = 7,
                                                        nwg = FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B= oneclass_cubic_nre_drsc_model_2011_2019)


save(oneclass_cubic_nre_drsc_model_2011_2019,
     twoclass_cubic_nre_drsc_model_2011_2019,
     threeclass_cubic_nre_drsc_model_2011_2019,
     fourclass_cubic_nre_drsc_model_2011_2019,
     fiveclass_cubic_nre_drsc_model_2011_2019,
     sixclass_cubic_nre_drsc_model_2011_2019,
     sevenclass_cubic_nre_drsc_model_2011_2019, file = "cubic_nre_drsc_model_2011_2019.RData")
### Random effects intercept ------------------------------------------------


oneclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    #mixture = ~1+year+I(year^2),
                                                                    random = ~1,
                                                                    ng = 1,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id")

twoclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 2,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2019)

threeclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 3,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2011_2019)

fourclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~1,
                                                                     ng = 4,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2011_2019,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_drsc_model_random_intercept_2011_2019)

fiveclass_linear_drsc_model_random_intercept_2011_2019<- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 5,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2019)

sixclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 6,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2011_2019,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2011_2019)

sevenclass_linear_drsc_model_random_intercept_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 7,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2011_2019,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2011_2019)

save(oneclass_linear_drsc_model_random_intercept_2011_2019,
     twoclass_linear_drsc_model_random_intercept_2011_2019,
     threeclass_linear_drsc_model_random_intercept_2011_2019,
     fourclass_linear_drsc_model_random_intercept_2011_2019,
     fiveclass_linear_drsc_model_random_intercept_2011_2019,
     sixclass_linear_drsc_model_random_intercept_2011_2019,
     sevenclass_linear_drsc_model_random_intercept_2011_2019, file= "linear_drsc_model_random_intercept_2011_2019.RData")
### Random effects slope ------------------------------------------


oneclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          #mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 1,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2019,
                                                                          subject = "id")

twoclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 2,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2019,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)

threeclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 3,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2011_2019,
                                                                            subject = "id",
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)


fourclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 4,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2011_2019,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)

fiveclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 5,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2011_2019,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)

sixclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 6,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2011_2019,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)

sevenclass_linear_drsc_model_random_intercept_slope_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 7,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2011_2019,
                                                                            subject = "id",  
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2011_2019)


save(oneclass_linear_drsc_model_random_intercept_slope_2011_2019,
     twoclass_linear_drsc_model_random_intercept_slope_2011_2019,
     threeclass_linear_drsc_model_random_intercept_slope_2011_2019,
     fourclass_linear_drsc_model_random_intercept_slope_2011_2019,
     fiveclass_linear_drsc_model_random_intercept_slope_2011_2019,
     sixclass_linear_drsc_model_random_intercept_slope_2011_2019,
     sevenclass_linear_drsc_model_random_intercept_slope_2011_2019, file= "linear_drsc_model_random_intercept_slope_2011_2019.RData")

## Quadratic random effects model Equal ---------------------------------------


oneclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~ 1 + year, 
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data=coverage_long_2011_2019,
                                                      subject = "id")


twoclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 2,                                 
                                                      nwg = FALSE,   
                                                      idiag=FALSE,
                                                      data = coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2019
)


threeclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 3,                                 
                                                        nwg = FALSE,
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2019
)



fourclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 4,                                 
                                                       nwg = FALSE, 
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2019
)


fiveclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 5,                                 
                                                       nwg = FALSE,  
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2019
)

sixclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 6,                                 
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data = coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2019
)

sevenclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 7,                                 
                                                        nwg = FALSE,    
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2019
)


save(oneclass_quadratic_drsc_model_2011_2019,
     twoclass_quadratic_drsc_model_2011_2019,
     threeclass_quadratic_drsc_model_2011_2019,
     fourclass_quadratic_drsc_model_2011_2019,
     fiveclass_quadratic_drsc_model_2011_2019,
     sixclass_quadratic_drsc_model_2011_2019,
     sevenclass_quadratic_drsc_model_2011_2019, file= "quadratic_drsc_model_2011_2019.RData")



## Quadratic random effects model Proportional ---------------------------------------


oneclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           #mixture = ~1+year+I(year^2),
                                                           random = ~ 1 + year, 
                                                           ng = 1,
                                                           #nwg = TRUE, 
                                                           idiag=FALSE,
                                                           data=coverage_long_2011_2019,
                                                           subject = "id")


twoclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 2,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2011_2019,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2011_2019
)



threeclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 3,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2011_2019,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2011_2019
)

fourclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 4,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2011_2019,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2011_2019
)



fiveclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 5,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2011_2019,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2011_2019
)

sixclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 6,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2011_2019,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2011_2019
)



sevenclass_quadratic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 7,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2011_2019,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2011_2019
)


save(oneclass_quadratic_prop_drsc_model_2011_2019,
     twoclass_quadratic_prop_drsc_model_2011_2019,
     threeclass_quadratic_prop_drsc_model_2011_2019,
     fourclass_quadratic_prop_drsc_model_2011_2019,
     fiveclass_quadratic_prop_drsc_model_2011_2019,
     sixclass_quadratic_prop_drsc_model_2011_2019,
     sevenclass_quadratic_prop_drsc_model_2011_2019, file= "quadratic_prop_drsc_model_2011_2019.RData")

## Cubic random effects model Equal ---------------------------------------


oneclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                  #mixture = ~1+year+I(year^2)+I(year^3),
                                                  random = ~ 1 + year, 
                                                  ng = 1,
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data=coverage_long_2011_2019,
                                                  subject = "id")


twoclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 2,                                 
                                                  nwg = FALSE,   
                                                  idiag=FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2011_2019
)


threeclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 3,                                 
                                                    nwg = FALSE,
                                                    idiag=FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2011_2019
)



fourclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 4,                                 
                                                   nwg = FALSE, 
                                                   idiag=FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2011_2019
)


fiveclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 5,                                 
                                                   nwg = FALSE,  
                                                   idiag=FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2011_2019
)

sixclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 6,                                 
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2011_2019
)

sevenclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 7,                                 
                                                    nwg = FALSE,    
                                                    idiag=FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2011_2019
)


save(oneclass_cubic_drsc_model_2011_2019,
     twoclass_cubic_drsc_model_2011_2019,
     threeclass_cubic_drsc_model_2011_2019,
     fourclass_cubic_drsc_model_2011_2019,
     fiveclass_cubic_drsc_model_2011_2019,
     sixclass_cubic_drsc_model_2011_2019,
     sevenclass_cubic_drsc_model_2011_2019, file= "cubic_drsc_model_2011_2019.RData")



## Cubic random effects model Proportional ---------------------------------------


oneclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                       #mixture = ~1+year+I(year^2)+I(year^3),
                                                       random = ~ 1 + year, 
                                                       ng = 1,
                                                       #nwg = TRUE, 
                                                       idiag=FALSE,
                                                       data=coverage_long_2011_2019,
                                                       subject = "id")


twoclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 2,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2011_2019
)



threeclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 3,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2011_2019,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2011_2019
)

fourclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 4,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2011_2019
)



fiveclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 5,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2011_2019
)

sixclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 6,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2011_2019
)



sevenclass_cubic_prop_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 7,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2011_2019,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2011_2019
)


save(oneclass_cubic_prop_drsc_model_2011_2019,
     twoclass_cubic_prop_drsc_model_2011_2019,
     threeclass_cubic_prop_drsc_model_2011_2019,
     fourclass_cubic_prop_drsc_model_2011_2019,
     fiveclass_cubic_prop_drsc_model_2011_2019,
     sixclass_cubic_prop_drsc_model_2011_2019,
     sevenclass_cubic_prop_drsc_model_2011_2019, file= "cubic_prop_drsc_model_2011_2019.RData")



# Step 1 - models 2020-2023 --------------------------------------------------


## Linear non random effects model (homocedastic) -----------------------------------------


oneclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    #mixture = ~1+year+I(year^2),
                                                                    random = ~-1,
                                                                    ng = 1,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id")

twoclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~-1,
                                                                    ng = 2,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)

threeclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 3,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)

fourclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~-1,
                                                                     ng = 4,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2020_2023,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)

fiveclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~-1,
                                                                     ng = 5,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2020_2023,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)

sixclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~-1,
                                                                    ng = 6,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)

sevenclass_linear_nre_homocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 7,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_homocedastic_drsc_model_2020_2023)



# Save specific models to a file
save(oneclass_linear_nre_homocedastic_drsc_model_2020_2023,
     twoclass_linear_nre_homocedastic_drsc_model_2020_2023,
     threeclass_linear_nre_homocedastic_drsc_model_2020_2023,
     fourclass_linear_nre_homocedastic_drsc_model_2020_2023,
     fiveclass_linear_nre_homocedastic_drsc_model_2020_2023,
     sixclass_linear_nre_homocedastic_drsc_model_2020_2023,
     sevenclass_linear_nre_homocedastic_drsc_model_2020_2023, file = "linear_nre_homocedastic_drsc_model_2020_2023.RData")



## Linear non random heterocedastic -----------------------------------------


oneclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      #mixture = ~1+year+I(year^2),
                                                                      random = ~-1,
                                                                      ng = 1,
                                                                      #nwg = TRUE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id")

twoclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 2,
                                                                      nwg = TRUE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)

threeclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                        mixture = drsc~1+year,
                                                                        random = ~-1,
                                                                        ng = 3,
                                                                        nwg = TRUE, 
                                                                        data=coverage_long_2020_2023,
                                                                        subject = "id", 
                                                                        B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)

fourclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                       mixture = drsc~1+year,
                                                                       random = ~-1,
                                                                       ng = 4,
                                                                       nwg = TRUE, 
                                                                       data=coverage_long_2020_2023,
                                                                       subject = "id", 
                                                                       B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)

fiveclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                       mixture = drsc~1+year,
                                                                       random = ~-1,
                                                                       ng = 5,
                                                                       nwg = TRUE, 
                                                                       data=coverage_long_2020_2023,
                                                                       subject = "id", 
                                                                       B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)

sixclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~-1,
                                                                      ng = 6,
                                                                      nwg = TRUE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)

sevenclass_linear_nre_heterocedastic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                        mixture = drsc~1+year,
                                                                        random = ~-1,
                                                                        ng = 7,
                                                                        nwg = TRUE, 
                                                                        data=coverage_long_2020_2023,
                                                                        subject = "id", 
                                                                        B= oneclass_linear_nre_heterocedastic_drsc_model_2020_2023)



save(oneclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     twoclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     threeclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     fourclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     fiveclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     sixclass_linear_nre_heterocedastic_drsc_model_2020_2023,
     sevenclass_linear_nre_heterocedastic_drsc_model_2020_2023, file = "linear_nre_heterocedastic_drsc_model_2020_2023.RData")

## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          #mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 1,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2020_2023,
                                                          subject = "id")

twoclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 2,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2020_2023,
                                                          subject = "id",
                                                          B=oneclass_quadratic_nre_drsc_model_2020_2023)


threeclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                            mixture = ~1+year+I(year^2),
                                                            random = ~-1,
                                                            ng = 3,
                                                            nwg = FALSE, 
                                                            data=coverage_long_2020_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_nre_drsc_model_2020_2023)


fourclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           mixture = ~1+year+I(year^2),
                                                           random = ~-1,
                                                           ng = 4,
                                                           nwg = FALSE, 
                                                           data=coverage_long_2020_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_nre_drsc_model_2020_2023)

fiveclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           mixture = ~1+year+I(year^2),
                                                           random = ~-1,
                                                           ng = 5,
                                                           nwg = FALSE, 
                                                           data=coverage_long_2020_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_nre_drsc_model_2020_2023)

sixclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                          mixture = ~1+year+I(year^2),
                                                          random = ~-1,
                                                          ng = 6,
                                                          nwg = FALSE, 
                                                          data=coverage_long_2020_2023,
                                                          subject = "id",
                                                          B=oneclass_quadratic_nre_drsc_model_2020_2023)

sevenclass_quadratic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                            mixture = ~1+year+I(year^2),
                                                            random = ~-1,
                                                            ng = 7,
                                                            nwg = FALSE, 
                                                            data=coverage_long_2020_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_nre_drsc_model_2020_2023)


save(oneclass_quadratic_nre_drsc_model_2020_2023,
     twoclass_quadratic_nre_drsc_model_2020_2023,
     threeclass_quadratic_nre_drsc_model_2020_2023,
     fourclass_quadratic_nre_drsc_model_2020_2023,
     fiveclass_quadratic_nre_drsc_model_2020_2023,
     sixclass_quadratic_nre_drsc_model_2020_2023,
     sevenclass_quadratic_nre_drsc_model_2020_2023, file = "quadratic_nre_drsc_model_2020_2023.RData")
## Cubic non random effects model ---------------------------------------


oneclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE,
                                                      data = coverage_long_2020_2023,
                                                      subject = "id")


twoclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE,
                                                      data = coverage_long_2020_2023,
                                                      subject = "id",
                                                      B= oneclass_cubic_nre_drsc_model_2020_2023)


threeclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B= oneclass_cubic_nre_drsc_model_2020_2023)

fourclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B= oneclass_cubic_nre_drsc_model_2020_2023)


fiveclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B= oneclass_cubic_nre_drsc_model_2020_2023)

sixclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE,
                                                      data = coverage_long_2020_2023,
                                                      subject = "id",
                                                      B= oneclass_cubic_nre_drsc_model_2020_2023)


sevenclass_cubic_nre_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                        random = ~-1,
                                                        ng = 7,
                                                        nwg = FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B= oneclass_cubic_nre_drsc_model_2020_2023)

save(oneclass_cubic_nre_drsc_model_2020_2023,
     twoclass_cubic_nre_drsc_model_2020_2023,
     threeclass_cubic_nre_drsc_model_2020_2023,
     fourclass_cubic_nre_drsc_model_2020_2023,
     fiveclass_cubic_nre_drsc_model_2020_2023,
     sixclass_cubic_nre_drsc_model_2020_2023,
     sevenclass_cubic_nre_drsc_model_2020_2023, file = "cubic_nre_drsc_model_2020_2023.RData")

### Random effects intercept ------------------------------------------------


oneclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    #mixture = ~1+year+I(year^2),
                                                                    random = ~1,
                                                                    ng = 1,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id")

twoclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 2,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2020_2023)

threeclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 3,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2020_2023)

fourclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                     mixture = drsc~1+year,
                                                                     random = ~1,
                                                                     ng = 4,
                                                                     nwg = FALSE, 
                                                                     data=coverage_long_2020_2023,
                                                                     subject = "id", 
                                                                     B= oneclass_linear_drsc_model_random_intercept_2020_2023)

fiveclass_linear_drsc_model_random_intercept_2020_2023<- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 5,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2020_2023)

sixclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                    mixture = drsc~1+year,
                                                                    random = ~1,
                                                                    ng = 6,
                                                                    nwg = FALSE, 
                                                                    data=coverage_long_2020_2023,
                                                                    subject = "id", 
                                                                    B= oneclass_linear_drsc_model_random_intercept_2020_2023)

sevenclass_linear_drsc_model_random_intercept_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                      mixture = drsc~1+year,
                                                                      random = ~1,
                                                                      ng = 7,
                                                                      nwg = FALSE, 
                                                                      data=coverage_long_2020_2023,
                                                                      subject = "id", 
                                                                      B= oneclass_linear_drsc_model_random_intercept_2020_2023)

save(oneclass_linear_drsc_model_random_intercept_2020_2023,
     twoclass_linear_drsc_model_random_intercept_2020_2023,
     threeclass_linear_drsc_model_random_intercept_2020_2023,
     fourclass_linear_drsc_model_random_intercept_2020_2023,
     fiveclass_linear_drsc_model_random_intercept_2020_2023,
     sixclass_linear_drsc_model_random_intercept_2020_2023,
     sevenclass_linear_drsc_model_random_intercept_2020_2023, file= "linear_drsc_model_random_intercept_2020_2023.RData")

### Random effects slope ------------------------------------------


oneclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          #mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 1,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2020_2023,
                                                                          subject = "id")

twoclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 2,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2020_2023,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)

threeclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 3,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2020_2023,
                                                                            subject = "id",
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)


fourclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 4,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2020_2023,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)

fiveclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                           mixture = drsc~1+year,
                                                                           random = ~1+year,
                                                                           ng = 5,
                                                                           nwg = FALSE, 
                                                                           data=coverage_long_2020_2023,
                                                                           subject = "id",
                                                                           B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)

sixclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                          mixture = drsc~1+year,
                                                                          random = ~1+year,
                                                                          ng = 6,
                                                                          nwg = FALSE, 
                                                                          data=coverage_long_2020_2023,
                                                                          subject = "id",
                                                                          B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)

sevenclass_linear_drsc_model_random_intercept_slope_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                                            mixture = drsc~1+year,
                                                                            random = ~1+year,
                                                                            ng = 7,
                                                                            nwg = FALSE, 
                                                                            data=coverage_long_2020_2023,
                                                                            subject = "id",  
                                                                            B= oneclass_linear_drsc_model_random_intercept_slope_2020_2023)


save(oneclass_linear_drsc_model_random_intercept_slope_2020_2023,
     twoclass_linear_drsc_model_random_intercept_slope_2020_2023,
     threeclass_linear_drsc_model_random_intercept_slope_2020_2023,
     fourclass_linear_drsc_model_random_intercept_slope_2020_2023,
     fiveclass_linear_drsc_model_random_intercept_slope_2020_2023,
     sixclass_linear_drsc_model_random_intercept_slope_2020_2023,
     sevenclass_linear_drsc_model_random_intercept_slope_2020_2023, file= "linear_drsc_model_random_intercept_slope_2020_2023.RData")

## Quadratic random effects model Equal ---------------------------------------


oneclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~ 1 + year, 
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data=coverage_long_2020_2023,
                                                      subject = "id")


twoclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 2,                                 
                                                      nwg = FALSE,   
                                                      idiag=FALSE,
                                                      data = coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2020_2023
)


threeclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 3,                                 
                                                        nwg = FALSE,
                                                        idiag=FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2020_2023
)



fourclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 4,                                 
                                                       nwg = FALSE, 
                                                       idiag=FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2020_2023
)


fiveclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                       mixture = ~ 1 + year + I(year^2), 
                                                       random = ~ 1 + year,       
                                                       ng = 5,                                 
                                                       nwg = FALSE,  
                                                       idiag=FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2020_2023
)

sixclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                      mixture = ~ 1 + year + I(year^2), 
                                                      random = ~ 1 + year,       
                                                      ng = 6,                                 
                                                      nwg = FALSE, 
                                                      idiag=FALSE,
                                                      data = coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2020_2023
)

sevenclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                        mixture = ~ 1 + year + I(year^2), 
                                                        random = ~ 1 + year,       
                                                        ng = 7,                                 
                                                        nwg = FALSE,    
                                                        idiag=FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2020_2023
)


save(oneclass_quadratic_drsc_model_2020_2023,
     twoclass_quadratic_drsc_model_2020_2023,
     threeclass_quadratic_drsc_model_2020_2023,
     fourclass_quadratic_drsc_model_2020_2023,
     fiveclass_quadratic_drsc_model_2020_2023,
     sixclass_quadratic_drsc_model_2020_2023,
     sevenclass_quadratic_drsc_model_2020_2023, file= "quadratic_drsc_model_2020_2023.RData")



## Quadratic random effects model Proportional ---------------------------------------


oneclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                           #mixture = ~1+year+I(year^2),
                                                           random = ~ 1 + year, 
                                                           ng = 1,
                                                           #nwg = TRUE, 
                                                           idiag=FALSE,
                                                           data=coverage_long_2020_2023,
                                                           subject = "id")


twoclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 2,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2020_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2020_2023
)



threeclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 3,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2020_2023,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2020_2023
)

fourclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 4,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2020_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2020_2023
)



fiveclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                            mixture = ~ 1 + year + I(year^2), 
                                                            random = ~ 1 + year,       
                                                            ng = 5,                                 
                                                            nwg = TRUE,   
                                                            idiag=FALSE,
                                                            data = coverage_long_2020_2023,
                                                            subject = "id",
                                                            B=oneclass_quadratic_prop_drsc_model_2020_2023
)

sixclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                           mixture = ~ 1 + year + I(year^2), 
                                                           random = ~ 1 + year,       
                                                           ng = 6,                                 
                                                           nwg = TRUE,   
                                                           idiag=FALSE,
                                                           data = coverage_long_2020_2023,
                                                           subject = "id",
                                                           B=oneclass_quadratic_prop_drsc_model_2020_2023
)



sevenclass_quadratic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2), 
                                                             mixture = ~ 1 + year + I(year^2), 
                                                             random = ~ 1 + year,       
                                                             ng = 7,                                 
                                                             nwg = TRUE,   
                                                             idiag=FALSE,
                                                             data = coverage_long_2020_2023,
                                                             subject = "id",
                                                             B=oneclass_quadratic_prop_drsc_model_2020_2023
)


save(oneclass_quadratic_prop_drsc_model_2020_2023,
     twoclass_quadratic_prop_drsc_model_2020_2023,
     threeclass_quadratic_prop_drsc_model_2020_2023,
     fourclass_quadratic_prop_drsc_model_2020_2023,
     fiveclass_quadratic_prop_drsc_model_2020_2023,
     sixclass_quadratic_prop_drsc_model_2020_2023,
     sevenclass_quadratic_prop_drsc_model_2020_2023, file= "quadratic_prop_drsc_model_2020_2023.RData")

## Cubic random effects model Equal ---------------------------------------


oneclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                  #mixture = ~1+year+I(year^2)+I(year^3),
                                                  random = ~ 1 + year, 
                                                  ng = 1,
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data=coverage_long_2020_2023,
                                                  subject = "id")


twoclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 2,                                 
                                                  nwg = FALSE,   
                                                  idiag=FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2020_2023
)


threeclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 3,                                 
                                                    nwg = FALSE,
                                                    idiag=FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2020_2023
)



fourclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 4,                                 
                                                   nwg = FALSE, 
                                                   idiag=FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2020_2023
)


fiveclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                   mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                   random = ~ 1 + year,       
                                                   ng = 5,                                 
                                                   nwg = FALSE,  
                                                   idiag=FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B=oneclass_cubic_drsc_model_2020_2023
)

sixclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                  mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                  random = ~ 1 + year,       
                                                  ng = 6,                                 
                                                  nwg = FALSE, 
                                                  idiag=FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B=oneclass_cubic_drsc_model_2020_2023
)

sevenclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                    mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                    random = ~ 1 + year,       
                                                    ng = 7,                                 
                                                    nwg = FALSE,    
                                                    idiag=FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B=oneclass_cubic_drsc_model_2020_2023
)



save(oneclass_cubic_drsc_model_2020_2023,
     twoclass_cubic_drsc_model_2020_2023,
     threeclass_cubic_drsc_model_2020_2023,
     fourclass_cubic_drsc_model_2020_2023,
     fiveclass_cubic_drsc_model_2020_2023,
     sixclass_cubic_drsc_model_2020_2023,
     sevenclass_cubic_drsc_model_2020_2023, file= "cubic_drsc_model_2020_2023.RData")


## Cubic random effects model Proportional ---------------------------------------


oneclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2)+I(year^3),
                                                       #mixture = ~1+year+I(year^2)+I(year^3),
                                                       random = ~ 1 + year, 
                                                       ng = 1,
                                                       #nwg = TRUE, 
                                                       idiag=FALSE,
                                                       data=coverage_long_2020_2023,
                                                       subject = "id")


twoclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 2,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2020_2023
)



threeclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 3,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2020_2023,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2020_2023
)

fourclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 4,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2020_2023
)



fiveclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                        mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                        random = ~ 1 + year,       
                                                        ng = 5,                                 
                                                        nwg = TRUE,   
                                                        idiag=FALSE,
                                                        data = coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_cubic_prop_drsc_model_2020_2023
)

sixclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                       mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                       random = ~ 1 + year,       
                                                       ng = 6,                                 
                                                       nwg = TRUE,   
                                                       idiag=FALSE,
                                                       data = coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_cubic_prop_drsc_model_2020_2023
)



sevenclass_cubic_prop_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2)+I(year^3), 
                                                         mixture = ~ 1 + year + I(year^2)+I(year^3), 
                                                         random = ~ 1 + year,       
                                                         ng = 7,                                 
                                                         nwg = TRUE,   
                                                         idiag=FALSE,
                                                         data = coverage_long_2020_2023,
                                                         subject = "id",
                                                         B=oneclass_cubic_prop_drsc_model_2020_2023
)

save(oneclass_cubic_prop_drsc_model_2020_2023,
     twoclass_cubic_prop_drsc_model_2020_2023,
     threeclass_cubic_prop_drsc_model_2020_2023,
     fourclass_cubic_prop_drsc_model_2020_2023,
     fiveclass_cubic_prop_drsc_model_2020_2023,
     sixclass_cubic_prop_drsc_model_2020_2023,
     sevenclass_cubic_prop_drsc_model_2020_2023, file= "cubic_prop_drsc_model_2020_2023.RData")















residualplot_step1 <- function(model, nameofoutcome, nameofage, data, 
                               ylimit = c(-50, 50), save_path = NULL, model_name = NULL) {
  require(dplyr)
  require(ggplot2)
  require(gridExtra)  # for grid.arrange
  
  if (is.null(model_name)) {
    stop("Model name must be provided.")
  }
  
  k <- model$ng
  preds <- model$pred
  names(preds)[6] <- nameofoutcome
  nameofid <- names(model$pred)[1]
  test <- dplyr::left_join(preds, model$pprob, by = nameofid)
  test <- dplyr::left_join(test, data, by = c(nameofid, nameofoutcome))
  
  plot_list <- list()  # List to store plots
  
  xlim_range <- range(data[[nameofage]], na.rm = TRUE)
  
  # Generate plots for each class
  for (i in 1:k) {
    newplotvalues <- test %>% 
      filter(class == i) %>% 
      mutate(Residuals = get(nameofoutcome) - eval(parse(text = paste0("pred_ss", i))))
    
    plotvaluessub <- newplotvalues
    
    p <- ggplot2::ggplot(data = plotvaluessub, aes(x = get(nameofage), y = Residuals, group = class)) +
      theme(axis.text = element_text(size = 8), text = element_text(size = 8)) + 
      geom_point(size = 0.1) + 
      stat_summary(fun = mean, geom = "line", size = 1, col = "CadetBlue", group = 1) + 
      ggtitle(paste("Residuals in class", i)) + 
      ylim(ylimit) + 
      xlim(xlim_range) + 
      labs(x = "Year") + 
      #coord_fixed(ratio = 1) + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
    
    plot_list[[i]] <- p
  }
  
  # Handle empty space if the last row has fewer than 4 plots
  total_plots <- length(plot_list)
  max_cols <- 7
  if (total_plots %% max_cols != 0) {
    empty_plots <- max_cols - (total_plots %% max_cols)
    for (j in seq_len(empty_plots)) {
      plot_list[[total_plots + j]] <- ggplot() + 
        theme_void() + 
        theme(panel.border = element_rect(colour = "white"))
    }
  }
  
  # Combine plots into a grid
  combined_plot <- grid.arrange(grobs = plot_list, ncol = max_cols, nrow = ceiling(total_plots / max_cols))
  
  # Save the combined plot if a save path is provided
  if (!is.null(save_path)) {
    g <- arrangeGrob(grobs = plot_list, ncol = max_cols, nrow = ceiling(total_plots / max_cols))
    ggsave(filename = file.path(save_path, paste0(model_name, ".jpeg")), plot = g,
           width = 21, height = 3, dpi = 1200)
  }
  
  return(combined_plot)  # Return the combined plot
}



# Apply the function to all models
linear_nre_homocedastic_drsc_model_2011_2023 <- list(
  oneclass_linear_nre_homocedastic_drsc_model_2011_2023,
  twoclass_linear_nre_homocedastic_drsc_model_2011_2023,
  threeclass_linear_nre_homocedastic_drsc_model_2011_2023,
  fourclass_linear_nre_homocedastic_drsc_model_2011_2023,
  fiveclass_linear_nre_homocedastic_drsc_model_2011_2023,
  sixclass_linear_nre_homocedastic_drsc_model_2011_2023,
  sevenclass_linear_nre_homocedastic_drsc_model_2011_2023
)

model_names_linear_nre_homocedastic_drsc_model_2011_2023 <- c(
  "oneclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "twoclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "threeclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "fourclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "fiveclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "sixclass_linear_nre_homocedastic_drsc_model_2011_2023",
  "sevenclass_linear_nre_homocedastic_drsc_model_2011_2023"
)

for (i in seq_along(linear_nre_homocedastic_drsc_model_2011_2023)) {
  residualplot_step1(
    model = linear_nre_homocedastic_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_homocedastic_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
linear_nre_heterocedastic_drsc_model_2011_2023 <- list(
  oneclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  twoclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  threeclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  fourclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  sixclass_linear_nre_heterocedastic_drsc_model_2011_2023,
  sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023
)

model_names_linear_nre_heterocedastic_drsc_model_2011_2023 <- c(
  "oneclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "twoclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "threeclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "fourclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "sixclass_linear_nre_heterocedastic_drsc_model_2011_2023",
  "sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023"
)

for (i in seq_along(linear_nre_heterocedastic_drsc_model_2011_2023)) {
  residualplot_step1(
    model = linear_nre_heterocedastic_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_heterocedastic_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
quadratic_nre_drsc_model_2011_2023 <- list(
  oneclass_quadratic_nre_drsc_model_2011_2023,
  twoclass_quadratic_nre_drsc_model_2011_2023,
  threeclass_quadratic_nre_drsc_model_2011_2023,
  fourclass_quadratic_nre_drsc_model_2011_2023,
  fiveclass_quadratic_nre_drsc_model_2011_2023,
  sixclass_quadratic_nre_drsc_model_2011_2023,
  sevenclass_quadratic_nre_drsc_model_2011_2023
)

model_names_quadratic_nre_drsc_model_2011_2023 <- c(
  "oneclass_quadratic_nre_drsc_model_2011_2023",
  "twoclass_quadratic_nre_drsc_model_2011_2023",
  "threeclass_quadratic_nre_drsc_model_2011_2023",
  "fourclass_quadratic_nre_drsc_model_2011_2023",
  "fiveclass_quadratic_nre_drsc_model_2011_2023",
  "sixclass_quadratic_nre_drsc_model_2011_2023",
  "sevenclass_quadratic_nre_drsc_model_2011_2023"
)

for (i in seq_along(quadratic_nre_drsc_model_2011_2023)) {
  residualplot_step1(
    model = quadratic_nre_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_nre_drsc_model_2011_2023[i]
  )
}


# Apply the function to all models
cubic_nre_drsc_model_2011_2023 <- list(
  oneclass_cubic_nre_drsc_model_2011_2023,
  twoclass_cubic_nre_drsc_model_2011_2023,
  threeclass_cubic_nre_drsc_model_2011_2023,
  fourclass_cubic_nre_drsc_model_2011_2023,
  fiveclass_cubic_nre_drsc_model_2011_2023,
  sixclass_cubic_nre_drsc_model_2011_2023,
  sevenclass_cubic_nre_drsc_model_2011_2023
)

model_names_cubic_nre_drsc_model_2011_2023 <- c(
  "oneclass_cubic_nre_drsc_model_2011_2023",
  "twoclass_cubic_nre_drsc_model_2011_2023",
  "threeclass_cubic_nre_drsc_model_2011_2023",
  "fourclass_cubic_nre_drsc_model_2011_2023",
  "fiveclass_cubic_nre_drsc_model_2011_2023",
  "sixclass_cubic_nre_drsc_model_2011_2023",
  "sevenclass_cubic_nre_drsc_model_2011_2023"
)

for (i in seq_along(cubic_nre_drsc_model_2011_2023)) {
  residualplot_step1(
    model = cubic_nre_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_nre_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
linear_drsc_model_random_intercept_2011_2023 <- list(
  oneclass_linear_drsc_model_random_intercept_2011_2023,
  twoclass_linear_drsc_model_random_intercept_2011_2023,
  threeclass_linear_drsc_model_random_intercept_2011_2023,
  fourclass_linear_drsc_model_random_intercept_2011_2023,
  fiveclass_linear_drsc_model_random_intercept_2011_2023,
  sixclass_linear_drsc_model_random_intercept_2011_2023,
  sevenclass_linear_drsc_model_random_intercept_2011_2023
)

model_names_linear_drsc_model_random_intercept_2011_2023 <- c(
  "oneclass_linear_drsc_model_random_intercept_2011_2023",
  "twoclass_linear_drsc_model_random_intercept_2011_2023",
  "threeclass_linear_drsc_model_random_intercept_2011_2023",
  "fourclass_linear_drsc_model_random_intercept_2011_2023",
  "fiveclass_linear_drsc_model_random_intercept_2011_2023",
  "sixclass_linear_drsc_model_random_intercept_2011_2023",
  "sevenclass_linear_drsc_model_random_intercept_2011_2023"
)

for (i in seq_along(linear_drsc_model_random_intercept_2011_2023)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_2011_2023[i]
  )
}




# Apply the function to all models
linear_drsc_model_random_intercept_slope_2011_2023 <- list(
  oneclass_linear_drsc_model_random_intercept_slope_2011_2023,
  twoclass_linear_drsc_model_random_intercept_slope_2011_2023,
  threeclass_linear_drsc_model_random_intercept_slope_2011_2023,
  fourclass_linear_drsc_model_random_intercept_slope_2011_2023,
  fiveclass_linear_drsc_model_random_intercept_slope_2011_2023,
  sixclass_linear_drsc_model_random_intercept_slope_2011_2023,
  sevenclass_linear_drsc_model_random_intercept_slope_2011_2023
)

model_names_linear_drsc_model_random_intercept_slope_2011_2023 <- c(
  "oneclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "twoclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "threeclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "fourclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "fiveclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "sixclass_linear_drsc_model_random_intercept_slope_2011_2023",
  "sevenclass_linear_drsc_model_random_intercept_slope_2011_2023"
)

for (i in seq_along(linear_drsc_model_random_intercept_slope_2011_2023)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_slope_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_slope_2011_2023[i]
  )
}


# Apply the function to all models
quadratic_drsc_model_2011_2023 <- list(
  oneclass_quadratic_drsc_model_2011_2023,
  twoclass_quadratic_drsc_model_2011_2023,
  threeclass_quadratic_drsc_model_2011_2023,
  fourclass_quadratic_drsc_model_2011_2023,
  fiveclass_quadratic_drsc_model_2011_2023,
  sixclass_quadratic_drsc_model_2011_2023,
  sevenclass_quadratic_drsc_model_2011_2023
)

model_names_quadratic_drsc_model_2011_2023 <- c(
  "oneclass_quadratic_drsc_model_2011_2023",
  "twoclass_quadratic_drsc_model_2011_2023",
  "threeclass_quadratic_drsc_model_2011_2023",
  "fourclass_quadratic_drsc_model_2011_2023",
  "fiveclass_quadratic_drsc_model_2011_2023",
  "sixclass_quadratic_drsc_model_2011_2023",
  "sevenclass_quadratic_drsc_model_2011_2023"
)

for (i in seq_along(quadratic_drsc_model_2011_2023)) {
  residualplot_step1(
    model = quadratic_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_drsc_model_2011_2023[i]
  )
}




# Apply the function to all models
quadratic_prop_drsc_model_2011_2023 <- list(
  oneclass_quadratic_prop_drsc_model_2011_2023,
  twoclass_quadratic_prop_drsc_model_2011_2023,
  threeclass_quadratic_prop_drsc_model_2011_2023,
  fourclass_quadratic_prop_drsc_model_2011_2023,
  fiveclass_quadratic_prop_drsc_model_2011_2023,
  sixclass_quadratic_prop_drsc_model_2011_2023,
  sevenclass_quadratic_prop_drsc_model_2011_2023
)

model_names_quadratic_prop_drsc_model_2011_2023 <- c(
  "oneclass_quadratic_prop_drsc_model_2011_2023",
  "twoclass_quadratic_prop_drsc_model_2011_2023",
  "threeclass_quadratic_prop_drsc_model_2011_2023",
  "fourclass_quadratic_prop_drsc_model_2011_2023",
  "fiveclass_quadratic_prop_drsc_model_2011_2023",
  "sixclass_quadratic_prop_drsc_model_2011_2023",
  "sevenclass_quadratic_prop_drsc_model_2011_2023"
)

for (i in seq_along(quadratic_prop_drsc_model_2011_2023)) {
  residualplot_step1(
    model = quadratic_prop_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_prop_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
cubic_drsc_model_2011_2023 <- list(
  oneclass_cubic_drsc_model_2011_2023,
  twoclass_cubic_drsc_model_2011_2023,
  threeclass_cubic_drsc_model_2011_2023,
  fourclass_cubic_drsc_model_2011_2023,
  fiveclass_cubic_drsc_model_2011_2023,
  sixclass_cubic_drsc_model_2011_2023,
  sevenclass_cubic_drsc_model_2011_2023
)

model_names_cubic_drsc_model_2011_2023 <- c(
  "oneclass_cubic_drsc_model_2011_2023",
  "twoclass_cubic_drsc_model_2011_2023",
  "threeclass_cubic_drsc_model_2011_2023",
  "fourclass_cubic_drsc_model_2011_2023",
  "fiveclass_cubic_drsc_model_2011_2023",
  "sixclass_cubic_drsc_model_2011_2023",
  "sevenclass_cubic_drsc_model_2011_2023"
)

for (i in seq_along(cubic_drsc_model_2011_2023)) {
  residualplot_step1(
    model = cubic_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
cubic_prop_drsc_model_2011_2023 <- list(
  oneclass_cubic_prop_drsc_model_2011_2023,
  twoclass_cubic_prop_drsc_model_2011_2023,
  threeclass_cubic_prop_drsc_model_2011_2023,
  fourclass_cubic_prop_drsc_model_2011_2023,
  fiveclass_cubic_prop_drsc_model_2011_2023,
  sixclass_cubic_prop_drsc_model_2011_2023,
  sevenclass_cubic_prop_drsc_model_2011_2023
)

model_names_cubic_prop_drsc_model_2011_2023 <- c(
  "oneclass_cubic_prop_drsc_model_2011_2023",
  "twoclass_cubic_prop_drsc_model_2011_2023",
  "threeclass_cubic_prop_drsc_model_2011_2023",
  "fourclass_cubic_prop_drsc_model_2011_2023",
  "fiveclass_cubic_prop_drsc_model_2011_2023",
  "sixclass_cubic_prop_drsc_model_2011_2023",
  "sevenclass_cubic_prop_drsc_model_2011_2023"
)

for (i in seq_along(cubic_prop_drsc_model_2011_2023)) {
  residualplot_step1(
    model = cubic_prop_drsc_model_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_prop_drsc_model_2011_2023[i]
  )
}



# Apply the function to all models
linear_nre_homocedastic_drsc_model_2011_2019 <- list(
  oneclass_linear_nre_homocedastic_drsc_model_2011_2019,
  twoclass_linear_nre_homocedastic_drsc_model_2011_2019,
  threeclass_linear_nre_homocedastic_drsc_model_2011_2019,
  fourclass_linear_nre_homocedastic_drsc_model_2011_2019,
  fiveclass_linear_nre_homocedastic_drsc_model_2011_2019,
  sixclass_linear_nre_homocedastic_drsc_model_2011_2019,
  sevenclass_linear_nre_homocedastic_drsc_model_2011_2019
)

model_names_linear_nre_homocedastic_drsc_model_2011_2019 <- c(
  "oneclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "twoclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "threeclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "fourclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "fiveclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "sixclass_linear_nre_homocedastic_drsc_model_2011_2019",
  "sevenclass_linear_nre_homocedastic_drsc_model_2011_2019"
)

for (i in seq_along(linear_nre_homocedastic_drsc_model_2011_2019)) {
  residualplot_step1(
    model = linear_nre_homocedastic_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_homocedastic_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
linear_nre_heterocedastic_drsc_model_2011_2019 <- list(
  oneclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  twoclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  threeclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  fourclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  sixclass_linear_nre_heterocedastic_drsc_model_2011_2019,
  sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019
)

model_names_linear_nre_heterocedastic_drsc_model_2011_2019 <- c(
  "oneclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "twoclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "threeclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "fourclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "sixclass_linear_nre_heterocedastic_drsc_model_2011_2019",
  "sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019"
)

for (i in seq_along(linear_nre_heterocedastic_drsc_model_2011_2019)) {
  residualplot_step1(
    model = linear_nre_heterocedastic_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_heterocedastic_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
quadratic_nre_drsc_model_2011_2019 <- list(
  oneclass_quadratic_nre_drsc_model_2011_2019,
  twoclass_quadratic_nre_drsc_model_2011_2019,
  threeclass_quadratic_nre_drsc_model_2011_2019,
  fourclass_quadratic_nre_drsc_model_2011_2019,
  fiveclass_quadratic_nre_drsc_model_2011_2019,
  sixclass_quadratic_nre_drsc_model_2011_2019,
  sevenclass_quadratic_nre_drsc_model_2011_2019
)

model_names_quadratic_nre_drsc_model_2011_2019 <- c(
  "oneclass_quadratic_nre_drsc_model_2011_2019",
  "twoclass_quadratic_nre_drsc_model_2011_2019",
  "threeclass_quadratic_nre_drsc_model_2011_2019",
  "fourclass_quadratic_nre_drsc_model_2011_2019",
  "fiveclass_quadratic_nre_drsc_model_2011_2019",
  "sixclass_quadratic_nre_drsc_model_2011_2019",
  "sevenclass_quadratic_nre_drsc_model_2011_2019"
)

for (i in seq_along(quadratic_nre_drsc_model_2011_2019)) {
  residualplot_step1(
    model = quadratic_nre_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_nre_drsc_model_2011_2019[i]
  )
}


# Apply the function to all models
cubic_nre_drsc_model_2011_2019 <- list(
  oneclass_cubic_nre_drsc_model_2011_2019,
  twoclass_cubic_nre_drsc_model_2011_2019,
  threeclass_cubic_nre_drsc_model_2011_2019,
  fourclass_cubic_nre_drsc_model_2011_2019,
  fiveclass_cubic_nre_drsc_model_2011_2019,
  sixclass_cubic_nre_drsc_model_2011_2019,
  sevenclass_cubic_nre_drsc_model_2011_2019
)

model_names_cubic_nre_drsc_model_2011_2019 <- c(
  "oneclass_cubic_nre_drsc_model_2011_2019",
  "twoclass_cubic_nre_drsc_model_2011_2019",
  "threeclass_cubic_nre_drsc_model_2011_2019",
  "fourclass_cubic_nre_drsc_model_2011_2019",
  "fiveclass_cubic_nre_drsc_model_2011_2019",
  "sixclass_cubic_nre_drsc_model_2011_2019",
  "sevenclass_cubic_nre_drsc_model_2011_2019"
)

for (i in seq_along(cubic_nre_drsc_model_2011_2019)) {
  residualplot_step1(
    model = cubic_nre_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_nre_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
linear_drsc_model_random_intercept_2011_2019 <- list(
  oneclass_linear_drsc_model_random_intercept_2011_2019,
  twoclass_linear_drsc_model_random_intercept_2011_2019,
  threeclass_linear_drsc_model_random_intercept_2011_2019,
  fourclass_linear_drsc_model_random_intercept_2011_2019,
  fiveclass_linear_drsc_model_random_intercept_2011_2019,
  sixclass_linear_drsc_model_random_intercept_2011_2019,
  sevenclass_linear_drsc_model_random_intercept_2011_2019
)

model_names_linear_drsc_model_random_intercept_2011_2019 <- c(
  "oneclass_linear_drsc_model_random_intercept_2011_2019",
  "twoclass_linear_drsc_model_random_intercept_2011_2019",
  "threeclass_linear_drsc_model_random_intercept_2011_2019",
  "fourclass_linear_drsc_model_random_intercept_2011_2019",
  "fiveclass_linear_drsc_model_random_intercept_2011_2019",
  "sixclass_linear_drsc_model_random_intercept_2011_2019",
  "sevenclass_linear_drsc_model_random_intercept_2011_2019"
)

for (i in seq_along(linear_drsc_model_random_intercept_2011_2019)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_2011_2019[i]
  )
}




# Apply the function to all models
linear_drsc_model_random_intercept_slope_2011_2019 <- list(
  oneclass_linear_drsc_model_random_intercept_slope_2011_2019,
  twoclass_linear_drsc_model_random_intercept_slope_2011_2019,
  threeclass_linear_drsc_model_random_intercept_slope_2011_2019,
  fourclass_linear_drsc_model_random_intercept_slope_2011_2019,
  fiveclass_linear_drsc_model_random_intercept_slope_2011_2019,
  sixclass_linear_drsc_model_random_intercept_slope_2011_2019,
  sevenclass_linear_drsc_model_random_intercept_slope_2011_2019
)

model_names_linear_drsc_model_random_intercept_slope_2011_2019 <- c(
  "oneclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "twoclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "threeclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "fourclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "fiveclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "sixclass_linear_drsc_model_random_intercept_slope_2011_2019",
  "sevenclass_linear_drsc_model_random_intercept_slope_2011_2019"
)

for (i in seq_along(linear_drsc_model_random_intercept_slope_2011_2019)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_slope_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_slope_2011_2019[i]
  )
}


# Apply the function to all models
quadratic_drsc_model_2011_2019 <- list(
  oneclass_quadratic_drsc_model_2011_2019,
  twoclass_quadratic_drsc_model_2011_2019,
  threeclass_quadratic_drsc_model_2011_2019,
  fourclass_quadratic_drsc_model_2011_2019,
  fiveclass_quadratic_drsc_model_2011_2019,
  sixclass_quadratic_drsc_model_2011_2019,
  sevenclass_quadratic_drsc_model_2011_2019
)

model_names_quadratic_drsc_model_2011_2019 <- c(
  "oneclass_quadratic_drsc_model_2011_2019",
  "twoclass_quadratic_drsc_model_2011_2019",
  "threeclass_quadratic_drsc_model_2011_2019",
  "fourclass_quadratic_drsc_model_2011_2019",
  "fiveclass_quadratic_drsc_model_2011_2019",
  "sixclass_quadratic_drsc_model_2011_2019",
  "sevenclass_quadratic_drsc_model_2011_2019"
)

for (i in seq_along(quadratic_drsc_model_2011_2019)) {
  residualplot_step1(
    model = quadratic_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_drsc_model_2011_2019[i]
  )
}




# Apply the function to all models
quadratic_prop_drsc_model_2011_2019 <- list(
  oneclass_quadratic_prop_drsc_model_2011_2019,
  twoclass_quadratic_prop_drsc_model_2011_2019,
  threeclass_quadratic_prop_drsc_model_2011_2019,
  fourclass_quadratic_prop_drsc_model_2011_2019,
  fiveclass_quadratic_prop_drsc_model_2011_2019,
  sixclass_quadratic_prop_drsc_model_2011_2019,
  sevenclass_quadratic_prop_drsc_model_2011_2019
)

model_names_quadratic_prop_drsc_model_2011_2019 <- c(
  "oneclass_quadratic_prop_drsc_model_2011_2019",
  "twoclass_quadratic_prop_drsc_model_2011_2019",
  "threeclass_quadratic_prop_drsc_model_2011_2019",
  "fourclass_quadratic_prop_drsc_model_2011_2019",
  "fiveclass_quadratic_prop_drsc_model_2011_2019",
  "sixclass_quadratic_prop_drsc_model_2011_2019",
  "sevenclass_quadratic_prop_drsc_model_2011_2019"
)

for (i in seq_along(quadratic_prop_drsc_model_2011_2019)) {
  residualplot_step1(
    model = quadratic_prop_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_prop_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
cubic_drsc_model_2011_2019 <- list(
  oneclass_cubic_drsc_model_2011_2019,
  twoclass_cubic_drsc_model_2011_2019,
  threeclass_cubic_drsc_model_2011_2019,
  fourclass_cubic_drsc_model_2011_2019,
  fiveclass_cubic_drsc_model_2011_2019,
  sixclass_cubic_drsc_model_2011_2019,
  sevenclass_cubic_drsc_model_2011_2019
)

model_names_cubic_drsc_model_2011_2019 <- c(
  "oneclass_cubic_drsc_model_2011_2019",
  "twoclass_cubic_drsc_model_2011_2019",
  "threeclass_cubic_drsc_model_2011_2019",
  "fourclass_cubic_drsc_model_2011_2019",
  "fiveclass_cubic_drsc_model_2011_2019",
  "sixclass_cubic_drsc_model_2011_2019",
  "sevenclass_cubic_drsc_model_2011_2019"
)

for (i in seq_along(cubic_drsc_model_2011_2019)) {
  residualplot_step1(
    model = cubic_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
cubic_prop_drsc_model_2011_2019 <- list(
  oneclass_cubic_prop_drsc_model_2011_2019,
  twoclass_cubic_prop_drsc_model_2011_2019,
  threeclass_cubic_prop_drsc_model_2011_2019,
  fourclass_cubic_prop_drsc_model_2011_2019,
  fiveclass_cubic_prop_drsc_model_2011_2019,
  sixclass_cubic_prop_drsc_model_2011_2019,
  sevenclass_cubic_prop_drsc_model_2011_2019
)

model_names_cubic_prop_drsc_model_2011_2019 <- c(
  "oneclass_cubic_prop_drsc_model_2011_2019",
  "twoclass_cubic_prop_drsc_model_2011_2019",
  "threeclass_cubic_prop_drsc_model_2011_2019",
  "fourclass_cubic_prop_drsc_model_2011_2019",
  "fiveclass_cubic_prop_drsc_model_2011_2019",
  "sixclass_cubic_prop_drsc_model_2011_2019",
  "sevenclass_cubic_prop_drsc_model_2011_2019"
)

for (i in seq_along(cubic_prop_drsc_model_2011_2019)) {
  residualplot_step1(
    model = cubic_prop_drsc_model_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_prop_drsc_model_2011_2019[i]
  )
}



# Apply the function to all models
linear_nre_homocedastic_drsc_model_2020_2023 <- list(
  oneclass_linear_nre_homocedastic_drsc_model_2020_2023,
  twoclass_linear_nre_homocedastic_drsc_model_2020_2023,
  threeclass_linear_nre_homocedastic_drsc_model_2020_2023,
  fourclass_linear_nre_homocedastic_drsc_model_2020_2023,
  fiveclass_linear_nre_homocedastic_drsc_model_2020_2023,
  sixclass_linear_nre_homocedastic_drsc_model_2020_2023,
  sevenclass_linear_nre_homocedastic_drsc_model_2020_2023
)

model_names_linear_nre_homocedastic_drsc_model_2020_2023 <- c(
  "oneclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "twoclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "threeclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "fourclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "fiveclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "sixclass_linear_nre_homocedastic_drsc_model_2020_2023",
  "sevenclass_linear_nre_homocedastic_drsc_model_2020_2023"
)

for (i in seq_along(linear_nre_homocedastic_drsc_model_2020_2023)) {
  residualplot_step1(
    model = linear_nre_homocedastic_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_homocedastic_drsc_model_2020_2023[i]
  )
}



# Apply the function to all models
linear_nre_heterocedastic_drsc_model_2020_2023 <- list(
  oneclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  twoclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  threeclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  fourclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  fiveclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  sixclass_linear_nre_heterocedastic_drsc_model_2020_2023,
  sevenclass_linear_nre_heterocedastic_drsc_model_2020_2023
)

model_names_linear_nre_heterocedastic_drsc_model_2020_2023 <- c(
  "oneclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "twoclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "threeclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "fourclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "fiveclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "sixclass_linear_nre_heterocedastic_drsc_model_2020_2023",
  "sevenclass_linear_nre_heterocedastic_drsc_model_2020_2023"
)

for (i in seq_along(linear_nre_heterocedastic_drsc_model_2020_2023)) {
  residualplot_step1(
    model = linear_nre_heterocedastic_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_nre_heterocedastic_drsc_model_2020_2023[i]
  )
}



# Apply the function to all models
quadratic_nre_drsc_model_2020_2023 <- list(
  oneclass_quadratic_nre_drsc_model_2020_2023,
  twoclass_quadratic_nre_drsc_model_2020_2023,
  threeclass_quadratic_nre_drsc_model_2020_2023,
  fourclass_quadratic_nre_drsc_model_2020_2023,
  fiveclass_quadratic_nre_drsc_model_2020_2023,
  sixclass_quadratic_nre_drsc_model_2020_2023,
  sevenclass_quadratic_nre_drsc_model_2020_2023
)

model_names_quadratic_nre_drsc_model_2020_2023 <- c(
  "oneclass_quadratic_nre_drsc_model_2020_2023",
  "twoclass_quadratic_nre_drsc_model_2020_2023",
  "threeclass_quadratic_nre_drsc_model_2020_2023",
  "fourclass_quadratic_nre_drsc_model_2020_2023",
  "fiveclass_quadratic_nre_drsc_model_2020_2023",
  "sixclass_quadratic_nre_drsc_model_2020_2023",
  "sevenclass_quadratic_nre_drsc_model_2020_2023"
)

for (i in seq_along(quadratic_nre_drsc_model_2020_2023)) {
  residualplot_step1(
    model = quadratic_nre_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_nre_drsc_model_2020_2023[i]
  )
}


# Apply the function to all models
cubic_nre_drsc_model_2020_2023 <- list(
  oneclass_cubic_nre_drsc_model_2020_2023,
  twoclass_cubic_nre_drsc_model_2020_2023,
  threeclass_cubic_nre_drsc_model_2020_2023,
  fourclass_cubic_nre_drsc_model_2020_2023,
  fiveclass_cubic_nre_drsc_model_2020_2023,
  sixclass_cubic_nre_drsc_model_2020_2023,
  sevenclass_cubic_nre_drsc_model_2020_2023
)

model_names_cubic_nre_drsc_model_2020_2023 <- c(
  "oneclass_cubic_nre_drsc_model_2020_2023",
  "twoclass_cubic_nre_drsc_model_2020_2023",
  "threeclass_cubic_nre_drsc_model_2020_2023",
  "fourclass_cubic_nre_drsc_model_2020_2023",
  "fiveclass_cubic_nre_drsc_model_2020_2023",
  "sixclass_cubic_nre_drsc_model_2020_2023",
  "sevenclass_cubic_nre_drsc_model_2020_2023"
)

for (i in seq_along(cubic_nre_drsc_model_2020_2023)) {
  residualplot_step1(
    model = cubic_nre_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_nre_drsc_model_2020_2023[i]
  )
}



# Apply the function to all models
linear_drsc_model_random_intercept_2020_2023 <- list(
  oneclass_linear_drsc_model_random_intercept_2020_2023,
  twoclass_linear_drsc_model_random_intercept_2020_2023,
  threeclass_linear_drsc_model_random_intercept_2020_2023,
  fourclass_linear_drsc_model_random_intercept_2020_2023,
  fiveclass_linear_drsc_model_random_intercept_2020_2023,
  sixclass_linear_drsc_model_random_intercept_2020_2023,
  sevenclass_linear_drsc_model_random_intercept_2020_2023
)

model_names_linear_drsc_model_random_intercept_2020_2023 <- c(
  "oneclass_linear_drsc_model_random_intercept_2020_2023",
  "twoclass_linear_drsc_model_random_intercept_2020_2023",
  "threeclass_linear_drsc_model_random_intercept_2020_2023",
  "fourclass_linear_drsc_model_random_intercept_2020_2023",
  "fiveclass_linear_drsc_model_random_intercept_2020_2023",
  "sixclass_linear_drsc_model_random_intercept_2020_2023",
  "sevenclass_linear_drsc_model_random_intercept_2020_2023"
)

for (i in seq_along(linear_drsc_model_random_intercept_2020_2023)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_2020_2023[i]
  )
}




# Apply the function to all models
linear_drsc_model_random_intercept_slope_2020_2023 <- list(
  oneclass_linear_drsc_model_random_intercept_slope_2020_2023,
  twoclass_linear_drsc_model_random_intercept_slope_2020_2023,
  threeclass_linear_drsc_model_random_intercept_slope_2020_2023,
  fourclass_linear_drsc_model_random_intercept_slope_2020_2023,
  fiveclass_linear_drsc_model_random_intercept_slope_2020_2023,
  sixclass_linear_drsc_model_random_intercept_slope_2020_2023,
  sevenclass_linear_drsc_model_random_intercept_slope_2020_2023
)

model_names_linear_drsc_model_random_intercept_slope_2020_2023 <- c(
  "oneclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "twoclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "threeclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "fourclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "fiveclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "sixclass_linear_drsc_model_random_intercept_slope_2020_2023",
  "sevenclass_linear_drsc_model_random_intercept_slope_2020_2023"
)

for (i in seq_along(linear_drsc_model_random_intercept_slope_2020_2023)) {
  residualplot_step1(
    model = linear_drsc_model_random_intercept_slope_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_linear_drsc_model_random_intercept_slope_2020_2023[i]
  )
}


# Apply the function to all models
quadratic_drsc_model_2020_2023 <- list(
  oneclass_quadratic_drsc_model_2020_2023,
  twoclass_quadratic_drsc_model_2020_2023,
  threeclass_quadratic_drsc_model_2020_2023,
  fourclass_quadratic_drsc_model_2020_2023,
  fiveclass_quadratic_drsc_model_2020_2023,
  sixclass_quadratic_drsc_model_2020_2023,
  sevenclass_quadratic_drsc_model_2020_2023
)

model_names_quadratic_drsc_model_2020_2023 <- c(
  "oneclass_quadratic_drsc_model_2020_2023",
  "twoclass_quadratic_drsc_model_2020_2023",
  "threeclass_quadratic_drsc_model_2020_2023",
  "fourclass_quadratic_drsc_model_2020_2023",
  "fiveclass_quadratic_drsc_model_2020_2023",
  "sixclass_quadratic_drsc_model_2020_2023",
  "sevenclass_quadratic_drsc_model_2020_2023"
)

for (i in seq_along(quadratic_drsc_model_2020_2023)) {
  residualplot_step1(
    model = quadratic_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_drsc_model_2020_2023[i]
  )
}




# Apply the function to all models
quadratic_prop_drsc_model_2020_2023 <- list(
  oneclass_quadratic_prop_drsc_model_2020_2023,
  twoclass_quadratic_prop_drsc_model_2020_2023,
  threeclass_quadratic_prop_drsc_model_2020_2023,
  fourclass_quadratic_prop_drsc_model_2020_2023,
  fiveclass_quadratic_prop_drsc_model_2020_2023,
  sixclass_quadratic_prop_drsc_model_2020_2023,
  sevenclass_quadratic_prop_drsc_model_2020_2023
)

model_names_quadratic_prop_drsc_model_2020_2023 <- c(
  "oneclass_quadratic_prop_drsc_model_2020_2023",
  "twoclass_quadratic_prop_drsc_model_2020_2023",
  "threeclass_quadratic_prop_drsc_model_2020_2023",
  "fourclass_quadratic_prop_drsc_model_2020_2023",
  "fiveclass_quadratic_prop_drsc_model_2020_2023",
  "sixclass_quadratic_prop_drsc_model_2020_2023",
  "sevenclass_quadratic_prop_drsc_model_2020_2023"
)

for (i in seq_along(quadratic_prop_drsc_model_2020_2023)) {
  residualplot_step1(
    model = quadratic_prop_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_quadratic_prop_drsc_model_2020_2023[i]
  )
}



# Apply the function to all models
cubic_drsc_model_2020_2023 <- list(
  oneclass_cubic_drsc_model_2020_2023,
  twoclass_cubic_drsc_model_2020_2023,
  threeclass_cubic_drsc_model_2020_2023,
  fourclass_cubic_drsc_model_2020_2023,
  fiveclass_cubic_drsc_model_2020_2023,
  sixclass_cubic_drsc_model_2020_2023,
  sevenclass_cubic_drsc_model_2020_2023
)

model_names_cubic_drsc_model_2020_2023 <- c(
  "oneclass_cubic_drsc_model_2020_2023",
  "twoclass_cubic_drsc_model_2020_2023",
  "threeclass_cubic_drsc_model_2020_2023",
  "fourclass_cubic_drsc_model_2020_2023",
  "fiveclass_cubic_drsc_model_2020_2023",
  "sixclass_cubic_drsc_model_2020_2023",
  "sevenclass_cubic_drsc_model_2020_2023"
)

for (i in seq_along(cubic_drsc_model_2020_2023)) {
  residualplot_step1(
    model = cubic_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_drsc_model_2020_2023[i]
  )
}



# Apply the function to all models
cubic_prop_drsc_model_2020_2023 <- list(
  oneclass_cubic_prop_drsc_model_2020_2023,
  twoclass_cubic_prop_drsc_model_2020_2023,
  threeclass_cubic_prop_drsc_model_2020_2023,
  fourclass_cubic_prop_drsc_model_2020_2023,
  fiveclass_cubic_prop_drsc_model_2020_2023,
  sixclass_cubic_prop_drsc_model_2020_2023,
  sevenclass_cubic_prop_drsc_model_2020_2023
)

model_names_cubic_prop_drsc_model_2020_2023 <- c(
  "oneclass_cubic_prop_drsc_model_2020_2023",
  "twoclass_cubic_prop_drsc_model_2020_2023",
  "threeclass_cubic_prop_drsc_model_2020_2023",
  "fourclass_cubic_prop_drsc_model_2020_2023",
  "fiveclass_cubic_prop_drsc_model_2020_2023",
  "sixclass_cubic_prop_drsc_model_2020_2023",
  "sevenclass_cubic_prop_drsc_model_2020_2023"
)

for (i in seq_along(cubic_prop_drsc_model_2020_2023)) {
  residualplot_step1(
    model = cubic_prop_drsc_model_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
      ylimit = c(-1, 1),
    save_path = getwd(),
    model_name = model_names_cubic_prop_drsc_model_2020_2023[i]
  )
}




# Step 2 - run a random effects model --------------------------------------------------

## Linear non random effects model DRSC 2011-2023 -----------------------------------------


summarytable_linear_nre_homocedastic_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              twoclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              threeclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              fourclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              fiveclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              sixclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                              sevenclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                                                                        which=c("G", 
                                                                                                "loglik", 
                                                                                                "conv", 
                                                                                                "npm", 
                                                                                                "AIC", 
                                                                                                "BIC", 
                                                                                                "SABIC", 
                                                                                                "entropy", 
                                                                                                "ICL", 
                                                                                                "ICL1", 
                                                                                                "ICL2", 
                                                                                                "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob(threeclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob(fourclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob(fiveclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob(sixclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob(sevenclass_linear_nre_homocedastic_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(twoclass_linear_nre_homocedastic_drsc_model_2011_2023)  
postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(threeclass_linear_nre_homocedastic_drsc_model_2011_2023)  
postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(fourclass_linear_nre_homocedastic_drsc_model_2011_2023)  
postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(fiveclass_linear_nre_homocedastic_drsc_model_2011_2023)  
postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(sixclass_linear_nre_homocedastic_drsc_model_2011_2023)
postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2023 <- postprob(sevenclass_linear_nre_homocedastic_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_nre_homocedastic_drsc_model_2011_2023 <- list(twoclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                               threeclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                               fourclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                               fiveclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                               sixclass_linear_nre_homocedastic_drsc_model_2011_2023,
                                               sevenclass_linear_nre_homocedastic_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_nre_homocedastic_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$loglik[7], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$npm[7], summarytable_linear_nre_homocedastic_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_nre_homocedastic_drsc_model_2011_2023 <- sapply(outputs_vllrt_linear_nre_homocedastic_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_linear_nre_homocedastic_drsc_model_2011_2023 <- summarytable_linear_nre_homocedastic_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_linear_nre_homocedastic_drsc_model_2011_2023$ns, min(postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_nre_homocedastic_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_nre_homocedastic_drsc_model_2011_2023, "summarytable_linear_nre_homocedastic_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_linear_nre_homocedastic_drsc_model_2011_2023)


## Linear non random effects model DRSC 2011-2023 -----------------------------------------


summarytable_linear_nre_heterocedastic_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               twoclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               threeclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               fourclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               sixclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                                               which=c("G", 
                                                                                                       "loglik", 
                                                                                                       "conv", 
                                                                                                       "npm", 
                                                                                                       "AIC", 
                                                                                                       "BIC", 
                                                                                                       "SABIC", 
                                                                                                       "entropy", 
                                                                                                       "ICL", 
                                                                                                       "ICL1", 
                                                                                                       "ICL2", 
                                                                                                       "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob(threeclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob(fourclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob(sixclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob(sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(twoclass_linear_nre_heterocedastic_drsc_model_2011_2023)  
postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(threeclass_linear_nre_heterocedastic_drsc_model_2011_2023)  
postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(fourclass_linear_nre_heterocedastic_drsc_model_2011_2023)  
postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023)  
postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(sixclass_linear_nre_heterocedastic_drsc_model_2011_2023)
postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023 <- postprob(sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_nre_heterocedastic_drsc_model_2011_2023 <- list(twoclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                threeclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                fourclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                sixclass_linear_nre_heterocedastic_drsc_model_2011_2023,
                                                                sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2023 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_nre_heterocedastic_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$loglik[7], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$npm[7], summarytable_linear_nre_heterocedastic_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_nre_heterocedastic_drsc_model_2011_2023 <- sapply(outputs_vllrt_linear_nre_heterocedastic_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_linear_nre_heterocedastic_drsc_model_2011_2023 <- summarytable_linear_nre_heterocedastic_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_linear_nre_heterocedastic_drsc_model_2011_2023$ns, min(postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_nre_heterocedastic_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_nre_heterocedastic_drsc_model_2011_2023, "summarytable_linear_nre_heterocedastic_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_linear_nre_heterocedastic_drsc_model_2011_2023)



## Linear random effects model DRSC 2011-2023 -----------------------------------------


summarytable_linear_drsc_model_random_intercept_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               twoclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               threeclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               fourclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               fiveclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               sixclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               sevenclass_linear_drsc_model_random_intercept_2011_2023,
                                                                                               which=c("G", 
                                                                                                       "loglik", 
                                                                                                       "conv", 
                                                                                                       "npm", 
                                                                                                       "AIC", 
                                                                                                       "BIC", 
                                                                                                       "SABIC", 
                                                                                                       "entropy", 
                                                                                                       "ICL", 
                                                                                                       "ICL1", 
                                                                                                       "ICL2", 
                                                                                                       "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_drsc_model_random_intercept_2011_2023)
postprob(threeclass_linear_drsc_model_random_intercept_2011_2023)
postprob(fourclass_linear_drsc_model_random_intercept_2011_2023)
postprob(fiveclass_linear_drsc_model_random_intercept_2011_2023)
postprob(sixclass_linear_drsc_model_random_intercept_2011_2023)
postprob(sevenclass_linear_drsc_model_random_intercept_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(twoclass_linear_drsc_model_random_intercept_2011_2023)  
postprob_threeclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(threeclass_linear_drsc_model_random_intercept_2011_2023)  
postprob_fourclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(fourclass_linear_drsc_model_random_intercept_2011_2023)  
postprob_fiveclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(fiveclass_linear_drsc_model_random_intercept_2011_2023)  
postprob_sixclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(sixclass_linear_drsc_model_random_intercept_2011_2023)
postprob_sevenclass_linear_drsc_model_random_intercept_2011_2023 <- postprob(sevenclass_linear_drsc_model_random_intercept_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_drsc_model_random_intercept_2011_2023 <- list(twoclass_linear_drsc_model_random_intercept_2011_2023,
                                                                threeclass_linear_drsc_model_random_intercept_2011_2023,
                                                                fourclass_linear_drsc_model_random_intercept_2011_2023,
                                                                fiveclass_linear_drsc_model_random_intercept_2011_2023,
                                                                sixclass_linear_drsc_model_random_intercept_2011_2023,
                                                                sevenclass_linear_drsc_model_random_intercept_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_drsc_model_random_intercept_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_drsc_model_random_intercept_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_drsc_model_random_intercept_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_drsc_model_random_intercept_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_drsc_model_random_intercept_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_drsc_model_random_intercept_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_drsc_model_random_intercept_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[1], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[1], summarytable_linear_drsc_model_random_intercept_2011_2023$G[1], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[2], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[2], summarytable_linear_drsc_model_random_intercept_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[2], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[2], summarytable_linear_drsc_model_random_intercept_2011_2023$G[2], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[3], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[3], summarytable_linear_drsc_model_random_intercept_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[3], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[3], summarytable_linear_drsc_model_random_intercept_2011_2023$G[3], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[4], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[4], summarytable_linear_drsc_model_random_intercept_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[4], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[4], summarytable_linear_drsc_model_random_intercept_2011_2023$G[4], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[5], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[5], summarytable_linear_drsc_model_random_intercept_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[5], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[5], summarytable_linear_drsc_model_random_intercept_2011_2023$G[5], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[6], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[6], summarytable_linear_drsc_model_random_intercept_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_linear_drsc_model_random_intercept_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[6], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[6], summarytable_linear_drsc_model_random_intercept_2011_2023$G[6], summarytable_linear_drsc_model_random_intercept_2011_2023$loglik[7], summarytable_linear_drsc_model_random_intercept_2011_2023$npm[7], summarytable_linear_drsc_model_random_intercept_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_drsc_model_random_intercept_2011_2023 <- sapply(outputs_vllrt_linear_drsc_model_random_intercept_2011_2023, function(x) gsub("^= ", "", x))



summarytable_linear_drsc_model_random_intercept_2011_2023 <- summarytable_linear_drsc_model_random_intercept_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,]), min(postprob_threeclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,]), min(postprob_fourclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,]), min(postprob_sixclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_linear_drsc_model_random_intercept_2011_2023$ns, min(postprob_twoclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]), min(postprob_threeclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]), min(postprob_fourclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]), min(postprob_sixclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_drsc_model_random_intercept_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_drsc_model_random_intercept_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_linear_drsc_model_random_intercept_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_drsc_model_random_intercept_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_drsc_model_random_intercept_2011_2023, "summarytable_linear_drsc_model_random_intercept_2011_2023.csv", row.names = F)

View(summarytable_linear_drsc_model_random_intercept_2011_2023)




## Linear random effects model DRSC 2011-2023 -----------------------------------------


summarytable_linear_drsc_model_random_intercept_slope_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               twoclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               threeclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               fourclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               fiveclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               sixclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               sevenclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                                               which=c("G", 
                                                                                                       "loglik", 
                                                                                                       "conv", 
                                                                                                       "npm", 
                                                                                                       "AIC", 
                                                                                                       "BIC", 
                                                                                                       "SABIC", 
                                                                                                       "entropy", 
                                                                                                       "ICL", 
                                                                                                       "ICL1", 
                                                                                                       "ICL2", 
                                                                                                       "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob(threeclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob(fourclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob(fiveclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob(sixclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob(sevenclass_linear_drsc_model_random_intercept_slope_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(twoclass_linear_drsc_model_random_intercept_slope_2011_2023)  
postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(threeclass_linear_drsc_model_random_intercept_slope_2011_2023)  
postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(fourclass_linear_drsc_model_random_intercept_slope_2011_2023)  
postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(fiveclass_linear_drsc_model_random_intercept_slope_2011_2023)  
postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(sixclass_linear_drsc_model_random_intercept_slope_2011_2023)
postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2023 <- postprob(sevenclass_linear_drsc_model_random_intercept_slope_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_drsc_model_random_intercept_slope_2011_2023 <- list(twoclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                threeclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                fourclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                fiveclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                sixclass_linear_drsc_model_random_intercept_slope_2011_2023,
                                                                sevenclass_linear_drsc_model_random_intercept_slope_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2023 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_drsc_model_random_intercept_slope_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$loglik[7], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$npm[7], summarytable_linear_drsc_model_random_intercept_slope_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_drsc_model_random_intercept_slope_2011_2023 <- sapply(outputs_vllrt_linear_drsc_model_random_intercept_slope_2011_2023, function(x) gsub("^= ", "", x))



summarytable_linear_drsc_model_random_intercept_slope_2011_2023 <- summarytable_linear_drsc_model_random_intercept_slope_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,]), min(postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,]), min(postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,]), min(postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_linear_drsc_model_random_intercept_slope_2011_2023$ns, min(postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]), min(postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]), min(postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]), min(postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_drsc_model_random_intercept_slope_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_drsc_model_random_intercept_slope_2011_2023, "summarytable_linear_drsc_model_random_intercept_slope_2011_2023.csv", row.names = F)

View(summarytable_linear_drsc_model_random_intercept_slope_2011_2023)



## Linear random effects model DRSC 2011-2023 -----------------------------------------


summarytable_quadratic_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_quadratic_drsc_model_2011_2023,
                                                                                                     twoclass_quadratic_drsc_model_2011_2023,
                                                                                                     threeclass_quadratic_drsc_model_2011_2023,
                                                                                                     fourclass_quadratic_drsc_model_2011_2023,
                                                                                                     fiveclass_quadratic_drsc_model_2011_2023,
                                                                                                     sixclass_quadratic_drsc_model_2011_2023,
                                                                                                     sevenclass_quadratic_drsc_model_2011_2023,
                                                                                                     which=c("G", 
                                                                                                             "loglik", 
                                                                                                             "conv", 
                                                                                                             "npm", 
                                                                                                             "AIC", 
                                                                                                             "BIC", 
                                                                                                             "SABIC", 
                                                                                                             "entropy", 
                                                                                                             "ICL", 
                                                                                                             "ICL1", 
                                                                                                             "ICL2", 
                                                                                                             "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_quadratic_drsc_model_2011_2023)
postprob(threeclass_quadratic_drsc_model_2011_2023)
postprob(fourclass_quadratic_drsc_model_2011_2023)
postprob(fiveclass_quadratic_drsc_model_2011_2023)
postprob(sixclass_quadratic_drsc_model_2011_2023)
postprob(sevenclass_quadratic_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_quadratic_drsc_model_2011_2023 <- postprob(twoclass_quadratic_drsc_model_2011_2023)  
postprob_threeclass_quadratic_drsc_model_2011_2023 <- postprob(threeclass_quadratic_drsc_model_2011_2023)  
postprob_fourclass_quadratic_drsc_model_2011_2023 <- postprob(fourclass_quadratic_drsc_model_2011_2023)  
postprob_fiveclass_quadratic_drsc_model_2011_2023 <- postprob(fiveclass_quadratic_drsc_model_2011_2023)  
postprob_sixclass_quadratic_drsc_model_2011_2023 <- postprob(sixclass_quadratic_drsc_model_2011_2023)
postprob_sevenclass_quadratic_drsc_model_2011_2023 <- postprob(sevenclass_quadratic_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_quadratic_drsc_model_2011_2023 <- list(twoclass_quadratic_drsc_model_2011_2023,
                                                                      threeclass_quadratic_drsc_model_2011_2023,
                                                                      fourclass_quadratic_drsc_model_2011_2023,
                                                                      fiveclass_quadratic_drsc_model_2011_2023,
                                                                      sixclass_quadratic_drsc_model_2011_2023,
                                                                      sevenclass_quadratic_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_quadratic_drsc_model_2011_2023 <- sapply(model_list_quadratic_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_quadratic_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_quadratic_drsc_model_2011_2023 <- sapply(model_list_quadratic_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_quadratic_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_quadratic_drsc_model_2011_2023 <- sapply(model_list_quadratic_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_quadratic_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_quadratic_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[1], summarytable_quadratic_drsc_model_2011_2023$npm[1], summarytable_quadratic_drsc_model_2011_2023$G[1], summarytable_quadratic_drsc_model_2011_2023$loglik[2], summarytable_quadratic_drsc_model_2011_2023$npm[2], summarytable_quadratic_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[2], summarytable_quadratic_drsc_model_2011_2023$npm[2], summarytable_quadratic_drsc_model_2011_2023$G[2], summarytable_quadratic_drsc_model_2011_2023$loglik[3], summarytable_quadratic_drsc_model_2011_2023$npm[3], summarytable_quadratic_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[3], summarytable_quadratic_drsc_model_2011_2023$npm[3], summarytable_quadratic_drsc_model_2011_2023$G[3], summarytable_quadratic_drsc_model_2011_2023$loglik[4], summarytable_quadratic_drsc_model_2011_2023$npm[4], summarytable_quadratic_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[4], summarytable_quadratic_drsc_model_2011_2023$npm[4], summarytable_quadratic_drsc_model_2011_2023$G[4], summarytable_quadratic_drsc_model_2011_2023$loglik[5], summarytable_quadratic_drsc_model_2011_2023$npm[5], summarytable_quadratic_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[5], summarytable_quadratic_drsc_model_2011_2023$npm[5], summarytable_quadratic_drsc_model_2011_2023$G[5], summarytable_quadratic_drsc_model_2011_2023$loglik[6], summarytable_quadratic_drsc_model_2011_2023$npm[6], summarytable_quadratic_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_quadratic_drsc_model_2011_2023$ns, summarytable_quadratic_drsc_model_2011_2023$loglik[6], summarytable_quadratic_drsc_model_2011_2023$npm[6], summarytable_quadratic_drsc_model_2011_2023$G[6], summarytable_quadratic_drsc_model_2011_2023$loglik[7], summarytable_quadratic_drsc_model_2011_2023$npm[7], summarytable_quadratic_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_quadratic_drsc_model_2011_2023 <- sapply(outputs_vllrt_quadratic_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_quadratic_drsc_model_2011_2023 <- summarytable_quadratic_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_quadratic_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_quadratic_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_quadratic_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_quadratic_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_quadratic_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_quadratic_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_quadratic_drsc_model_2011_2023$ns, min(postprob_twoclass_quadratic_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_quadratic_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_quadratic_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_quadratic_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_quadratic_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_quadratic_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_quadratic_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_quadratic_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_quadratic_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_quadratic_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_quadratic_drsc_model_2011_2023, "summarytable_quadratic_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_quadratic_drsc_model_2011_2023)





## Quadratics effects prop model DRSC 2011-2023 -----------------------------------------


summarytable_quadratic_prop_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 twoclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 threeclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 fourclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 fiveclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 sixclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 sevenclass_quadratic_prop_drsc_model_2011_2023,
                                                                                 which=c("G", 
                                                                                         "loglik", 
                                                                                         "conv", 
                                                                                         "npm", 
                                                                                         "AIC", 
                                                                                         "BIC", 
                                                                                         "SABIC", 
                                                                                         "entropy", 
                                                                                         "ICL", 
                                                                                         "ICL1", 
                                                                                         "ICL2", 
                                                                                         "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_quadratic_prop_drsc_model_2011_2023)
postprob(threeclass_quadratic_prop_drsc_model_2011_2023)
postprob(fourclass_quadratic_prop_drsc_model_2011_2023)
postprob(fiveclass_quadratic_prop_drsc_model_2011_2023)
postprob(sixclass_quadratic_prop_drsc_model_2011_2023)
postprob(sevenclass_quadratic_prop_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_quadratic_prop_drsc_model_2011_2023 <- postprob(twoclass_quadratic_prop_drsc_model_2011_2023)  
postprob_threeclass_quadratic_prop_drsc_model_2011_2023 <- postprob(threeclass_quadratic_prop_drsc_model_2011_2023)  
postprob_fourclass_quadratic_prop_drsc_model_2011_2023 <- postprob(fourclass_quadratic_prop_drsc_model_2011_2023)  
postprob_fiveclass_quadratic_prop_drsc_model_2011_2023 <- postprob(fiveclass_quadratic_prop_drsc_model_2011_2023)  
postprob_sixclass_quadratic_prop_drsc_model_2011_2023 <- postprob(sixclass_quadratic_prop_drsc_model_2011_2023)
postprob_sevenclass_quadratic_prop_drsc_model_2011_2023 <- postprob(sevenclass_quadratic_prop_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_quadratic_prop_drsc_model_2011_2023 <- list(twoclass_quadratic_prop_drsc_model_2011_2023,
                                                  threeclass_quadratic_prop_drsc_model_2011_2023,
                                                  fourclass_quadratic_prop_drsc_model_2011_2023,
                                                  fiveclass_quadratic_prop_drsc_model_2011_2023,
                                                  sixclass_quadratic_prop_drsc_model_2011_2023,
                                                  sevenclass_quadratic_prop_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_quadratic_prop_drsc_model_2011_2023 <- sapply(model_list_quadratic_prop_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_quadratic_prop_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_quadratic_prop_drsc_model_2011_2023 <- sapply(model_list_quadratic_prop_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_quadratic_prop_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_quadratic_prop_drsc_model_2011_2023 <- sapply(model_list_quadratic_prop_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_quadratic_prop_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_quadratic_prop_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[1], summarytable_quadratic_prop_drsc_model_2011_2023$npm[1], summarytable_quadratic_prop_drsc_model_2011_2023$G[1], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[2], summarytable_quadratic_prop_drsc_model_2011_2023$npm[2], summarytable_quadratic_prop_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[2], summarytable_quadratic_prop_drsc_model_2011_2023$npm[2], summarytable_quadratic_prop_drsc_model_2011_2023$G[2], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[3], summarytable_quadratic_prop_drsc_model_2011_2023$npm[3], summarytable_quadratic_prop_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[3], summarytable_quadratic_prop_drsc_model_2011_2023$npm[3], summarytable_quadratic_prop_drsc_model_2011_2023$G[3], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[4], summarytable_quadratic_prop_drsc_model_2011_2023$npm[4], summarytable_quadratic_prop_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[4], summarytable_quadratic_prop_drsc_model_2011_2023$npm[4], summarytable_quadratic_prop_drsc_model_2011_2023$G[4], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[5], summarytable_quadratic_prop_drsc_model_2011_2023$npm[5], summarytable_quadratic_prop_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[5], summarytable_quadratic_prop_drsc_model_2011_2023$npm[5], summarytable_quadratic_prop_drsc_model_2011_2023$G[5], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[6], summarytable_quadratic_prop_drsc_model_2011_2023$npm[6], summarytable_quadratic_prop_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_quadratic_prop_drsc_model_2011_2023$ns, summarytable_quadratic_prop_drsc_model_2011_2023$loglik[6], summarytable_quadratic_prop_drsc_model_2011_2023$npm[6], summarytable_quadratic_prop_drsc_model_2011_2023$G[6], summarytable_quadratic_prop_drsc_model_2011_2023$loglik[7], summarytable_quadratic_prop_drsc_model_2011_2023$npm[7], summarytable_quadratic_prop_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_quadratic_prop_drsc_model_2011_2023 <- sapply(outputs_vllrt_quadratic_prop_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_quadratic_prop_drsc_model_2011_2023 <- summarytable_quadratic_prop_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_quadratic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_quadratic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_quadratic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_quadratic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_quadratic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_quadratic_prop_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_quadratic_prop_drsc_model_2011_2023$ns, min(postprob_twoclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_quadratic_prop_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_quadratic_prop_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_quadratic_prop_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_quadratic_prop_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_quadratic_prop_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_quadratic_prop_drsc_model_2011_2023, "summarytable_quadratic_prop_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_quadratic_prop_drsc_model_2011_2023)



# Cubic effects prop model DRSC 2011-2023 -----------------------------------------


summarytable_cubic_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_cubic_drsc_model_2011_2023,
                                                                                      twoclass_cubic_drsc_model_2011_2023,
                                                                                      threeclass_cubic_drsc_model_2011_2023,
                                                                                      fourclass_cubic_drsc_model_2011_2023,
                                                                                      fiveclass_cubic_drsc_model_2011_2023,
                                                                                      sixclass_cubic_drsc_model_2011_2023,
                                                                                      sevenclass_cubic_drsc_model_2011_2023,
                                                                                      which=c("G", 
                                                                                              "loglik", 
                                                                                              "conv", 
                                                                                              "npm", 
                                                                                              "AIC", 
                                                                                              "BIC", 
                                                                                              "SABIC", 
                                                                                              "entropy", 
                                                                                              "ICL", 
                                                                                              "ICL1", 
                                                                                              "ICL2", 
                                                                                              "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_cubic_drsc_model_2011_2023)
postprob(threeclass_cubic_drsc_model_2011_2023)
postprob(fourclass_cubic_drsc_model_2011_2023)
postprob(fiveclass_cubic_drsc_model_2011_2023)
postprob(sixclass_cubic_drsc_model_2011_2023)
postprob(sevenclass_cubic_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_cubic_drsc_model_2011_2023 <- postprob(twoclass_cubic_drsc_model_2011_2023)  
postprob_threeclass_cubic_drsc_model_2011_2023 <- postprob(threeclass_cubic_drsc_model_2011_2023)  
postprob_fourclass_cubic_drsc_model_2011_2023 <- postprob(fourclass_cubic_drsc_model_2011_2023)  
postprob_fiveclass_cubic_drsc_model_2011_2023 <- postprob(fiveclass_cubic_drsc_model_2011_2023)  
postprob_sixclass_cubic_drsc_model_2011_2023 <- postprob(sixclass_cubic_drsc_model_2011_2023)
postprob_sevenclass_cubic_drsc_model_2011_2023 <- postprob(sevenclass_cubic_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_cubic_drsc_model_2011_2023 <- list(twoclass_cubic_drsc_model_2011_2023,
                                                       threeclass_cubic_drsc_model_2011_2023,
                                                       fourclass_cubic_drsc_model_2011_2023,
                                                       fiveclass_cubic_drsc_model_2011_2023,
                                                       sixclass_cubic_drsc_model_2011_2023,
                                                       sevenclass_cubic_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_cubic_drsc_model_2011_2023 <- sapply(model_list_cubic_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_cubic_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_cubic_drsc_model_2011_2023 <- sapply(model_list_cubic_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_cubic_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_cubic_drsc_model_2011_2023 <- sapply(model_list_cubic_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_cubic_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_cubic_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[1], summarytable_cubic_drsc_model_2011_2023$npm[1], summarytable_cubic_drsc_model_2011_2023$G[1], summarytable_cubic_drsc_model_2011_2023$loglik[2], summarytable_cubic_drsc_model_2011_2023$npm[2], summarytable_cubic_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[2], summarytable_cubic_drsc_model_2011_2023$npm[2], summarytable_cubic_drsc_model_2011_2023$G[2], summarytable_cubic_drsc_model_2011_2023$loglik[3], summarytable_cubic_drsc_model_2011_2023$npm[3], summarytable_cubic_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[3], summarytable_cubic_drsc_model_2011_2023$npm[3], summarytable_cubic_drsc_model_2011_2023$G[3], summarytable_cubic_drsc_model_2011_2023$loglik[4], summarytable_cubic_drsc_model_2011_2023$npm[4], summarytable_cubic_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[4], summarytable_cubic_drsc_model_2011_2023$npm[4], summarytable_cubic_drsc_model_2011_2023$G[4], summarytable_cubic_drsc_model_2011_2023$loglik[5], summarytable_cubic_drsc_model_2011_2023$npm[5], summarytable_cubic_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[5], summarytable_cubic_drsc_model_2011_2023$npm[5], summarytable_cubic_drsc_model_2011_2023$G[5], summarytable_cubic_drsc_model_2011_2023$loglik[6], summarytable_cubic_drsc_model_2011_2023$npm[6], summarytable_cubic_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_cubic_drsc_model_2011_2023$ns, summarytable_cubic_drsc_model_2011_2023$loglik[6], summarytable_cubic_drsc_model_2011_2023$npm[6], summarytable_cubic_drsc_model_2011_2023$G[6], summarytable_cubic_drsc_model_2011_2023$loglik[7], summarytable_cubic_drsc_model_2011_2023$npm[7], summarytable_cubic_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_cubic_drsc_model_2011_2023 <- sapply(outputs_vllrt_cubic_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_cubic_drsc_model_2011_2023 <- summarytable_cubic_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_cubic_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_cubic_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_cubic_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_cubic_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_cubic_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_cubic_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_cubic_drsc_model_2011_2023$ns, min(postprob_twoclass_cubic_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_cubic_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_cubic_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_cubic_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_cubic_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_cubic_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_cubic_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_cubic_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_cubic_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_cubic_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_cubic_drsc_model_2011_2023, "summarytable_cubic_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_cubic_drsc_model_2011_2023)




# cubic prop effects prop model DRSC 2011-2023 -----------------------------------------


summarytable_cubic_prop_drsc_model_2011_2023 <-  as.data.frame(lcmm::summarytable(oneclass_cubic_prop_drsc_model_2011_2023,
                                                                             twoclass_cubic_prop_drsc_model_2011_2023,
                                                                             threeclass_cubic_prop_drsc_model_2011_2023,
                                                                             fourclass_cubic_prop_drsc_model_2011_2023,
                                                                             fiveclass_cubic_prop_drsc_model_2011_2023,
                                                                             sixclass_cubic_prop_drsc_model_2011_2023,
                                                                             sevenclass_cubic_prop_drsc_model_2011_2023,
                                                                             which=c("G", 
                                                                                     "loglik", 
                                                                                     "conv", 
                                                                                     "npm", 
                                                                                     "AIC", 
                                                                                     "BIC", 
                                                                                     "SABIC", 
                                                                                     "entropy", 
                                                                                     "ICL", 
                                                                                     "ICL1", 
                                                                                     "ICL2", 
                                                                                     "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_cubic_prop_drsc_model_2011_2023)
postprob(threeclass_cubic_prop_drsc_model_2011_2023)
postprob(fourclass_cubic_prop_drsc_model_2011_2023)
postprob(fiveclass_cubic_prop_drsc_model_2011_2023)
postprob(sixclass_cubic_prop_drsc_model_2011_2023)
postprob(sevenclass_cubic_prop_drsc_model_2011_2023)


# Assuming postprob() returns a structured list
postprob_twoclass_cubic_prop_drsc_model_2011_2023 <- postprob(twoclass_cubic_prop_drsc_model_2011_2023)  
postprob_threeclass_cubic_prop_drsc_model_2011_2023 <- postprob(threeclass_cubic_prop_drsc_model_2011_2023)  
postprob_fourclass_cubic_prop_drsc_model_2011_2023 <- postprob(fourclass_cubic_prop_drsc_model_2011_2023)  
postprob_fiveclass_cubic_prop_drsc_model_2011_2023 <- postprob(fiveclass_cubic_prop_drsc_model_2011_2023)  
postprob_sixclass_cubic_prop_drsc_model_2011_2023 <- postprob(sixclass_cubic_prop_drsc_model_2011_2023)
postprob_sevenclass_cubic_prop_drsc_model_2011_2023 <- postprob(sevenclass_cubic_prop_drsc_model_2011_2023)  



# Assuming you have a list of your model objects --------------------------


model_list_cubic_prop_drsc_model_2011_2023 <- list(twoclass_cubic_prop_drsc_model_2011_2023,
                                              threeclass_cubic_prop_drsc_model_2011_2023,
                                              fourclass_cubic_prop_drsc_model_2011_2023,
                                              fiveclass_cubic_prop_drsc_model_2011_2023,
                                              sixclass_cubic_prop_drsc_model_2011_2023,
                                              sevenclass_cubic_prop_drsc_model_2011_2023)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_cubic_prop_drsc_model_2011_2023 <- sapply(model_list_cubic_prop_drsc_model_2011_2023, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_cubic_prop_drsc_model_2011_2023)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_cubic_prop_drsc_model_2011_2023 <- sapply(model_list_cubic_prop_drsc_model_2011_2023, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_cubic_prop_drsc_model_2011_2023)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_cubic_prop_drsc_model_2011_2023 <- sapply(model_list_cubic_prop_drsc_model_2011_2023, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_cubic_prop_drsc_model_2011_2023)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_cubic_prop_drsc_model_2011_2023 <- list(
  extract_p_value_from_lrt(oneclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[1], summarytable_cubic_prop_drsc_model_2011_2023$npm[1], summarytable_cubic_prop_drsc_model_2011_2023$G[1], summarytable_cubic_prop_drsc_model_2011_2023$loglik[2], summarytable_cubic_prop_drsc_model_2011_2023$npm[2], summarytable_cubic_prop_drsc_model_2011_2023$G[2]),
  extract_p_value_from_lrt(twoclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[2], summarytable_cubic_prop_drsc_model_2011_2023$npm[2], summarytable_cubic_prop_drsc_model_2011_2023$G[2], summarytable_cubic_prop_drsc_model_2011_2023$loglik[3], summarytable_cubic_prop_drsc_model_2011_2023$npm[3], summarytable_cubic_prop_drsc_model_2011_2023$G[3]),
  extract_p_value_from_lrt(threeclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[3], summarytable_cubic_prop_drsc_model_2011_2023$npm[3], summarytable_cubic_prop_drsc_model_2011_2023$G[3], summarytable_cubic_prop_drsc_model_2011_2023$loglik[4], summarytable_cubic_prop_drsc_model_2011_2023$npm[4], summarytable_cubic_prop_drsc_model_2011_2023$G[4]),
  extract_p_value_from_lrt(fourclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[4], summarytable_cubic_prop_drsc_model_2011_2023$npm[4], summarytable_cubic_prop_drsc_model_2011_2023$G[4], summarytable_cubic_prop_drsc_model_2011_2023$loglik[5], summarytable_cubic_prop_drsc_model_2011_2023$npm[5], summarytable_cubic_prop_drsc_model_2011_2023$G[5]),
  extract_p_value_from_lrt(fiveclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[5], summarytable_cubic_prop_drsc_model_2011_2023$npm[5], summarytable_cubic_prop_drsc_model_2011_2023$G[5], summarytable_cubic_prop_drsc_model_2011_2023$loglik[6], summarytable_cubic_prop_drsc_model_2011_2023$npm[6], summarytable_cubic_prop_drsc_model_2011_2023$G[6]),
  extract_p_value_from_lrt(sixclass_cubic_prop_drsc_model_2011_2023$ns, summarytable_cubic_prop_drsc_model_2011_2023$loglik[6], summarytable_cubic_prop_drsc_model_2011_2023$npm[6], summarytable_cubic_prop_drsc_model_2011_2023$G[6], summarytable_cubic_prop_drsc_model_2011_2023$loglik[7], summarytable_cubic_prop_drsc_model_2011_2023$npm[7], summarytable_cubic_prop_drsc_model_2011_2023$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_cubic_prop_drsc_model_2011_2023 <- sapply(outputs_vllrt_cubic_prop_drsc_model_2011_2023, function(x) gsub("^= ", "", x))



summarytable_cubic_prop_drsc_model_2011_2023 <- summarytable_cubic_prop_drsc_model_2011_2023 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_cubic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_threeclass_cubic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_fourclass_cubic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_fiveclass_cubic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_sixclass_cubic_prop_drsc_model_2011_2023[[1]][2,]), min(postprob_sevenclass_cubic_prop_drsc_model_2011_2023[[1]][2,])),
         smallest_class_count = c(oneclass_cubic_prop_drsc_model_2011_2023$ns, min(postprob_twoclass_cubic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_threeclass_cubic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_fourclass_cubic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_fiveclass_cubic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_sixclass_cubic_prop_drsc_model_2011_2023[[1]][1,]), min(postprob_sevenclass_cubic_prop_drsc_model_2011_2023[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_cubic_prop_drsc_model_2011_2023),
         "Highest MMV" =c(NA, highest_mismatch_values_cubic_prop_drsc_model_2011_2023),
         "Lowest OCC" = c(NA, lower_occ_values_cubic_prop_drsc_model_2011_2023),
         VLMRLRT = c(NA, values_vllrt_outputs_cubic_prop_drsc_model_2011_2023)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_cubic_prop_drsc_model_2011_2023, "summarytable_cubic_prop_drsc_model_2011_2023.csv", row.names = F)

View(summarytable_cubic_prop_drsc_model_2011_2023)








# Step 2 - run a random effects model --------------------------------------------------

## Linear non random effects model DRSC 2011-2019 -----------------------------------------


summarytable_linear_nre_homocedastic_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               twoclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               threeclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               fourclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               fiveclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               sixclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               sevenclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                                               which=c("G", 
                                                                                                       "loglik", 
                                                                                                       "conv", 
                                                                                                       "npm", 
                                                                                                       "AIC", 
                                                                                                       "BIC", 
                                                                                                       "SABIC", 
                                                                                                       "entropy", 
                                                                                                       "ICL", 
                                                                                                       "ICL1", 
                                                                                                       "ICL2", 
                                                                                                       "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob(threeclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob(fourclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob(fiveclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob(sixclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob(sevenclass_linear_nre_homocedastic_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(twoclass_linear_nre_homocedastic_drsc_model_2011_2019)  
postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(threeclass_linear_nre_homocedastic_drsc_model_2011_2019)  
postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(fourclass_linear_nre_homocedastic_drsc_model_2011_2019)  
postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(fiveclass_linear_nre_homocedastic_drsc_model_2011_2019)  
postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(sixclass_linear_nre_homocedastic_drsc_model_2011_2019)
postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2019 <- postprob(sevenclass_linear_nre_homocedastic_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_nre_homocedastic_drsc_model_2011_2019 <- list(twoclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                threeclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                fourclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                fiveclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                sixclass_linear_nre_homocedastic_drsc_model_2011_2019,
                                                                sevenclass_linear_nre_homocedastic_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_homocedastic_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_nre_homocedastic_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[1], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[2], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[3], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[4], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[5], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[6], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$loglik[7], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$npm[7], summarytable_linear_nre_homocedastic_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_nre_homocedastic_drsc_model_2011_2019 <- sapply(outputs_vllrt_linear_nre_homocedastic_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_linear_nre_homocedastic_drsc_model_2011_2019 <- summarytable_linear_nre_homocedastic_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_linear_nre_homocedastic_drsc_model_2011_2019$ns, min(postprob_twoclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_linear_nre_homocedastic_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_nre_homocedastic_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_nre_homocedastic_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_linear_nre_homocedastic_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_nre_homocedastic_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_nre_homocedastic_drsc_model_2011_2019, "summarytable_linear_nre_homocedastic_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_linear_nre_homocedastic_drsc_model_2011_2019)


## Linear non random effects model DRSC 2011-2019 -----------------------------------------


summarytable_linear_nre_heterocedastic_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 twoclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 threeclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 fourclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 sixclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                                                 which=c("G", 
                                                                                                         "loglik", 
                                                                                                         "conv", 
                                                                                                         "npm", 
                                                                                                         "AIC", 
                                                                                                         "BIC", 
                                                                                                         "SABIC", 
                                                                                                         "entropy", 
                                                                                                         "ICL", 
                                                                                                         "ICL1", 
                                                                                                         "ICL2", 
                                                                                                         "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob(threeclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob(fourclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob(sixclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob(sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(twoclass_linear_nre_heterocedastic_drsc_model_2011_2019)  
postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(threeclass_linear_nre_heterocedastic_drsc_model_2011_2019)  
postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(fourclass_linear_nre_heterocedastic_drsc_model_2011_2019)  
postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019)  
postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(sixclass_linear_nre_heterocedastic_drsc_model_2011_2019)
postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019 <- postprob(sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_nre_heterocedastic_drsc_model_2011_2019 <- list(twoclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                  threeclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                  fourclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                  fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                  sixclass_linear_nre_heterocedastic_drsc_model_2011_2019,
                                                                  sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2019 <- sapply(model_list_linear_nre_heterocedastic_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_nre_heterocedastic_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[1], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[2], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[3], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[4], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[5], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[6], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$loglik[7], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$npm[7], summarytable_linear_nre_heterocedastic_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_nre_heterocedastic_drsc_model_2011_2019 <- sapply(outputs_vllrt_linear_nre_heterocedastic_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_linear_nre_heterocedastic_drsc_model_2011_2019 <- summarytable_linear_nre_heterocedastic_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_linear_nre_heterocedastic_drsc_model_2011_2019$ns, min(postprob_twoclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_linear_nre_heterocedastic_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_nre_heterocedastic_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_nre_heterocedastic_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_linear_nre_heterocedastic_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_nre_heterocedastic_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_nre_heterocedastic_drsc_model_2011_2019, "summarytable_linear_nre_heterocedastic_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_linear_nre_heterocedastic_drsc_model_2011_2019)



## Linear random effects model DRSC 2011-2019 -----------------------------------------


summarytable_linear_drsc_model_random_intercept_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               twoclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               threeclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               fourclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               fiveclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               sixclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               sevenclass_linear_drsc_model_random_intercept_2011_2019,
                                                                                               which=c("G", 
                                                                                                       "loglik", 
                                                                                                       "conv", 
                                                                                                       "npm", 
                                                                                                       "AIC", 
                                                                                                       "BIC", 
                                                                                                       "SABIC", 
                                                                                                       "entropy", 
                                                                                                       "ICL", 
                                                                                                       "ICL1", 
                                                                                                       "ICL2", 
                                                                                                       "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_drsc_model_random_intercept_2011_2019)
postprob(threeclass_linear_drsc_model_random_intercept_2011_2019)
postprob(fourclass_linear_drsc_model_random_intercept_2011_2019)
postprob(fiveclass_linear_drsc_model_random_intercept_2011_2019)
postprob(sixclass_linear_drsc_model_random_intercept_2011_2019)
postprob(sevenclass_linear_drsc_model_random_intercept_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(twoclass_linear_drsc_model_random_intercept_2011_2019)  
postprob_threeclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(threeclass_linear_drsc_model_random_intercept_2011_2019)  
postprob_fourclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(fourclass_linear_drsc_model_random_intercept_2011_2019)  
postprob_fiveclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(fiveclass_linear_drsc_model_random_intercept_2011_2019)  
postprob_sixclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(sixclass_linear_drsc_model_random_intercept_2011_2019)
postprob_sevenclass_linear_drsc_model_random_intercept_2011_2019 <- postprob(sevenclass_linear_drsc_model_random_intercept_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_drsc_model_random_intercept_2011_2019 <- list(twoclass_linear_drsc_model_random_intercept_2011_2019,
                                                                threeclass_linear_drsc_model_random_intercept_2011_2019,
                                                                fourclass_linear_drsc_model_random_intercept_2011_2019,
                                                                fiveclass_linear_drsc_model_random_intercept_2011_2019,
                                                                sixclass_linear_drsc_model_random_intercept_2011_2019,
                                                                sevenclass_linear_drsc_model_random_intercept_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_drsc_model_random_intercept_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_drsc_model_random_intercept_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_drsc_model_random_intercept_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_drsc_model_random_intercept_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_drsc_model_random_intercept_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_drsc_model_random_intercept_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_drsc_model_random_intercept_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[1], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[1], summarytable_linear_drsc_model_random_intercept_2011_2019$G[1], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[2], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[2], summarytable_linear_drsc_model_random_intercept_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[2], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[2], summarytable_linear_drsc_model_random_intercept_2011_2019$G[2], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[3], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[3], summarytable_linear_drsc_model_random_intercept_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[3], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[3], summarytable_linear_drsc_model_random_intercept_2011_2019$G[3], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[4], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[4], summarytable_linear_drsc_model_random_intercept_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[4], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[4], summarytable_linear_drsc_model_random_intercept_2011_2019$G[4], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[5], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[5], summarytable_linear_drsc_model_random_intercept_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[5], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[5], summarytable_linear_drsc_model_random_intercept_2011_2019$G[5], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[6], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[6], summarytable_linear_drsc_model_random_intercept_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_linear_drsc_model_random_intercept_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[6], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[6], summarytable_linear_drsc_model_random_intercept_2011_2019$G[6], summarytable_linear_drsc_model_random_intercept_2011_2019$loglik[7], summarytable_linear_drsc_model_random_intercept_2011_2019$npm[7], summarytable_linear_drsc_model_random_intercept_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_drsc_model_random_intercept_2011_2019 <- sapply(outputs_vllrt_linear_drsc_model_random_intercept_2011_2019, function(x) gsub("^= ", "", x))



summarytable_linear_drsc_model_random_intercept_2011_2019 <- summarytable_linear_drsc_model_random_intercept_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,]), min(postprob_threeclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,]), min(postprob_fourclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,]), min(postprob_sixclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_linear_drsc_model_random_intercept_2011_2019$ns, min(postprob_twoclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]), min(postprob_threeclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]), min(postprob_fourclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]), min(postprob_sixclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_drsc_model_random_intercept_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_drsc_model_random_intercept_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_linear_drsc_model_random_intercept_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_drsc_model_random_intercept_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_drsc_model_random_intercept_2011_2019, "summarytable_linear_drsc_model_random_intercept_2011_2019.csv", row.names = F)

View(summarytable_linear_drsc_model_random_intercept_2011_2019)




## Linear random effects model DRSC 2011-2019 -----------------------------------------


summarytable_linear_drsc_model_random_intercept_slope_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     twoclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     threeclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     fourclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     fiveclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     sixclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     sevenclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                                                     which=c("G", 
                                                                                                             "loglik", 
                                                                                                             "conv", 
                                                                                                             "npm", 
                                                                                                             "AIC", 
                                                                                                             "BIC", 
                                                                                                             "SABIC", 
                                                                                                             "entropy", 
                                                                                                             "ICL", 
                                                                                                             "ICL1", 
                                                                                                             "ICL2", 
                                                                                                             "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob(threeclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob(fourclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob(fiveclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob(sixclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob(sevenclass_linear_drsc_model_random_intercept_slope_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(twoclass_linear_drsc_model_random_intercept_slope_2011_2019)  
postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(threeclass_linear_drsc_model_random_intercept_slope_2011_2019)  
postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(fourclass_linear_drsc_model_random_intercept_slope_2011_2019)  
postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(fiveclass_linear_drsc_model_random_intercept_slope_2011_2019)  
postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(sixclass_linear_drsc_model_random_intercept_slope_2011_2019)
postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2019 <- postprob(sevenclass_linear_drsc_model_random_intercept_slope_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_linear_drsc_model_random_intercept_slope_2011_2019 <- list(twoclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                      threeclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                      fourclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                      fiveclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                      sixclass_linear_drsc_model_random_intercept_slope_2011_2019,
                                                                      sevenclass_linear_drsc_model_random_intercept_slope_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2019 <- sapply(model_list_linear_drsc_model_random_intercept_slope_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_linear_drsc_model_random_intercept_slope_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[1], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[2], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[3], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[4], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[5], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[6], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$loglik[7], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$npm[7], summarytable_linear_drsc_model_random_intercept_slope_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_linear_drsc_model_random_intercept_slope_2011_2019 <- sapply(outputs_vllrt_linear_drsc_model_random_intercept_slope_2011_2019, function(x) gsub("^= ", "", x))



summarytable_linear_drsc_model_random_intercept_slope_2011_2019 <- summarytable_linear_drsc_model_random_intercept_slope_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,]), min(postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,]), min(postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,]), min(postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_linear_drsc_model_random_intercept_slope_2011_2019$ns, min(postprob_twoclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]), min(postprob_threeclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]), min(postprob_fourclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]), min(postprob_fiveclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]), min(postprob_sixclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]), min(postprob_sevenclass_linear_drsc_model_random_intercept_slope_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_linear_drsc_model_random_intercept_slope_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_linear_drsc_model_random_intercept_slope_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_linear_drsc_model_random_intercept_slope_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_linear_drsc_model_random_intercept_slope_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_linear_drsc_model_random_intercept_slope_2011_2019, "summarytable_linear_drsc_model_random_intercept_slope_2011_2019.csv", row.names = F)

View(summarytable_linear_drsc_model_random_intercept_slope_2011_2019)



## Linear random effects model DRSC 2011-2019 -----------------------------------------


summarytable_quadratic_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_quadratic_drsc_model_2011_2019,
                                                                                 twoclass_quadratic_drsc_model_2011_2019,
                                                                                 threeclass_quadratic_drsc_model_2011_2019,
                                                                                 fourclass_quadratic_drsc_model_2011_2019,
                                                                                 fiveclass_quadratic_drsc_model_2011_2019,
                                                                                 sixclass_quadratic_drsc_model_2011_2019,
                                                                                 sevenclass_quadratic_drsc_model_2011_2019,
                                                                                 which=c("G", 
                                                                                         "loglik", 
                                                                                         "conv", 
                                                                                         "npm", 
                                                                                         "AIC", 
                                                                                         "BIC", 
                                                                                         "SABIC", 
                                                                                         "entropy", 
                                                                                         "ICL", 
                                                                                         "ICL1", 
                                                                                         "ICL2", 
                                                                                         "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_quadratic_drsc_model_2011_2019)
postprob(threeclass_quadratic_drsc_model_2011_2019)
postprob(fourclass_quadratic_drsc_model_2011_2019)
postprob(fiveclass_quadratic_drsc_model_2011_2019)
postprob(sixclass_quadratic_drsc_model_2011_2019)
postprob(sevenclass_quadratic_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_quadratic_drsc_model_2011_2019 <- postprob(twoclass_quadratic_drsc_model_2011_2019)  
postprob_threeclass_quadratic_drsc_model_2011_2019 <- postprob(threeclass_quadratic_drsc_model_2011_2019)  
postprob_fourclass_quadratic_drsc_model_2011_2019 <- postprob(fourclass_quadratic_drsc_model_2011_2019)  
postprob_fiveclass_quadratic_drsc_model_2011_2019 <- postprob(fiveclass_quadratic_drsc_model_2011_2019)  
postprob_sixclass_quadratic_drsc_model_2011_2019 <- postprob(sixclass_quadratic_drsc_model_2011_2019)
postprob_sevenclass_quadratic_drsc_model_2011_2019 <- postprob(sevenclass_quadratic_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_quadratic_drsc_model_2011_2019 <- list(twoclass_quadratic_drsc_model_2011_2019,
                                                  threeclass_quadratic_drsc_model_2011_2019,
                                                  fourclass_quadratic_drsc_model_2011_2019,
                                                  fiveclass_quadratic_drsc_model_2011_2019,
                                                  sixclass_quadratic_drsc_model_2011_2019,
                                                  sevenclass_quadratic_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_quadratic_drsc_model_2011_2019 <- sapply(model_list_quadratic_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_quadratic_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_quadratic_drsc_model_2011_2019 <- sapply(model_list_quadratic_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_quadratic_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_quadratic_drsc_model_2011_2019 <- sapply(model_list_quadratic_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_quadratic_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_quadratic_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[1], summarytable_quadratic_drsc_model_2011_2019$npm[1], summarytable_quadratic_drsc_model_2011_2019$G[1], summarytable_quadratic_drsc_model_2011_2019$loglik[2], summarytable_quadratic_drsc_model_2011_2019$npm[2], summarytable_quadratic_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[2], summarytable_quadratic_drsc_model_2011_2019$npm[2], summarytable_quadratic_drsc_model_2011_2019$G[2], summarytable_quadratic_drsc_model_2011_2019$loglik[3], summarytable_quadratic_drsc_model_2011_2019$npm[3], summarytable_quadratic_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[3], summarytable_quadratic_drsc_model_2011_2019$npm[3], summarytable_quadratic_drsc_model_2011_2019$G[3], summarytable_quadratic_drsc_model_2011_2019$loglik[4], summarytable_quadratic_drsc_model_2011_2019$npm[4], summarytable_quadratic_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[4], summarytable_quadratic_drsc_model_2011_2019$npm[4], summarytable_quadratic_drsc_model_2011_2019$G[4], summarytable_quadratic_drsc_model_2011_2019$loglik[5], summarytable_quadratic_drsc_model_2011_2019$npm[5], summarytable_quadratic_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[5], summarytable_quadratic_drsc_model_2011_2019$npm[5], summarytable_quadratic_drsc_model_2011_2019$G[5], summarytable_quadratic_drsc_model_2011_2019$loglik[6], summarytable_quadratic_drsc_model_2011_2019$npm[6], summarytable_quadratic_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_quadratic_drsc_model_2011_2019$ns, summarytable_quadratic_drsc_model_2011_2019$loglik[6], summarytable_quadratic_drsc_model_2011_2019$npm[6], summarytable_quadratic_drsc_model_2011_2019$G[6], summarytable_quadratic_drsc_model_2011_2019$loglik[7], summarytable_quadratic_drsc_model_2011_2019$npm[7], summarytable_quadratic_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_quadratic_drsc_model_2011_2019 <- sapply(outputs_vllrt_quadratic_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_quadratic_drsc_model_2011_2019 <- summarytable_quadratic_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_quadratic_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_quadratic_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_quadratic_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_quadratic_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_quadratic_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_quadratic_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_quadratic_drsc_model_2011_2019$ns, min(postprob_twoclass_quadratic_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_quadratic_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_quadratic_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_quadratic_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_quadratic_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_quadratic_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_quadratic_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_quadratic_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_quadratic_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_quadratic_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_quadratic_drsc_model_2011_2019, "summarytable_quadratic_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_quadratic_drsc_model_2011_2019)





## Quadratics effects prop model DRSC 2011-2019 -----------------------------------------


summarytable_quadratic_prop_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      twoclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      threeclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      fourclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      fiveclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      sixclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      sevenclass_quadratic_prop_drsc_model_2011_2019,
                                                                                      which=c("G", 
                                                                                              "loglik", 
                                                                                              "conv", 
                                                                                              "npm", 
                                                                                              "AIC", 
                                                                                              "BIC", 
                                                                                              "SABIC", 
                                                                                              "entropy", 
                                                                                              "ICL", 
                                                                                              "ICL1", 
                                                                                              "ICL2", 
                                                                                              "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_quadratic_prop_drsc_model_2011_2019)
postprob(threeclass_quadratic_prop_drsc_model_2011_2019)
postprob(fourclass_quadratic_prop_drsc_model_2011_2019)
postprob(fiveclass_quadratic_prop_drsc_model_2011_2019)
postprob(sixclass_quadratic_prop_drsc_model_2011_2019)
postprob(sevenclass_quadratic_prop_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_quadratic_prop_drsc_model_2011_2019 <- postprob(twoclass_quadratic_prop_drsc_model_2011_2019)  
postprob_threeclass_quadratic_prop_drsc_model_2011_2019 <- postprob(threeclass_quadratic_prop_drsc_model_2011_2019)  
postprob_fourclass_quadratic_prop_drsc_model_2011_2019 <- postprob(fourclass_quadratic_prop_drsc_model_2011_2019)  
postprob_fiveclass_quadratic_prop_drsc_model_2011_2019 <- postprob(fiveclass_quadratic_prop_drsc_model_2011_2019)  
postprob_sixclass_quadratic_prop_drsc_model_2011_2019 <- postprob(sixclass_quadratic_prop_drsc_model_2011_2019)
postprob_sevenclass_quadratic_prop_drsc_model_2011_2019 <- postprob(sevenclass_quadratic_prop_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_quadratic_prop_drsc_model_2011_2019 <- list(twoclass_quadratic_prop_drsc_model_2011_2019,
                                                       threeclass_quadratic_prop_drsc_model_2011_2019,
                                                       fourclass_quadratic_prop_drsc_model_2011_2019,
                                                       fiveclass_quadratic_prop_drsc_model_2011_2019,
                                                       sixclass_quadratic_prop_drsc_model_2011_2019,
                                                       sevenclass_quadratic_prop_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_quadratic_prop_drsc_model_2011_2019 <- sapply(model_list_quadratic_prop_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_quadratic_prop_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_quadratic_prop_drsc_model_2011_2019 <- sapply(model_list_quadratic_prop_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_quadratic_prop_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_quadratic_prop_drsc_model_2011_2019 <- sapply(model_list_quadratic_prop_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_quadratic_prop_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_quadratic_prop_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[1], summarytable_quadratic_prop_drsc_model_2011_2019$npm[1], summarytable_quadratic_prop_drsc_model_2011_2019$G[1], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[2], summarytable_quadratic_prop_drsc_model_2011_2019$npm[2], summarytable_quadratic_prop_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[2], summarytable_quadratic_prop_drsc_model_2011_2019$npm[2], summarytable_quadratic_prop_drsc_model_2011_2019$G[2], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[3], summarytable_quadratic_prop_drsc_model_2011_2019$npm[3], summarytable_quadratic_prop_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[3], summarytable_quadratic_prop_drsc_model_2011_2019$npm[3], summarytable_quadratic_prop_drsc_model_2011_2019$G[3], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[4], summarytable_quadratic_prop_drsc_model_2011_2019$npm[4], summarytable_quadratic_prop_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[4], summarytable_quadratic_prop_drsc_model_2011_2019$npm[4], summarytable_quadratic_prop_drsc_model_2011_2019$G[4], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[5], summarytable_quadratic_prop_drsc_model_2011_2019$npm[5], summarytable_quadratic_prop_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[5], summarytable_quadratic_prop_drsc_model_2011_2019$npm[5], summarytable_quadratic_prop_drsc_model_2011_2019$G[5], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[6], summarytable_quadratic_prop_drsc_model_2011_2019$npm[6], summarytable_quadratic_prop_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_quadratic_prop_drsc_model_2011_2019$ns, summarytable_quadratic_prop_drsc_model_2011_2019$loglik[6], summarytable_quadratic_prop_drsc_model_2011_2019$npm[6], summarytable_quadratic_prop_drsc_model_2011_2019$G[6], summarytable_quadratic_prop_drsc_model_2011_2019$loglik[7], summarytable_quadratic_prop_drsc_model_2011_2019$npm[7], summarytable_quadratic_prop_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_quadratic_prop_drsc_model_2011_2019 <- sapply(outputs_vllrt_quadratic_prop_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_quadratic_prop_drsc_model_2011_2019 <- summarytable_quadratic_prop_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_quadratic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_quadratic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_quadratic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_quadratic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_quadratic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_quadratic_prop_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_quadratic_prop_drsc_model_2011_2019$ns, min(postprob_twoclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_quadratic_prop_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_quadratic_prop_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_quadratic_prop_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_quadratic_prop_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_quadratic_prop_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_quadratic_prop_drsc_model_2011_2019, "summarytable_quadratic_prop_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_quadratic_prop_drsc_model_2011_2019)



# Cubic effects prop model DRSC 2011-2019 -----------------------------------------


summarytable_cubic_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_cubic_drsc_model_2011_2019,
                                                                             twoclass_cubic_drsc_model_2011_2019,
                                                                             threeclass_cubic_drsc_model_2011_2019,
                                                                             fourclass_cubic_drsc_model_2011_2019,
                                                                             fiveclass_cubic_drsc_model_2011_2019,
                                                                             sixclass_cubic_drsc_model_2011_2019,
                                                                             sevenclass_cubic_drsc_model_2011_2019,
                                                                             which=c("G", 
                                                                                     "loglik", 
                                                                                     "conv", 
                                                                                     "npm", 
                                                                                     "AIC", 
                                                                                     "BIC", 
                                                                                     "SABIC", 
                                                                                     "entropy", 
                                                                                     "ICL", 
                                                                                     "ICL1", 
                                                                                     "ICL2", 
                                                                                     "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_cubic_drsc_model_2011_2019)
postprob(threeclass_cubic_drsc_model_2011_2019)
postprob(fourclass_cubic_drsc_model_2011_2019)
postprob(fiveclass_cubic_drsc_model_2011_2019)
postprob(sixclass_cubic_drsc_model_2011_2019)
postprob(sevenclass_cubic_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_cubic_drsc_model_2011_2019 <- postprob(twoclass_cubic_drsc_model_2011_2019)  
postprob_threeclass_cubic_drsc_model_2011_2019 <- postprob(threeclass_cubic_drsc_model_2011_2019)  
postprob_fourclass_cubic_drsc_model_2011_2019 <- postprob(fourclass_cubic_drsc_model_2011_2019)  
postprob_fiveclass_cubic_drsc_model_2011_2019 <- postprob(fiveclass_cubic_drsc_model_2011_2019)  
postprob_sixclass_cubic_drsc_model_2011_2019 <- postprob(sixclass_cubic_drsc_model_2011_2019)
postprob_sevenclass_cubic_drsc_model_2011_2019 <- postprob(sevenclass_cubic_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_cubic_drsc_model_2011_2019 <- list(twoclass_cubic_drsc_model_2011_2019,
                                              threeclass_cubic_drsc_model_2011_2019,
                                              fourclass_cubic_drsc_model_2011_2019,
                                              fiveclass_cubic_drsc_model_2011_2019,
                                              sixclass_cubic_drsc_model_2011_2019,
                                              sevenclass_cubic_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_cubic_drsc_model_2011_2019 <- sapply(model_list_cubic_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_cubic_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_cubic_drsc_model_2011_2019 <- sapply(model_list_cubic_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_cubic_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_cubic_drsc_model_2011_2019 <- sapply(model_list_cubic_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_cubic_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_cubic_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[1], summarytable_cubic_drsc_model_2011_2019$npm[1], summarytable_cubic_drsc_model_2011_2019$G[1], summarytable_cubic_drsc_model_2011_2019$loglik[2], summarytable_cubic_drsc_model_2011_2019$npm[2], summarytable_cubic_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[2], summarytable_cubic_drsc_model_2011_2019$npm[2], summarytable_cubic_drsc_model_2011_2019$G[2], summarytable_cubic_drsc_model_2011_2019$loglik[3], summarytable_cubic_drsc_model_2011_2019$npm[3], summarytable_cubic_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[3], summarytable_cubic_drsc_model_2011_2019$npm[3], summarytable_cubic_drsc_model_2011_2019$G[3], summarytable_cubic_drsc_model_2011_2019$loglik[4], summarytable_cubic_drsc_model_2011_2019$npm[4], summarytable_cubic_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[4], summarytable_cubic_drsc_model_2011_2019$npm[4], summarytable_cubic_drsc_model_2011_2019$G[4], summarytable_cubic_drsc_model_2011_2019$loglik[5], summarytable_cubic_drsc_model_2011_2019$npm[5], summarytable_cubic_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[5], summarytable_cubic_drsc_model_2011_2019$npm[5], summarytable_cubic_drsc_model_2011_2019$G[5], summarytable_cubic_drsc_model_2011_2019$loglik[6], summarytable_cubic_drsc_model_2011_2019$npm[6], summarytable_cubic_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_cubic_drsc_model_2011_2019$ns, summarytable_cubic_drsc_model_2011_2019$loglik[6], summarytable_cubic_drsc_model_2011_2019$npm[6], summarytable_cubic_drsc_model_2011_2019$G[6], summarytable_cubic_drsc_model_2011_2019$loglik[7], summarytable_cubic_drsc_model_2011_2019$npm[7], summarytable_cubic_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_cubic_drsc_model_2011_2019 <- sapply(outputs_vllrt_cubic_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_cubic_drsc_model_2011_2019 <- summarytable_cubic_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_cubic_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_cubic_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_cubic_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_cubic_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_cubic_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_cubic_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_cubic_drsc_model_2011_2019$ns, min(postprob_twoclass_cubic_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_cubic_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_cubic_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_cubic_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_cubic_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_cubic_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_cubic_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_cubic_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_cubic_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_cubic_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_cubic_drsc_model_2011_2019, "summarytable_cubic_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_cubic_drsc_model_2011_2019)




# cubic prop effects prop model DRSC 2011-2019 -----------------------------------------


summarytable_cubic_prop_drsc_model_2011_2019 <-  as.data.frame(lcmm::summarytable(oneclass_cubic_prop_drsc_model_2011_2019,
                                                                                  twoclass_cubic_prop_drsc_model_2011_2019,
                                                                                  threeclass_cubic_prop_drsc_model_2011_2019,
                                                                                  fourclass_cubic_prop_drsc_model_2011_2019,
                                                                                  fiveclass_cubic_prop_drsc_model_2011_2019,
                                                                                  sixclass_cubic_prop_drsc_model_2011_2019,
                                                                                  sevenclass_cubic_prop_drsc_model_2011_2019,
                                                                                  which=c("G", 
                                                                                          "loglik", 
                                                                                          "conv", 
                                                                                          "npm", 
                                                                                          "AIC", 
                                                                                          "BIC", 
                                                                                          "SABIC", 
                                                                                          "entropy", 
                                                                                          "ICL", 
                                                                                          "ICL1", 
                                                                                          "ICL2", 
                                                                                          "%class")))




## Model Adequacy -----------------------------------------------------

## average latent class posterior probability ------------------------------

postprob(twoclass_cubic_prop_drsc_model_2011_2019)
postprob(threeclass_cubic_prop_drsc_model_2011_2019)
postprob(fourclass_cubic_prop_drsc_model_2011_2019)
postprob(fiveclass_cubic_prop_drsc_model_2011_2019)
postprob(sixclass_cubic_prop_drsc_model_2011_2019)
postprob(sevenclass_cubic_prop_drsc_model_2011_2019)


# Assuming postprob() returns a structured list
postprob_twoclass_cubic_prop_drsc_model_2011_2019 <- postprob(twoclass_cubic_prop_drsc_model_2011_2019)  
postprob_threeclass_cubic_prop_drsc_model_2011_2019 <- postprob(threeclass_cubic_prop_drsc_model_2011_2019)  
postprob_fourclass_cubic_prop_drsc_model_2011_2019 <- postprob(fourclass_cubic_prop_drsc_model_2011_2019)  
postprob_fiveclass_cubic_prop_drsc_model_2011_2019 <- postprob(fiveclass_cubic_prop_drsc_model_2011_2019)  
postprob_sixclass_cubic_prop_drsc_model_2011_2019 <- postprob(sixclass_cubic_prop_drsc_model_2011_2019)
postprob_sevenclass_cubic_prop_drsc_model_2011_2019 <- postprob(sevenclass_cubic_prop_drsc_model_2011_2019)  



# Assuming you have a list of your model objects --------------------------


model_list_cubic_prop_drsc_model_2011_2019 <- list(twoclass_cubic_prop_drsc_model_2011_2019,
                                                   threeclass_cubic_prop_drsc_model_2011_2019,
                                                   fourclass_cubic_prop_drsc_model_2011_2019,
                                                   fiveclass_cubic_prop_drsc_model_2011_2019,
                                                   sixclass_cubic_prop_drsc_model_2011_2019,
                                                   sevenclass_cubic_prop_drsc_model_2011_2019)  

## OCC ---------------------------------------------------------------------

# Extract the lower OCC values for each model
lower_occ_values_cubic_prop_drsc_model_2011_2019 <- sapply(model_list_cubic_prop_drsc_model_2011_2019, function(model) {
  occ_values <- LCTMtoolkit(model)$occ
  min(as.numeric(occ_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_occ_values_cubic_prop_drsc_model_2011_2019)


## APPA ---------------------------------------------------------------------


# Extract the lower APPA values for each model
lower_appa_values_cubic_prop_drsc_model_2011_2019 <- sapply(model_list_cubic_prop_drsc_model_2011_2019, function(model) {
  appa_values <- LCTMtoolkit(model)$appa
  min(as.numeric(appa_values[1,]), na.rm = TRUE)  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(lower_appa_values_cubic_prop_drsc_model_2011_2019)


# Mismatch ----------------------------------------------------------------


# Extract the lower OCC values for each model
highest_mismatch_values_cubic_prop_drsc_model_2011_2019 <- sapply(model_list_cubic_prop_drsc_model_2011_2019, function(model) {
  mismatch_values <- LCTMtoolkit(model)$mismatch
  max(as.numeric(mismatch_values[1,]))  # Assuming OCC is in the first row of the OCC matrix
})

# Print the lower OCC values
print(highest_mismatch_values_cubic_prop_drsc_model_2011_2019)


## VLLRT test --------------------------------------------------------------


# Define a function to extract p-value from calc_lrt output
extract_p_value_from_lrt <- function(ns, loglik1, npm1, G1, loglik2, npm2, G2) {
  # Capture the output of calc_lrt
  output <- capture.output(tidyLPA::calc_lrt(ns, loglik1, npm1, G1, loglik2, npm2, G2))
  
  # Combine output into a single string
  output_text <- paste(output, collapse = " ")
  
  # Extract the p-value part from the output
  p_value <- str_extract(output_text, "(?<=p\\s).*")
  return(str_trim(p_value))
}

# Sample inputs for calc_lrt function calls Lo–Mendell–Rubin adjusted LRT  ---------------
outputs_vllrt_cubic_prop_drsc_model_2011_2019 <- list(
  extract_p_value_from_lrt(oneclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[1], summarytable_cubic_prop_drsc_model_2011_2019$npm[1], summarytable_cubic_prop_drsc_model_2011_2019$G[1], summarytable_cubic_prop_drsc_model_2011_2019$loglik[2], summarytable_cubic_prop_drsc_model_2011_2019$npm[2], summarytable_cubic_prop_drsc_model_2011_2019$G[2]),
  extract_p_value_from_lrt(twoclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[2], summarytable_cubic_prop_drsc_model_2011_2019$npm[2], summarytable_cubic_prop_drsc_model_2011_2019$G[2], summarytable_cubic_prop_drsc_model_2011_2019$loglik[3], summarytable_cubic_prop_drsc_model_2011_2019$npm[3], summarytable_cubic_prop_drsc_model_2011_2019$G[3]),
  extract_p_value_from_lrt(threeclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[3], summarytable_cubic_prop_drsc_model_2011_2019$npm[3], summarytable_cubic_prop_drsc_model_2011_2019$G[3], summarytable_cubic_prop_drsc_model_2011_2019$loglik[4], summarytable_cubic_prop_drsc_model_2011_2019$npm[4], summarytable_cubic_prop_drsc_model_2011_2019$G[4]),
  extract_p_value_from_lrt(fourclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[4], summarytable_cubic_prop_drsc_model_2011_2019$npm[4], summarytable_cubic_prop_drsc_model_2011_2019$G[4], summarytable_cubic_prop_drsc_model_2011_2019$loglik[5], summarytable_cubic_prop_drsc_model_2011_2019$npm[5], summarytable_cubic_prop_drsc_model_2011_2019$G[5]),
  extract_p_value_from_lrt(fiveclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[5], summarytable_cubic_prop_drsc_model_2011_2019$npm[5], summarytable_cubic_prop_drsc_model_2011_2019$G[5], summarytable_cubic_prop_drsc_model_2011_2019$loglik[6], summarytable_cubic_prop_drsc_model_2011_2019$npm[6], summarytable_cubic_prop_drsc_model_2011_2019$G[6]),
  extract_p_value_from_lrt(sixclass_cubic_prop_drsc_model_2011_2019$ns, summarytable_cubic_prop_drsc_model_2011_2019$loglik[6], summarytable_cubic_prop_drsc_model_2011_2019$npm[6], summarytable_cubic_prop_drsc_model_2011_2019$G[6], summarytable_cubic_prop_drsc_model_2011_2019$loglik[7], summarytable_cubic_prop_drsc_model_2011_2019$npm[7], summarytable_cubic_prop_drsc_model_2011_2019$G[7])
)

# Assuming outputs is already defined as your object
values_vllrt_outputs_cubic_prop_drsc_model_2011_2019 <- sapply(outputs_vllrt_cubic_prop_drsc_model_2011_2019, function(x) gsub("^= ", "", x))



summarytable_cubic_prop_drsc_model_2011_2019 <- summarytable_cubic_prop_drsc_model_2011_2019 %>% 
  mutate(smallest_class_size_perc = c(100, min(postprob_twoclass_cubic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_threeclass_cubic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_fourclass_cubic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_fiveclass_cubic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_sixclass_cubic_prop_drsc_model_2011_2019[[1]][2,]), min(postprob_sevenclass_cubic_prop_drsc_model_2011_2019[[1]][2,])),
         smallest_class_count = c(oneclass_cubic_prop_drsc_model_2011_2019$ns, min(postprob_twoclass_cubic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_threeclass_cubic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_fourclass_cubic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_fiveclass_cubic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_sixclass_cubic_prop_drsc_model_2011_2019[[1]][1,]), min(postprob_sevenclass_cubic_prop_drsc_model_2011_2019[[1]][1,]) ),
         "Lowest APPA" = c(NA, lower_appa_values_cubic_prop_drsc_model_2011_2019),
         "Highest MMV" =c(NA, highest_mismatch_values_cubic_prop_drsc_model_2011_2019),
         "Lowest OCC" = c(NA, lower_occ_values_cubic_prop_drsc_model_2011_2019),
         VLMRLRT = c(NA, values_vllrt_outputs_cubic_prop_drsc_model_2011_2019)
  ) %>% 
  tibble::as.tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::select(-rowname) 

write.csv(summarytable_cubic_prop_drsc_model_2011_2019, "summarytable_cubic_prop_drsc_model_2011_2019.csv", row.names = F)

View(summarytable_cubic_prop_drsc_model_2011_2019)




summarytable_linear_nre_homocedastic_drsc_model_2011_2023 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_nre_heterocedastic_drsc_model_2011_2023 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_drsc_model_random_intercept_2011_2023 %>% 
  filter(smallest_class_count > 9)%>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_drsc_model_random_intercept_slope_2011_2023 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_quadratic_drsc_model_2011_2023 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_quadratic_prop_drsc_model_2011_2023 %>% 
  filter(smallest_class_count > 9)%>% 
  arrange(BIC)%>% as.data.frame()

summarytable_cubic_drsc_model_2011_2023 %>% 
  #filter(smallest_class_count > 9) %>% 
  arrange(BIC) %>% as.data.frame()

summarytable_cubic_prop_drsc_model_2011_2023 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC) %>% as.data.frame()



summarytable_linear_nre_homocedastic_drsc_model_2011_2019 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_nre_heterocedastic_drsc_model_2011_2019 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_drsc_model_random_intercept_2011_2019 %>% 
  filter(smallest_class_count > 9)%>% 
  arrange(BIC)%>% as.data.frame()

summarytable_linear_drsc_model_random_intercept_slope_2011_2019 %>% 
  filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_quadratic_drsc_model_2011_2019 %>% 
  #filter(smallest_class_count > 9) %>% 
  arrange(BIC)%>% as.data.frame()

summarytable_quadratic_prop_drsc_model_2011_2019 %>% 
  #filter(smallest_class_count > 9)%>% 
  arrange(BIC)%>% as.data.frame()

summarytable_cubic_drsc_model_2011_2019 %>% 
  #filter(smallest_class_count > 9) %>% 
  arrange(BIC) %>% as.data.frame()

summarytable_cubic_prop_drsc_model_2011_2019 %>% 
  #filter(smallest_class_count > 9) %>% 
  arrange(BIC) %>% as.data.frame()


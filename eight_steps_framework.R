
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
  select(comuna2, ano, drs_coverage, dm_coverage) %>% 
  arrange(comuna2, -ano) %>% 
  mutate(ano = max(ano) - ano + 1) %>% 
  dplyr::rename(id = comuna2, 
         year = ano,
         drsc = drs_coverage,
         dgcc = dm_coverage)

# Second period: 2011-2019 (actual calendar years 2011 to 2019)
coverage_long_2011_2019 <- coverage_long %>%
  filter(year %in% 1:9)  # This corresponds to years 1 to 9 (2011-2019)

# Third period: 2020-2023 (actual calendar years 2020 to 2023)
coverage_long_2020_2023 <- coverage_long %>%
  filter(year %in% 10:13)  # This corresponds to years 10 to 23 (2020-2023)


# Step 1 - scoping model --------------------------------------------------



# scoping model 2011-2023 drsc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                              #mixture = ~1+year+I(year^2),
                              random = ~-1,
                              ng = 1,
                              nwg = FALSE, 
                              data=coverage_long_2011_2023,
                              subject = "id")

twoclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                         mixture = drsc~1+year,
                                         random = ~-1,
                                         ng = 2,
                                         nwg = FALSE, 
                                         data=coverage_long_2011_2023,
                                         subject = "id", 
                                         B= oneclass_linear_drsc_model_2011_2023)

threeclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 3,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2023)

fourclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2023)

fiveclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2023)

sixclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2023)

sevenclass_linear_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 7,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2023)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                            #mixture = ~1+year+I(year^2),
                                            random = ~-1,
                                            ng = 1,
                                            nwg = FALSE, 
                                            data=coverage_long_2011_2023,
                                            subject = "id")

twoclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2023)


threeclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 3,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2023)


fourclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 4,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2023)

fiveclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2023)

sixclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 6,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2023)

sevenclass_quadratic_drsc_model_2011_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2023)




## Cubic non random effects model ---------------------------------------


oneclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id")


twoclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2011_2023)


threeclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 3,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2011_2023)

fourclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 4,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2011_2023)


fiveclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_drsc_model_2011_2023)

sixclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2011_2023)


sevenclass_cubic_drsc_model_2011_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 7,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2011_2023)




# scoping model 2011-2019 drsc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                   #mixture = ~1+year+I(year^2),
                                                   random = ~-1,
                                                   ng = 1,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id")

twoclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 2,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2019)

threeclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 3,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2011_2019,
                                                     subject = "id", 
                                                     B= oneclass_linear_drsc_model_2011_2019)

fourclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 4,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2011_2019,
                                                    subject = "id", 
                                                    B= oneclass_linear_drsc_model_2011_2019)

fiveclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2011_2019,
                                                    subject = "id", 
                                                    B= oneclass_linear_drsc_model_2011_2019)

sixclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2011_2019)

sevenclass_linear_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 7,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2011_2019,
                                                     subject = "id", 
                                                     B= oneclass_linear_drsc_model_2011_2019)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id")

twoclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2019)


threeclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2019)


fourclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2019)

fiveclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2011_2019)

sixclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2011_2019)

sevenclass_quadratic_drsc_model_2011_2019 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 6,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2011_2019)




## Cubic non random effects model ---------------------------------------


oneclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id")


twoclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2011_2019)


threeclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 3,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B= oneclass_cubic_drsc_model_2011_2019)

fourclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2011_2019)


fiveclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2011_2019)

sixclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 6,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2011_2019)


sevenclass_cubic_drsc_model_2011_2019 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 7,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B= oneclass_cubic_drsc_model_2011_2019)


# scoping model 2020-2023 drsc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   #mixture = ~1+year+I(year^2),
                                                   random = ~-1,
                                                   ng = 1,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id")

twoclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 2,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2020_2023)

threeclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 3,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2020_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_drsc_model_2020_2023)

fourclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 4,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2020_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_drsc_model_2020_2023)

fiveclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                    mixture = drsc~1+year,
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2020_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_drsc_model_2020_2023)

sixclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                   mixture = drsc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_drsc_model_2020_2023)

sevenclass_linear_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year,
                                                     mixture = drsc~1+year,
                                                     random = ~-1,
                                                     ng = 7,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2020_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_drsc_model_2020_2023)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id")

twoclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2020_2023)


threeclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2020_2023)


fourclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2020_2023)

fiveclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_drsc_model_2020_2023)

sixclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_drsc_model_2020_2023)

sevenclass_quadratic_drsc_model_2020_2023 <- lcmm::hlme(fixed=drsc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 6,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_drsc_model_2020_2023)




## Cubic non random effects model ---------------------------------------


oneclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id")


twoclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2020_2023)


threeclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 3,
                                                    nwg = FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_drsc_model_2020_2023)

fourclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2020_2023)


fiveclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_drsc_model_2020_2023)

sixclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 6,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_drsc_model_2020_2023)


sevenclass_cubic_drsc_model_2020_2023 <- lcmm::hlme(fixed = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = drsc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 7,
                                                    nwg = FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_drsc_model_2020_2023)




#sacar el summary table
#probar si los modelos convergen sin necesidad de darle los parameters


residualplot_step1 <- function(model, nameofoutcome, nameofage, data, 
                               ylimit = c(-50, 50)) 
{
  require(dplyr)
  require(ggplot2)
  require(gridExtra)  # for grid.arrange
  
  k <- model$ng
  preds <- model$pred
  names(preds)[6] <- nameofoutcome
  nameofid <- names(model$pred)[1]
  test <- dplyr::left_join(preds, model$pprob, by = nameofid)
  test <- dplyr::left_join(test, data, by = c(nameofid, nameofoutcome))
  
  plot_list <- list()  # List to store plots
  
  # Loop through each class
  for (i in 1:k) {
    newplotvalues <- test %>% filter(class == i) %>% 
      mutate(Residuals = get(nameofoutcome) - eval(parse(text = paste0("pred_ss", i))))
    
    plotvaluessub <- newplotvalues
    pname <- paste0("p", i)
    
    # Create plot
    p <- ggplot2::ggplot(data = plotvaluessub, aes(x = get(nameofage), y = Residuals, group = class)) +
      theme(axis.text = element_text(size = 16), text = element_text(size = 16)) + 
      geom_point() + 
      stat_summary(fun.y = mean, geom = "line", size = 3, col = "CadetBlue", group = 1) + 
      ggtitle(paste("Residuals in class", i)) + 
      ylim(ylimit) + 
      labs(x = "Year")
    
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  
  # After the loop, display all the plots together in one figure
  if (length(plot_list) > 1) {
    grid.arrange(grobs = plot_list, ncol = 2, nrow=1)  # Arrange plots in a grid (2 columns)
  } else {
    print(plot_list[[1]])  # If only one plot, display it directly
  }
  
  return(plot_list)  # Return list of plots
}


plot_oneclass_linear_drsc_model_2011_2023 <- residualplot_step1(oneclass_linear_drsc_model_2011_2023, 
                    nameofoutcome="drsc",  
                    nameofage = "year",
                    data = coverage_long_2011_2023,
                    ylimit=c(-5,5)) 


residualplot_step1(twoclass_linear_drsc_model_2011_2023, 
                   nameofoutcome="drsc",  
                   nameofage = "year",
                   data = coverage_long_2011_2023,
                   ylimit=c(-5,5)) 

residualplot_step1(threeclass_linear_drsc_model_2011_2023, 
                   nameofoutcome="drsc",  
                   nameofage = "year",
                   data = coverage_long_2011_2023,
                   ylimit=c(-5,5)) 


plot_twoclass_linear_drsc_model_2011_2023







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
                                                      ng = 7,
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
                                                        ng = 7,
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
                                                        ng = 7,
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




# scoping model 2011-2023 dgcc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                              #mixture = ~1+year+I(year^2),
                              random = ~-1,
                              ng = 1,
                              nwg = FALSE, 
                              data=coverage_long_2011_2023,
                              subject = "id")

twoclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                         mixture = dgcc~1+year,
                                         random = ~-1,
                                         ng = 2,
                                         nwg = FALSE, 
                                         data=coverage_long_2011_2023,
                                         subject = "id", 
                                         B= oneclass_linear_dgcc_model_2011_2023)

threeclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 3,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2023)

fourclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2023)

fiveclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2023)

sixclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2023)

sevenclass_linear_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 7,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2023)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                            #mixture = ~1+year+I(year^2),
                                            random = ~-1,
                                            ng = 1,
                                            nwg = FALSE, 
                                            data=coverage_long_2011_2023,
                                            subject = "id")

twoclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2011_2023)


threeclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 3,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2011_2023)


fourclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 4,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_dgcc_model_2011_2023)

fiveclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2011_2023)

sixclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 6,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2011_2023)

sevenclass_quadratic_dgcc_model_2011_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 7,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2011_2023)




## Cubic non random effects model ---------------------------------------


oneclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id")


twoclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2011_2023)


threeclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 3,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2011_2023)

fourclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 4,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2011_2023)


fiveclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_dgcc_model_2011_2023)

sixclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2011_2023)


sevenclass_cubic_dgcc_model_2011_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 7,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2011_2023)



# scoping model 2011-2019 dgcc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   #mixture = ~1+year+I(year^2),
                                                   random = ~-1,
                                                   ng = 1,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id")

twoclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 2,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2019)

threeclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                     mixture = dgcc~1+year,
                                                     random = ~-1,
                                                     ng = 3,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2011_2019,
                                                     subject = "id", 
                                                     B= oneclass_linear_dgcc_model_2011_2019)

fourclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                    mixture = dgcc~1+year,
                                                    random = ~-1,
                                                    ng = 4,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2011_2019,
                                                    subject = "id", 
                                                    B= oneclass_linear_dgcc_model_2011_2019)

fiveclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                    mixture = dgcc~1+year,
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2011_2019,
                                                    subject = "id", 
                                                    B= oneclass_linear_dgcc_model_2011_2019)

sixclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2011_2019,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2011_2019)

sevenclass_linear_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year,
                                                     mixture = dgcc~1+year,
                                                     random = ~-1,
                                                     ng = 7,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2011_2019,
                                                     subject = "id", 
                                                     B= oneclass_linear_dgcc_model_2011_2019)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id")

twoclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2011_2019)


threeclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_dgcc_model_2011_2019)


fourclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2011_2019)

fiveclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2011_2019,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2011_2019)

sixclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2011_2019,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2011_2019)

sevenclass_quadratic_dgcc_model_2011_2019 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 7,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2011_2019,
                                                        subject = "id",
                                                        B=oneclass_quadratic_dgcc_model_2011_2019)




## Cubic non random effects model ---------------------------------------


oneclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id")


twoclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2011_2019)


threeclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 3,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B= oneclass_cubic_dgcc_model_2011_2019)

fourclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2011_2019)


fiveclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE,
                                                   data = coverage_long_2011_2019,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2011_2019)

sixclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 6,
                                                  nwg = FALSE,
                                                  data = coverage_long_2011_2019,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2011_2019)


sevenclass_cubic_dgcc_model_2011_2019 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 7,
                                                    nwg = FALSE,
                                                    data = coverage_long_2011_2019,
                                                    subject = "id",
                                                    B= oneclass_cubic_dgcc_model_2011_2019)


# scoping model 2020-2023 dgcc -------------------------------------------------



## Linear non random effects model -----------------------------------------


oneclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   #mixture = ~1+year+I(year^2),
                                                   random = ~-1,
                                                   ng = 1,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id")

twoclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 2,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2020_2023)

threeclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                     mixture = dgcc~1+year,
                                                     random = ~-1,
                                                     ng = 3,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2020_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_dgcc_model_2020_2023)

fourclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                    mixture = dgcc~1+year,
                                                    random = ~-1,
                                                    ng = 4,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2020_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_dgcc_model_2020_2023)

fiveclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                    mixture = dgcc~1+year,
                                                    random = ~-1,
                                                    ng = 5,
                                                    nwg = FALSE, 
                                                    data=coverage_long_2020_2023,
                                                    subject = "id", 
                                                    B= oneclass_linear_dgcc_model_2020_2023)

sixclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                   mixture = dgcc~1+year,
                                                   random = ~-1,
                                                   ng = 6,
                                                   nwg = FALSE, 
                                                   data=coverage_long_2020_2023,
                                                   subject = "id", 
                                                   B= oneclass_linear_dgcc_model_2020_2023)

sevenclass_linear_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year,
                                                     mixture = dgcc~1+year,
                                                     random = ~-1,
                                                     ng = 7,
                                                     nwg = FALSE, 
                                                     data=coverage_long_2020_2023,
                                                     subject = "id", 
                                                     B= oneclass_linear_dgcc_model_2020_2023)


## Quadratic non random effects model ---------------------------------------


oneclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      #mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 1,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id")

twoclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 2,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2020_2023)


threeclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 3,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_dgcc_model_2020_2023)


fourclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 4,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2020_2023)

fiveclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                       mixture = ~1+year+I(year^2),
                                                       random = ~-1,
                                                       ng = 5,
                                                       nwg = FALSE, 
                                                       data=coverage_long_2020_2023,
                                                       subject = "id",
                                                       B=oneclass_quadratic_dgcc_model_2020_2023)

sixclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                      mixture = ~1+year+I(year^2),
                                                      random = ~-1,
                                                      ng = 6,
                                                      nwg = FALSE, 
                                                      data=coverage_long_2020_2023,
                                                      subject = "id",
                                                      B=oneclass_quadratic_dgcc_model_2020_2023)

sevenclass_quadratic_dgcc_model_2020_2023 <- lcmm::hlme(fixed=dgcc~1+year+I(year^2),
                                                        mixture = ~1+year+I(year^2),
                                                        random = ~-1,
                                                        ng = 6,
                                                        nwg = FALSE, 
                                                        data=coverage_long_2020_2023,
                                                        subject = "id",
                                                        B=oneclass_quadratic_dgcc_model_2020_2023)

## Cubic non random effects model ---------------------------------------


oneclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 1,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id")


twoclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 2,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2020_2023)


threeclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 3,
                                                    nwg = FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_dgcc_model_2020_2023)

fourclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 4,
                                                   nwg = FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2020_2023)


fiveclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                   random = ~-1,
                                                   ng = 5,
                                                   nwg = FALSE,
                                                   data = coverage_long_2020_2023,
                                                   subject = "id",
                                                   B= oneclass_cubic_dgcc_model_2020_2023)

sixclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                  random = ~-1,
                                                  ng = 6,
                                                  nwg = FALSE,
                                                  data = coverage_long_2020_2023,
                                                  subject = "id",
                                                  B= oneclass_cubic_dgcc_model_2020_2023)


sevenclass_cubic_dgcc_model_2020_2023 <- lcmm::hlme(fixed = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    mixture = dgcc ~ 1 + year + I(year^2) + I(year^3),
                                                    random = ~-1,
                                                    ng = 7,
                                                    nwg = FALSE,
                                                    data = coverage_long_2020_2023,
                                                    subject = "id",
                                                    B= oneclass_cubic_dgcc_model_2020_2023)








residualplot_step1 <- function(model, nameofoutcome, nameofage, data, 
                               ylimit = c(-50, 50), save_path = NULL, model_name = NULL) {
  require(dplyr)
  require(ggplot2)
  require(gridExtra)  # for grid.arrange
  
  # Ensure model_name is provided
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
  
  # Define x-axis limits based on data range
  xlim_range <- range(data[[nameofage]], na.rm = TRUE)
  
  # Loop through each class
  for (i in 1:k) {
    newplotvalues <- test %>% 
      filter(class == i) %>% 
      mutate(Residuals = get(nameofoutcome) - eval(parse(text = paste0("pred_ss", i))))
    
    plotvaluessub <- newplotvalues
    
    # Create the plot
    p <- ggplot2::ggplot(data = plotvaluessub, aes(x = get(nameofage), y = Residuals, group = class)) +
      theme(axis.text = element_text(size = 8), text = element_text(size = 8)) + 
      geom_point(size = 0.1) + 
      stat_summary(fun = mean, geom = "line", size = 1, col = "CadetBlue", group = 1) + 
      ggtitle(paste("Residuals in class", i)) + 
      ylim(ylimit) + 
      xlim(xlim_range) +  # Standardise x-axis limits
      labs(x = "Year") + 
      coord_fixed(ratio = 1) +  # Ensure consistent aspect ratio
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
    
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  
  # Combine all the plots into one grid
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 4, nrow = ceiling(k / 4))
  
  # Save the combined plot if a save path is provided
  if (!is.null(save_path)) {
    # Convert combined plot to a grob
    g <- arrangeGrob(grobs = plot_list, ncol = 4, nrow = ceiling(k / 4))
    ggsave(filename = file.path(save_path, paste0(model_name, ".jpeg")), plot = g,
           width = 14, height = 8, dpi = 300)
  }
  
  return(combined_plot)  # Return the combined plot
}


# Directory where to save the plots
save_directory <- getwd()  # Change this path if needed

# Ensure the directory exists
if (!dir.exists(save_directory)) {
  dir.create(save_directory, recursive = TRUE)
}

# Apply the function to all models
linear_drsc_models_2011_2023 <- list(
  oneclass_linear_drsc_model_2011_2023,
  twoclass_linear_drsc_model_2011_2023,
  threeclass_linear_drsc_model_2011_2023,
  fourclass_linear_drsc_model_2011_2023,
  fiveclass_linear_drsc_model_2011_2023,
  sixclass_linear_drsc_model_2011_2023,
  sevenclass_linear_drsc_model_2011_2023
)

model_names_linear_drsc_model_2011_2023 <- c(
  "oneclass_linear_drsc_model_2011_2023",
  "twoclass_linear_drsc_model_2011_2023",
  "threeclass_linear_drsc_model_2011_2023",
  "fourclass_linear_drsc_model_2011_2023",
  "fiveclass_linear_drsc_model_2011_2023",
  "sixclass_linear_drsc_model_2011_2023",
  "sevenclass_linear_drsc_model_2011_2023"
)

for (i in seq_along(linear_drsc_models_2011_2023)) {
  residualplot_step1(
    model = linear_drsc_models_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_drsc_model_2011_2023[i]
  )
}


# Apply the function to all models
quadratic_drsc_models_2011_2023 <- list(
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

for (i in seq_along(quadratic_drsc_models_2011_2023)) {
  residualplot_step1(
    model = quadratic_drsc_models_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_drsc_model_2011_2023[i]
  )
}

# Apply the function to all models
cubic_drsc_models_2011_2023 <- list(
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

for (i in seq_along(cubic_drsc_models_2011_2023)) {
  residualplot_step1(
    model = cubic_drsc_models_2011_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_drsc_model_2011_2023[i]
  )
}


# Apply the function to all models
linear_drsc_models_2011_2019 <- list(
  oneclass_linear_drsc_model_2011_2019,
  twoclass_linear_drsc_model_2011_2019,
  threeclass_linear_drsc_model_2011_2019,
  fourclass_linear_drsc_model_2011_2019,
  fiveclass_linear_drsc_model_2011_2019,
  sixclass_linear_drsc_model_2011_2019,
  sevenclass_linear_drsc_model_2011_2019
)

model_names_linear_drsc_model_2011_2019 <- c(
  "oneclass_linear_drsc_model_2011_2019",
  "twoclass_linear_drsc_model_2011_2019",
  "threeclass_linear_drsc_model_2011_2019",
  "fourclass_linear_drsc_model_2011_2019",
  "fiveclass_linear_drsc_model_2011_2019",
  "sixclass_linear_drsc_model_2011_2019",
  "sevenclass_linear_drsc_model_2011_2019"
)

for (i in seq_along(linear_drsc_models_2011_2019)) {
  residualplot_step1(
    model = linear_drsc_models_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_drsc_model_2011_2019[i]
  )
}


# Apply the function to all models
quadratic_drsc_models_2011_2019 <- list(
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

for (i in seq_along(quadratic_drsc_models_2011_2019)) {
  residualplot_step1(
    model = quadratic_drsc_models_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_drsc_model_2011_2019[i]
  )
}

# Apply the function to all models
cubic_drsc_models_2011_2019 <- list(
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

for (i in seq_along(cubic_drsc_models_2011_2019)) {
  residualplot_step1(
    model = cubic_drsc_models_2011_2019[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_drsc_model_2011_2019[i]
  )
}

# Apply the function to all models
linear_drsc_models_2020_2023 <- list(
  oneclass_linear_drsc_model_2020_2023,
  twoclass_linear_drsc_model_2020_2023,
  threeclass_linear_drsc_model_2020_2023,
  fourclass_linear_drsc_model_2020_2023,
  fiveclass_linear_drsc_model_2020_2023,
  sixclass_linear_drsc_model_2020_2023,
  sevenclass_linear_drsc_model_2020_2023
)

model_names_linear_drsc_model_2020_2023 <- c(
  "oneclass_linear_drsc_model_2020_2023",
  "twoclass_linear_drsc_model_2020_2023",
  "threeclass_linear_drsc_model_2020_2023",
  "fourclass_linear_drsc_model_2020_2023",
  "fiveclass_linear_drsc_model_2020_2023",
  "sixclass_linear_drsc_model_2020_2023",
  "sevenclass_linear_drsc_model_2020_2023"
)

for (i in seq_along(linear_drsc_models_2020_2023)) {
  residualplot_step1(
    model = linear_drsc_models_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_drsc_model_2020_2023[i]
  )
}


# Apply the function to all models
quadratic_drsc_models_2020_2023 <- list(
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

for (i in seq_along(quadratic_drsc_models_2020_2023)) {
  residualplot_step1(
    model = quadratic_drsc_models_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_drsc_model_2020_2023[i]
  )
}

# Apply the function to all models
cubic_drsc_models_2020_2023 <- list(
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

for (i in seq_along(cubic_drsc_models_2020_2023)) {
  residualplot_step1(
    model = cubic_drsc_models_2020_2023[[i]],
    nameofoutcome = "drsc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_drsc_model_2020_2023[i]
  )
}


# Apply the function to all models
linear_dgcc_models_2011_2023 <- list(
  oneclass_linear_dgcc_model_2011_2023,
  twoclass_linear_dgcc_model_2011_2023,
  threeclass_linear_dgcc_model_2011_2023,
  fourclass_linear_dgcc_model_2011_2023,
  fiveclass_linear_dgcc_model_2011_2023,
  sixclass_linear_dgcc_model_2011_2023,
  sevenclass_linear_dgcc_model_2011_2023
)

model_names_linear_dgcc_model_2011_2023 <- c(
  "oneclass_linear_dgcc_model_2011_2023",
  "twoclass_linear_dgcc_model_2011_2023",
  "threeclass_linear_dgcc_model_2011_2023",
  "fourclass_linear_dgcc_model_2011_2023",
  "fiveclass_linear_dgcc_model_2011_2023",
  "sixclass_linear_dgcc_model_2011_2023",
  "sevenclass_linear_dgcc_model_2011_2023"
)

for (i in seq_along(linear_dgcc_models_2011_2023)) {
  residualplot_step1(
    model = linear_dgcc_models_2011_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_dgcc_model_2011_2023[i]
  )
}


# Apply the function to all models
quadratic_dgcc_models_2011_2023 <- list(
  oneclass_quadratic_dgcc_model_2011_2023,
  twoclass_quadratic_dgcc_model_2011_2023,
  threeclass_quadratic_dgcc_model_2011_2023,
  fourclass_quadratic_dgcc_model_2011_2023,
  fiveclass_quadratic_dgcc_model_2011_2023,
  sixclass_quadratic_dgcc_model_2011_2023,
  sevenclass_quadratic_dgcc_model_2011_2023
)

model_names_quadratic_dgcc_model_2011_2023 <- c(
  "oneclass_quadratic_dgcc_model_2011_2023",
  "twoclass_quadratic_dgcc_model_2011_2023",
  "threeclass_quadratic_dgcc_model_2011_2023",
  "fourclass_quadratic_dgcc_model_2011_2023",
  "fiveclass_quadratic_dgcc_model_2011_2023",
  "sixclass_quadratic_dgcc_model_2011_2023",
  "sevenclass_quadratic_dgcc_model_2011_2023"
)

for (i in seq_along(quadratic_dgcc_models_2011_2023)) {
  residualplot_step1(
    model = quadratic_dgcc_models_2011_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_dgcc_model_2011_2023[i]
  )
}

# Apply the function to all models
cubic_dgcc_models_2011_2023 <- list(
  oneclass_cubic_dgcc_model_2011_2023,
  twoclass_cubic_dgcc_model_2011_2023,
  threeclass_cubic_dgcc_model_2011_2023,
  fourclass_cubic_dgcc_model_2011_2023,
  fiveclass_cubic_dgcc_model_2011_2023,
  sixclass_cubic_dgcc_model_2011_2023,
  sevenclass_cubic_dgcc_model_2011_2023
)

model_names_cubic_dgcc_model_2011_2023 <- c(
  "oneclass_cubic_dgcc_model_2011_2023",
  "twoclass_cubic_dgcc_model_2011_2023",
  "threeclass_cubic_dgcc_model_2011_2023",
  "fourclass_cubic_dgcc_model_2011_2023",
  "fiveclass_cubic_dgcc_model_2011_2023",
  "sixclass_cubic_dgcc_model_2011_2023",
  "sevenclass_cubic_dgcc_model_2011_2023"
)

for (i in seq_along(cubic_dgcc_models_2011_2023)) {
  residualplot_step1(
    model = cubic_dgcc_models_2011_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_dgcc_model_2011_2023[i]
  )
}


# Apply the function to all models
linear_dgcc_models_2011_2019 <- list(
  oneclass_linear_dgcc_model_2011_2019,
  twoclass_linear_dgcc_model_2011_2019,
  threeclass_linear_dgcc_model_2011_2019,
  fourclass_linear_dgcc_model_2011_2019,
  fiveclass_linear_dgcc_model_2011_2019,
  sixclass_linear_dgcc_model_2011_2019,
  sevenclass_linear_dgcc_model_2011_2019
)

model_names_linear_dgcc_model_2011_2019 <- c(
  "oneclass_linear_dgcc_model_2011_2019",
  "twoclass_linear_dgcc_model_2011_2019",
  "threeclass_linear_dgcc_model_2011_2019",
  "fourclass_linear_dgcc_model_2011_2019",
  "fiveclass_linear_dgcc_model_2011_2019",
  "sixclass_linear_dgcc_model_2011_2019",
  "sevenclass_linear_dgcc_model_2011_2019"
)

for (i in seq_along(linear_dgcc_models_2011_2019)) {
  residualplot_step1(
    model = linear_dgcc_models_2011_2019[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_dgcc_model_2011_2019[i]
  )
}


# Apply the function to all models
quadratic_dgcc_models_2011_2019 <- list(
  oneclass_quadratic_dgcc_model_2011_2019,
  twoclass_quadratic_dgcc_model_2011_2019,
  threeclass_quadratic_dgcc_model_2011_2019,
  fourclass_quadratic_dgcc_model_2011_2019,
  fiveclass_quadratic_dgcc_model_2011_2019,
  sixclass_quadratic_dgcc_model_2011_2019,
  sevenclass_quadratic_dgcc_model_2011_2019
)

model_names_quadratic_dgcc_model_2011_2019 <- c(
  "oneclass_quadratic_dgcc_model_2011_2019",
  "twoclass_quadratic_dgcc_model_2011_2019",
  "threeclass_quadratic_dgcc_model_2011_2019",
  "fourclass_quadratic_dgcc_model_2011_2019",
  "fiveclass_quadratic_dgcc_model_2011_2019",
  "sixclass_quadratic_dgcc_model_2011_2019",
  "sevenclass_quadratic_dgcc_model_2011_2019"
)

for (i in seq_along(quadratic_dgcc_models_2011_2019)) {
  residualplot_step1(
    model = quadratic_dgcc_models_2011_2019[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_dgcc_model_2011_2019[i]
  )
}

# Apply the function to all models
cubic_dgcc_models_2011_2019 <- list(
  oneclass_cubic_dgcc_model_2011_2019,
  twoclass_cubic_dgcc_model_2011_2019,
  threeclass_cubic_dgcc_model_2011_2019,
  fourclass_cubic_dgcc_model_2011_2019,
  fiveclass_cubic_dgcc_model_2011_2019,
  sixclass_cubic_dgcc_model_2011_2019,
  sevenclass_cubic_dgcc_model_2011_2019
)

model_names_cubic_dgcc_model_2011_2019 <- c(
  "oneclass_cubic_dgcc_model_2011_2019",
  "twoclass_cubic_dgcc_model_2011_2019",
  "threeclass_cubic_dgcc_model_2011_2019",
  "fourclass_cubic_dgcc_model_2011_2019",
  "fiveclass_cubic_dgcc_model_2011_2019",
  "sixclass_cubic_dgcc_model_2011_2019",
  "sevenclass_cubic_dgcc_model_2011_2019"
)

for (i in seq_along(cubic_dgcc_models_2011_2019)) {
  residualplot_step1(
    model = cubic_dgcc_models_2011_2019[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2011_2019,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_dgcc_model_2011_2019[i]
  )
}

# Apply the function to all models
linear_dgcc_models_2020_2023 <- list(
  oneclass_linear_dgcc_model_2020_2023,
  twoclass_linear_dgcc_model_2020_2023,
  threeclass_linear_dgcc_model_2020_2023,
  fourclass_linear_dgcc_model_2020_2023,
  fiveclass_linear_dgcc_model_2020_2023,
  sixclass_linear_dgcc_model_2020_2023,
  sevenclass_linear_dgcc_model_2020_2023
)

model_names_linear_dgcc_model_2020_2023 <- c(
  "oneclass_linear_dgcc_model_2020_2023",
  "twoclass_linear_dgcc_model_2020_2023",
  "threeclass_linear_dgcc_model_2020_2023",
  "fourclass_linear_dgcc_model_2020_2023",
  "fiveclass_linear_dgcc_model_2020_2023",
  "sixclass_linear_dgcc_model_2020_2023",
  "sevenclass_linear_dgcc_model_2020_2023"
)

for (i in seq_along(linear_dgcc_models_2020_2023)) {
  residualplot_step1(
    model = linear_dgcc_models_2020_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_linear_dgcc_model_2020_2023[i]
  )
}


# Apply the function to all models
quadratic_dgcc_models_2020_2023 <- list(
  oneclass_quadratic_dgcc_model_2020_2023,
  twoclass_quadratic_dgcc_model_2020_2023,
  threeclass_quadratic_dgcc_model_2020_2023,
  fourclass_quadratic_dgcc_model_2020_2023,
  fiveclass_quadratic_dgcc_model_2020_2023,
  sixclass_quadratic_dgcc_model_2020_2023,
  sevenclass_quadratic_dgcc_model_2020_2023
)

model_names_quadratic_dgcc_model_2020_2023 <- c(
  "oneclass_quadratic_dgcc_model_2020_2023",
  "twoclass_quadratic_dgcc_model_2020_2023",
  "threeclass_quadratic_dgcc_model_2020_2023",
  "fourclass_quadratic_dgcc_model_2020_2023",
  "fiveclass_quadratic_dgcc_model_2020_2023",
  "sixclass_quadratic_dgcc_model_2020_2023",
  "sevenclass_quadratic_dgcc_model_2020_2023"
)

for (i in seq_along(quadratic_dgcc_models_2020_2023)) {
  residualplot_step1(
    model = quadratic_dgcc_models_2020_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_quadratic_dgcc_model_2020_2023[i]
  )
}

# Apply the function to all models
cubic_dgcc_models_2020_2023 <- list(
  oneclass_cubic_dgcc_model_2020_2023,
  twoclass_cubic_dgcc_model_2020_2023,
  threeclass_cubic_dgcc_model_2020_2023,
  fourclass_cubic_dgcc_model_2020_2023,
  fiveclass_cubic_dgcc_model_2020_2023,
  sixclass_cubic_dgcc_model_2020_2023,
  sevenclass_cubic_dgcc_model_2020_2023
)

model_names_cubic_dgcc_model_2020_2023 <- c(
  "oneclass_cubic_dgcc_model_2020_2023",
  "twoclass_cubic_dgcc_model_2020_2023",
  "threeclass_cubic_dgcc_model_2020_2023",
  "fourclass_cubic_dgcc_model_2020_2023",
  "fiveclass_cubic_dgcc_model_2020_2023",
  "sixclass_cubic_dgcc_model_2020_2023",
  "sevenclass_cubic_dgcc_model_2020_2023"
)

for (i in seq_along(cubic_dgcc_models_2020_2023)) {
  residualplot_step1(
    model = cubic_dgcc_models_2020_2023[[i]],
    nameofoutcome = "dgcc",
    nameofage = "year",
    data = coverage_long_2020_2023,
    ylimit = c(-5, 5),
    save_path = save_directory,
    model_name = model_names_cubic_dgcc_model_2020_2023[i]
  )
}


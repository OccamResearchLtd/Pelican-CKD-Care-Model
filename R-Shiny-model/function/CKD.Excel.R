############################
# CKD model validation     #
# Programmed by Ziyi Lin   #
############################
# Equivalent to Excel v.4.2
#rm(list=ls())

# control base case or sensitivity analysis
# 1: base
# 2: prob
# 2n+1: lower
# 2n+2: upper
# n = 1, 2, ..., 12 for 12 sets of variables
# 1. Cost of CKD health states
# 2. Utility of CKD health states
# 3. RR HF - CV Meds
# 4. RR MI - CV Meds
# 5. RR Stroke - CV Meds
# 6. Annual price - CV Meds
# 7. RR HF - SGLT2i
# 8. RR MI - SGLT2i
# 9. RR Stroke - SGLT2i
# 10. eGFR effect - SGLT2i
# 11. Annual price - SGLT2i
# 12. Others



CKD_model <- function(analysis, input_timehorizon, 
                      input_update_sub, input_egfrstage, input_albstage, input_age,
                      input_dia, input_sex, input_hyp,
                      intervention, comparator){ 
  # intervention and comparator:
  # Usual care
  # CV risk factor management
  # CV risk factor management for all eligible
  # CV risk factor management for all eligible and 50% CKD management
  # CV risk factor management for all eligible and full CKD management
  # CV risk factor management for all eligible and full CKD management, no delay

  #script_path <- rstudioapi::getActiveDocumentContext()$path
  #setwd(dirname(script_path)) #my current Working Directory
  
  ## model set-up ----
  # library(readr)
  # library(tidyverse)
  # library(rstudioapi)
  # library(htmlTable)
  # analysis <- 1
  # input_timehorizon <- "10 years"
  # ## options: 1 year, 2 years, 3 years, 4 years, 5 years, 10 years, 15 years, 20 years, 40 years, and 70 years (lifetime)
  # ## baseline eGFR albuminuria and diabetes
  # input_egfrstage <- "G1 (≥90)" # unit: mL/min/1.73m^2
  # input_albstage <- "A1 (<30)" # unit: mg/g
  # input_age <- 57
  # input_dia <- 0.445 # % diabetes
  # input_sex <- 0.409 # % female
  # input_hyp <- 0.861 # % hypertensive
  # 
  # ## Other function inputs (for debug)
  # intervention <- "CV risk factor management for all eligible and full CKD management, no delay"
  # comparator <- "Usual care"
  # 
  # input_update_sub <- T
  
  effect_sglt_RR <- T # how the effect of SGLT2 inhibitors are implemented
  ## if T, RRs applied to MI, stroke, and HF risk from a meta-analysis (Kaze et al.) 
  ## if F, risk of these events is impacted by treated eGFR

  disc.outcomes <- 0.03 # to be added to the interface
  disc.cost <- 0.03 # to be added to the interface
  
  n_cohort <- 1
  n_cycles <- 71
  timehorizon_year <- as.integer(sub(" year| years| years \\(lifetime\\)", "", input_timehorizon))
  n_cycles_short <- timehorizon_year + 1

   
  # The function for parameters following the gamma distribution
  ## to generate the lower and upper bounds for one-way sensitivity analysis
  ## and to generate the random value for probabilistic sensitivity analysis
  source("function/SA_gamma.R")
  
  # The function for parameters following the beta distribution
  ## to generate the lower and upper bounds for one-way sensitivity analysis
  ## and to generate the random value for probabilistic sensitivity analysis
  source("function/SA_beta.R")
  
  # The function for parameters following the normal distribution
  ## to generate the lower and upper bounds for one-way sensitivity analysis
  ## and to generate the random value for probabilistic sensitivity analysis
  source("function/SA_norm.R")
  
  # The function for parameters following the log-normal distribution
  ## to generate the lower and upper bounds for one-way sensitivity analysis
  ## and to generate the random value for probabilistic sensitivity analysis
  source("function/SA_lognorm.R")
  
  # Cost inputs ----
  # Consumer Price Index for Medical Care
  # Source: https://www.bls.gov/cpi/data.htm
  CPItable <- read_csv("data/CPI.csv", col_types = cols(Year = col_integer(), 
                                                        CPI = col_double()))
  CPItable$IFactor <- CPItable$CPI[CPItable$Year==2023] / CPItable$CPI
  
  infl_to_latest <- function(cost,year,to_year=2023){
    IF <- CPItable$IFactor[CPItable$Year==year]
    return(cost*IF)
  }
  
  # Source: Optum Data 2022
  cost_ckdcvd <- read_csv("data/cost_ckdcvd.csv", 
                          col_types = cols(State = col_character(), 
                                           FirstYear = col_character(), 
                                           Age = col_character(), 
                                           Stage = col_character(), 
                                           Diabetes = col_character(), 
                                           N = col_integer(), 
                                           Mean = col_double(), 
                                           SD = col_double()))
  
  l_state <- unique(cost_ckdcvd$State)
  n_state <- length(l_state)

  l_first_year <- unique(cost_ckdcvd$FirstYear)
  n_fy <- length(l_first_year)

  l_age65 <- unique(cost_ckdcvd$Age)
  n_age65 <- length(l_age65)

  l_ckdstage <- unique(cost_ckdcvd$Stage)
  n_ckdstage <- length(l_ckdstage)

  l_t2d <- unique(cost_ckdcvd$Diabetes)
  n_t2d <- length(l_t2d)
  
  cost_ckdcvd <- transform(cost_ckdcvd,
                           State = factor(State, levels = l_state),
                           FirstYear = factor(FirstYear, levels = l_first_year),
                           Age = factor(Age, levels = l_age65),
                           Stage = factor(Stage, levels = l_ckdstage),
                           Diabetes = factor(Diabetes, levels = l_t2d),
                           SE = SD / sqrt(N),
                           Year = as.integer(2022),
                           N = NULL,
                           SD = NULL)
  
  
  # Source: Nichols et al., 2020. supplement
  cost_esrd <- read_csv("data/cost_esrd.csv", 
                        col_types = cols(Condition = col_character(),
                                         G5 = col_double(),
                                         G4 = col_double()))
  
  cost_esrd <- transform(cost_esrd,
                         State = ifelse(
                           Condition %in% c("ESRD alone", "ESRD & DM"), "ckd only",
                           ifelse(
                             Condition %in% c("MI", "MI + DM"), "mi", 
                             ifelse(
                               Condition %in% c("Stroke", "Stroke + DM"), "str", 
                               "hf"))),
                         Diabetes = ifelse(
                           Condition %in% c("ESRD & DM", "MI + DM", "Stroke + DM", "HF + DM"), "baseline t2d",
                           "no baseline t2d"),
                         Stage = "g4",
                         Multiplier = G5 / G4,
                         Condition = NULL,
                         G5 = NULL,
                         G4 = NULL
                         )
  
  cost_esrd_lt65_y1 <- transform(cost_esrd, 
                                 Age = "lt 65", 
                                 FirstYear = "first")
  cost_esrd_ge65_y1 <- transform(cost_esrd, 
                                 Age = "ge 65", 
                                 FirstYear = "first")
  cost_esrd_lt65_sub <- transform(cost_esrd, 
                                  Age = "lt 65", 
                                  FirstYear = "subsequent")
  cost_esrd_ge65_sub <- transform(cost_esrd, 
                                  Age = "ge 65", 
                                  FirstYear = "subsequent")

  cost_esrd_new <- rbind(cost_esrd_lt65_y1, 
                         cost_esrd_ge65_y1, 
                         cost_esrd_lt65_sub, 
                         cost_esrd_ge65_sub)
  
  cost_esrd_new2 <- merge(cost_ckdcvd, cost_esrd_new, 
                          by = c("State", "Diabetes", "Stage", "Age", "FirstYear"), 
                          all = FALSE)
  cost_esrd_new2 <- transform(cost_esrd_new2,
                              Mean = Mean * Multiplier,
                              SE = SE * Multiplier,
                              Multiplier = NULL,
                              Stage = "g5")
  
  cost_ckdcvd <- rbind(cost_ckdcvd, cost_esrd_new2)
  
  # sort to crosscheck with the Excel model
  # View(cost_ckdcvd[order(cost_ckdcvd$FirstYear, cost_ckdcvd$Age, 
  #                        cost_ckdcvd$Diabetes, cost_ckdcvd$State),])
  
  
  # adjust for inflation
  cost_ckdcvd$Mean <- mapply(infl_to_latest, cost = cost_ckdcvd$Mean, year = cost_ckdcvd$Year)
  cost_ckdcvd$SE <- mapply(infl_to_latest, cost = cost_ckdcvd$SE, year = cost_ckdcvd$Year)
  cost_ckdcvd$Year <- NULL
  
  
  rm(cost_esrd,
     cost_esrd_lt65_y1, 
     cost_esrd_ge65_y1, 
     cost_esrd_lt65_sub, 
     cost_esrd_ge65_sub,
     cost_esrd_new,
     cost_esrd_new2)
  
  
  # Source: Nichols et al., 2020
  # NB: This is just the mean "no CKD" cost in their data, 
  # hard to identify cost for stage 2
  c_g12_input <- c("Mean" = 5631, "Lower" = 5497, "Upper" = 5768, "Year" = 2019)
  
  # adjust for inflation
  c_g12_input <- mapply(infl_to_latest, 
                        cost = c_g12_input[c("Mean", "Lower", "Upper")],
                        year = c_g12_input["Year"])
  
  # Source: Hoerger et al., 2010
  c_add_a3_input <- c("Mean" = 4854, "SE" = 4854 / 10, "Year" = 2006)

  # Source: Ward et al., 2014
  c_mi_30d_input <- c("Mean" = 31482, "SE" = 31482 / 20, "Year" = 2012)
  c_str_30d_input <- c("Mean" = 18835, "SE" = 18835 / 20, "Year" = 2012)
  c_hf_30d_input <- c("Mean" = 16950, "SE" = 16950 / 20, "Year" = 2012)
  
  # Source: Betts et al., 2021
  c_esrd_d_4m_input <- c("Mean" = 87538, "SE" = 6691, "Year" = 2020)
  c_esrd_t_4m_input <- c("Mean" = 124271, "SE" = 37052, "Year" = 2020)
  
  c_esrd_d_sub4m_input <- c("Mean" = 49573, "SE" = 6691, "Year" = 2020)
  c_esrd_t_sub4m_input <- c("Mean" = 7079, "SE" = 1422, "Year" = 2020)
  
  # adjust for inflation
  c_ms_vectors_input <- list(c_add_a3_input = c_add_a3_input,
                             
                             c_mi_30d_input = c_mi_30d_input, 
                             c_str_30d_input = c_str_30d_input, 
                             c_hf_30d_input = c_hf_30d_input,
                             
                             c_esrd_d_4m_input = c_esrd_d_4m_input, 
                             c_esrd_t_4m_input = c_esrd_t_4m_input,
                             c_esrd_d_sub4m_input = c_esrd_d_sub4m_input, 
                             c_esrd_t_sub4m_input = c_esrd_t_sub4m_input)
  
  c_ms_vectors_input <- lapply(c_ms_vectors_input, function(x) {
    x["Mean"] <- infl_to_latest(cost = x["Mean"], year = x["Year"])
    x["SE"] <- infl_to_latest(cost = x["SE"], year = x["Year"])
    return(x[c("Mean", "SE")])
  })
  
  # The live value (in use for base or SA)
  
  c_g12_live <- SA_gamma(vec.mean_l_u = c_g12_input, n = 1)[analysis] # gamma distribution, n = 1
  
  c_ms_vectors_live <- lapply(c_ms_vectors_input, function(x) {
    x <- SA_gamma(vec.mean_se = x, n = 1)[analysis]
    return(x)
  }) # gamma distribution, n = 1
  
  # Replace "input" with "live" in the names
  names_c_ms_vectors_live <- gsub("input", "live", names(c_ms_vectors_live))
  
  # Assign the new names back to the list
  names(c_ms_vectors_live) <- names_c_ms_vectors_live
  
  list2env(c_ms_vectors_live, .GlobalEnv)
  rm(names_c_ms_vectors_live)
  
  cost_ckdcvd$Live <- sapply(1:nrow(cost_ckdcvd), 
                             function(i) 
                               SA_gamma(vec.mean_se = as.vector(
                                 cost_ckdcvd[i,c("Mean", "SE")]),
                                 n = 1))[analysis,]  # gamma distribution, n = 1
  
  
  ####################### Smoothing function for COST ##########################
  coefficients_cost_CKD <- read_csv("data/Model_outputs/coefficients_cost_CKD.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_cost_CKD <- read_csv("data/Model_outputs/vcov_cost_CKD.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_cost_HF <- read_csv("data/Model_outputs/coefficients_cost_HF.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_cost_HF <- read_csv("data/Model_outputs/vcov_cost_HF.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_cost_MI <- read_csv("data/Model_outputs/coefficients_cost_MI.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_cost_MI <- read_csv("data/Model_outputs/vcov_cost_MI.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_cost_Stroke <- read_csv("data/Model_outputs/coefficients_cost_Stroke.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_cost_Stroke <- read_csv("data/Model_outputs/vcov_cost_Stroke.csv") %>%
    rename(parameter = `...1`)
  
  fn_predict_cost <- function(event = "CKD", 
                              eGFR = c(45,50), 
                              UACR = 150,
                              Age_LT65 = 1, 
                              Baseline_T2D = 1, 
                              Year1 = 1,
                              to_year = 2023) {
    inveGFR = 1 / eGFR
    if (event == "CKD") {
      coefficients_used <- coefficients_cost_CKD %>%
        pull(coef)
      
      value_used <- tibble(
        1, inveGFR, Age_LT65, Baseline_T2D, Year1, 
        inveGFR*Age_LT65, inveGFR*Year1, Age_LT65*Baseline_T2D,
        Age_LT65*Year1, 
        inveGFR*Age_LT65*Year1
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "HF") {
      coefficients_used <- coefficients_cost_HF %>%
        pull(coef) %>%
        as.matrix()
      
      value_used <- tibble(
        1, inveGFR, Age_LT65, Baseline_T2D, Year1, 
        inveGFR*Age_LT65, inveGFR*Year1, Age_LT65*Baseline_T2D,
        Age_LT65*Year1, Baseline_T2D*Year1,
        inveGFR*Age_LT65*Year1, Age_LT65*Baseline_T2D*Year1
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "MI") {
      coefficients_used <- coefficients_cost_MI %>%
        pull(coef) %>%
        as.matrix()
      
      value_used <- tibble(
        1, inveGFR, Age_LT65, Baseline_T2D, Year1, 
        inveGFR*Age_LT65, inveGFR*Baseline_T2D, inveGFR*Year1, 
        Age_LT65*Baseline_T2D,
        Age_LT65*Year1, Baseline_T2D*Year1,
        inveGFR*Age_LT65*Baseline_T2D, inveGFR*Age_LT65*Year1
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "Stroke") {
      coefficients_used <- coefficients_cost_Stroke %>%
        pull(coef)
      
      value_used <- tibble(
        1, inveGFR, Age_LT65, Baseline_T2D, Year1, 
        inveGFR*Age_LT65, inveGFR*Year1, 
        Age_LT65*Baseline_T2D,
        Age_LT65*Year1, Baseline_T2D*Year1,
        Age_LT65*Baseline_T2D*Year1
      ) %>%
        as.matrix() %>%
        unname()
    }

    predicted_cost <- 
      (eGFR >= 60) * 
      c_g12_live + 
      
      (eGFR < 60) * 
      infl_to_latest(
        cost = exp(t(coefficients_used) %*% t(value_used)),
        year = 2022,
        to_year = to_year
      )
       
    predicted_cost <- c(predicted_cost)
    
    if (UACR >= 300) {
      predicted_cost <- predicted_cost + c_add_a3_live
    }
    
    return(unname(predicted_cost))
  }
  
  # eGFR_seq <- seq(from = 15, to = 70, by = 1)
  # plot(x = eGFR_seq, y = fn_predict_cost("CKD", eGFR_seq, 360, 1, 1, 1))
  # fn_predict_cost("CKD", 45.3, 360, 1, 1, 1)
  
  
  ####################### Smoothing function for PROB ##########################
  coefficients_prob_HF <- read_csv("data/Model_outputs/coefficients_prob_HF.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_prob_HF <- read_csv("data/Model_outputs/vcov_prob_HF.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_prob_CVD <- read_csv("data/Model_outputs/coefficients_prob_CVD.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_prob_CVD <- read_csv("data/Model_outputs/vcov_prob_CVD.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_prob_ESRD <- read_csv("data/Model_outputs/coefficients_prob_ESRD.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_prob_ESRD <- read_csv("data/Model_outputs/vcov_prob_ESRD.csv") %>%
    rename(parameter = `...1`)
  
  coefficients_prob_Death <- read_csv("data/Model_outputs/coefficients_prob_Death.csv") %>%
    rename(parameter = `...1`, coef = x)
  
  vcov_prob_Death <- read_csv("data/Model_outputs/vcov_prob_Death.csv") %>%
    rename(parameter = `...1`)
  
  fn_predict_prob <- function(event = "HF", 
                              eGFR = c(45,50), 
                              UACR = 150,
                              Diabetes = 1) {
    inveGFR = 1 / eGFR
    lnUACR = log(UACR)
    if (event == "HF") {
      coefficients_used <- coefficients_prob_HF %>%
        pull(coef)
      
      value_used <- tibble(
        1, inveGFR, Diabetes, lnUACR, inveGFR*lnUACR, Diabetes*lnUACR
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "CVD") {
      coefficients_used <- coefficients_prob_CVD %>%
        pull(coef) %>%
        as.matrix()
      
      value_used <- tibble(
        1, inveGFR, Diabetes, lnUACR, inveGFR*lnUACR
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "ESRD") {
      coefficients_used <- coefficients_prob_ESRD %>%
        pull(coef) %>%
        as.matrix()
      
      value_used <- tibble(
        1, inveGFR, Diabetes, lnUACR
      ) %>%
        as.matrix() %>%
        unname()
    }else if (event == "Death") {
      coefficients_used <- coefficients_prob_Death %>%
        pull(coef)
      
      value_used <- tibble(
        1, inveGFR, Diabetes, lnUACR, inveGFR*lnUACR
      ) %>%
        as.matrix() %>%
        unname()
    }
    
    predicted_odds_5y <- 
      exp(t(coefficients_used) %*% t(value_used))
    
    predicted_prob_1y <- 1 - (1 + predicted_odds_5y) ^ (-1/5)
    
    predicted_prob_1y <- c(predicted_prob_1y)

    return(unname(predicted_prob_1y))
  }
  
  # eGFR_seq <- seq(from = 15, to = 90, by = 1)
  # UACR_seq <- seq(from = 15, to = 600, by = 5)
  # plot(x = eGFR_seq, y = fn_predict_prob("HF", eGFR_seq, 60, 1))
  # plot(x = UACR_seq, y = fn_predict_prob("HF", 25, UACR_seq, 1))
  # fn_predict_prob("HF", 60, 60, 1)
  
  
  cat_egfr_names <- c("G1/2", "G3a", "G3b", "G4/5") # >=60 45-59, 30-44, <30
  cat_egfr_6l_names <- c("G1", "G2", "G3a", "G3b", "G4", "G5") # >=90 60-89, 45-59, 30-44, 15-29, <15
  cat_alb_names <- c("A1", "A2", "A3") # <30, 30-299, >=300
  cat_diab_names <- c("With diabetes", "Without diabetes")
  
  n_cat_egfr <- length(cat_egfr_names)
  n_cat_egfr_6l <- length(cat_egfr_6l_names)
  n_cat_alb <- length(cat_alb_names)
  n_cat_diab <- length(cat_diab_names)
  
  cat_egfr_input_names <- c("G1 (≥90)",
                            "G2 (60-89)", 
                            "G3a (45-59)", 
                            "G3b (30-44)",
                            "G4 (15-29)",
                            "G5 (<15)")
  mid_egfr_input <- c(95,
                      75,
                      51.7,
                      37.8,
                      24.6,
                      7.5) 
  
  input_egfr <- mid_egfr_input[cat_egfr_input_names == input_egfrstage]
  
  cat_alb_input_names <- c(
    "A1 (<30)",
    "A2 (30-300)", 
    "A3 (>300)")
  
  mid_alb_input <- c(25.4,
                     90.6,
                     301.5)
  
  input_alb <- mid_alb_input[cat_alb_input_names == input_albstage]
  
  ## payoffs from core model ----
  # Set up the "CKD only" cost matrix ----
  c_mat_noevent <- array(dim = c(n_cat_egfr_6l, n_cat_alb, n_age65, n_fy),
                         dimnames = list(eGFR = cat_egfr_6l_names,
                                         Alb = cat_alb_names,
                                         Age = l_age65,
                                         Year = l_first_year))

  ## For G1 and G2 ----
  c_mat_noevent[1:2,1:2,,] <- c_g12_live
  c_mat_noevent[1:2,3,,] <- c_g12_live + c_add_a3_live
  
  
  ## Function to set up the costs used in the core model for G3a to G5 ----
  source("function/cost_set_up.R")
  
  c_noevent_g3a_g5 <- cost_set_up_g3a_g5(data = cost_ckdcvd,
                                         cat1 = cat_egfr_6l_names[-(1:2)],
                                         cat2 = l_age65,
                                         cat3 = l_first_year,
                                         event = "ckd only")
  c_mat_noevent[-(1:2),1,,] <- c_noevent_g3a_g5[[1]] * (1 - input_dia) + c_noevent_g3a_g5[[2]] * input_dia
                                                   
  ### cost of A2 = cost of A1
  ### additional cost of A3 on top of A1/A2
  c_mat_noevent[-(1:2),2,,] <- c_mat_noevent[-(1:2),1,,]
  c_mat_noevent[-(1:2),3,,] <- c_mat_noevent[-(1:2),1,,] + c_add_a3_live
  
  ## utility
  u_g12_input <- c("Mean" = 0.85, "Lower" = 1, "Upper" = 0.7)
  u_g3a_input <- c("Mean" = 0.8, "Lower" = 1, "Upper" = 0.69)
  u_g3b_input <- c("Mean" = 0.8, "Lower" = 1, "Upper" = 0.68)
  u_g4_input <- c("Mean" = ((423/(423+65))*0.74)+((1-(423/(423+65)))*0.55),
                   "SE" = ((423/(423+65))*(0.85-0.62)/(2 * qnorm(0.975)))+((1-(423/(423+65)))*0.34))
  u_g5_input <- c("Mean" = ((28/(28+75))*0.54)+((1-(28/(28+75)))*0.73),
                  "SE" = ((28/(28+75))*0.36)+((1-(28/(28+75)))*(1-0.72)/(2 * qnorm(0.975))))
  
  disu_a3_input <- c("Mean" = 0.048, "SE" = 0.022)
  
  disu_hf_3m_input <- c("Mean" = 0.105, "SE" = (0.63 - 0.54)/(2 * qnorm(0.975)))
  disu_mi_3m_input <- c("Mean" = 0.12, "SE" = (0.58 - 0.47)/(2 * qnorm(0.975)))
  disu_str_3m_input <- c("Mean" = 0.2, "SE" = (0.63 - 0.46)/(2 * qnorm(0.975)))
  disu_hf_sub_input <- c("Mean" = 0.073, "Lower" = 0.053, "Upper" = 0.094)
  disu_mi_sub_input <- c("Mean" = 0.043, "Lower" = 0.034, "Upper" = 0.051)
  disu_str_sub_input <- c("Mean" = 0.081, "Lower" = 0.066, "Upper" = 0.096)
  
  u_esrd_d_4m_input <-  c("Mean" = (1767*0.75+99*0.44+200*0.78+128*0.6+51*0.69+271*0.54+64*0.53+38*0.54+185*0.74+185*0.58+185*0.67)/(1767+99+200+128+51+271+64+38+185+185+185),
                          "SE" = ((1767*(0.25/sqrt(1767))+
                                     99*(0.32/sqrt(99))+
                                     200*(0.24/sqrt(200))+
                                     128*((0.94-0.55)/3.96)+
                                     51*((0.76-0.63)/3.96)+
                                     271*(0.31/sqrt(271))+
                                     64*(0.34/sqrt(64))+
                                     38*(0.31/sqrt(38))+
                                     185*(0.2/sqrt(185))+
                                     185*(0.26/sqrt(185))+
                                     185*(0.13/sqrt(185))))/(1767+99+200+128+51+271+64+38+185+185+185))
  
  u_esrd_d_sub_input <- c("Mean" = 0.8689*0.62+(1-0.8689)*0.67,
                          "SE" = 0.8689*0.03+(1-0.8689)*0.046)
  u_esrd_t_4m_input <- c("Mean" = (386*0.79+172*0.87+172*0.75+51*0.87+51*0.74+19*0.82+19*0.67+209*0.71+126*0.77+80*0.76)/(386+172+172+51+51+19+19+209+126+80),
                         "SE" = ((386*0.25+172*0.14+172*0.26+51*0.1+51*0.22+19*0.12+19*0.33+209*0.27))/(386+172+172+51+51+19+19+209))
  u_esrd_t_sub_input <- c("Mean" = 0.82, "SE" = 0.02)
  
  # The live value (in use for base or SA)
  u_vectors_input <- list(u_g12_input = u_g12_input,
                          
                          u_g3a_input = u_g3a_input,
                          u_g3b_input = u_g3b_input,
                          u_g4_input = u_g4_input,
                          u_g5_input = u_g5_input,
                          
                          disu_a3_input = disu_a3_input,
                          
                          disu_hf_3m_input = disu_hf_3m_input,
                          disu_mi_3m_input = disu_mi_3m_input,
                          disu_str_3m_input = disu_str_3m_input,
                          disu_hf_sub_input = disu_hf_sub_input,
                          disu_mi_sub_input = disu_mi_sub_input,
                          disu_str_sub_input = disu_str_sub_input,
                          
                          u_esrd_d_4m_input = u_esrd_d_4m_input,
                          
                          u_esrd_d_sub_input = u_esrd_d_sub_input,
                          u_esrd_t_4m_input = u_esrd_t_4m_input,
                          u_esrd_t_sub_input = u_esrd_t_sub_input)
  
  u_vectors_live <- lapply(u_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_beta(vec.mean_se = x, n = 2)[analysis]
    }else if (length(x) == 3){
      x <- SA_beta(vec.mean_l_u = x, n = 2)[analysis]
    }
    return(x)
  }) # beta distribution, n = 2
  
  # Replace "input" with "live" in the names
  names_u_vectors_live <- gsub("input", "live", names(u_vectors_live))
  
  # Assign the new names back to the list
  names(u_vectors_live) <- names_u_vectors_live
  
  list2env(u_vectors_live, .GlobalEnv)
  rm(names_u_vectors_live)

  # Set up the "CKD only" utility matrix ----
  u_mat_noevent <- matrix(data = rep(c(u_g12_live,u_g12_live,u_g3a_live,u_g3b_live,u_g4_live,u_g5_live),3),
                          ncol = n_cat_alb,
                          dimnames = list(eGFR = cat_egfr_6l_names,
                                          Alb = cat_alb_names))
  u_mat_noevent[,"A3"] <- u_mat_noevent[,"A3"] - disu_a3_live
  
  ## treatments
  # statin
  RR_sta_hf_input <- c("Mean" = 0.92, "Lower" = 0.85, "Upper" = 0.99)
  RR_sta_mi_input <- c("Mean" = 0.77^1.09, "Lower" = 0.73, "Upper" = 0.79)
  RR_sta_str_input <- c("Mean" = 0.83^1.09, "Lower" = 0.81, "Upper" = 0.9)
  c_sta_input <- c("Mean" = infl_to_latest(cost=72.1,year=2023), 
                   "SE" = 0.2 * infl_to_latest(cost=72.1,year=2023))
  p_sta_diabetes_input <- c("Mean" = 0.005*(1-input_dia), "SE" = NA)
  p_sta_mae_input <- c("Mean" = 0.047, "SE" = NA)
  p_sta_sae_input <- c("Mean" = 0.00006, "SE" = NA)
  c_sta_diabetes_input <- c("Mean" = infl_to_latest(cost=9600,year=2017), "SE" = NA)
  c_sta_mae_input <- c("Mean" = infl_to_latest(cost=178,year=2019), "SE" = NA)
  c_sta_sae_input <- c("Mean" = infl_to_latest(cost=7033,year=2019), "SE" = NA)
  disu_sta_diabetes_input <- c("Mean" = 0.0351, "SE" = NA)
  disu_sta_mae_input <- c("Mean" = 0.1*(2/365.25), "SE" = NA)
  disu_sta_sae_input <- c("Mean" = 0.2*(2/52), "SE" = NA)
  
  tx_vectors_input <- list(p_sta_diabetes_input = p_sta_diabetes_input,
                           p_sta_mae_input = p_sta_mae_input,
                           p_sta_sae_input = p_sta_sae_input,
                           c_sta_diabetes_input = c_sta_diabetes_input,
                           c_sta_mae_input = c_sta_mae_input,
                           c_sta_sae_input = c_sta_sae_input,
                           disu_sta_diabetes_input = disu_sta_diabetes_input,
                           disu_sta_mae_input = disu_sta_mae_input,
                           disu_sta_sae_input = disu_sta_sae_input)
  
  tx_vectors_input <- lapply(tx_vectors_input, function(x) {
    x["SE"] <- abs(x["Mean"] / 10) # Assuming the SE is 0.1 times of the Mean
    return(x)
  })
  list2env(tx_vectors_input, .GlobalEnv)
  rm(tx_vectors_input)
  
  c_sta_live <- SA_gamma(vec.mean_se = c_sta_input, n = 6)[analysis] # gamma distribution, n = 6
  RR_sta_hf_live = SA_lognorm(vec.mean_l_u = RR_sta_hf_input, n = 3)[analysis] # log-normal distribution, n = 3
  RR_sta_mi_live = SA_lognorm(vec.mean_l_u = RR_sta_mi_input, n = 4)[analysis] # log-normal distribution, n = 4
  RR_sta_str_live = SA_lognorm(vec.mean_l_u = RR_sta_str_input, n = 5)[analysis] # log-normal distribution, n = 5
  
  tx_b_vectors_input <- list(p_sta_diabetes_input = p_sta_diabetes_input,
                             p_sta_mae_input = p_sta_mae_input,
                             p_sta_sae_input = p_sta_sae_input) # beta distribution
  
  tx_g_vectors_input <- list(c_sta_diabetes_input = c_sta_diabetes_input,
                             c_sta_mae_input = c_sta_mae_input,
                             c_sta_sae_input = c_sta_sae_input) # gamma distribution
  
  tx_n_vectors_input <- list(disu_sta_diabetes_input = disu_sta_diabetes_input,
                             disu_sta_mae_input = disu_sta_mae_input,
                             disu_sta_sae_input = disu_sta_sae_input) # normal distribution
  
  tx_b_vectors_live <- lapply(tx_b_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_beta(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_beta(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_g_vectors_live <- lapply(tx_g_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_gamma(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_gamma(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_n_vectors_live <- lapply(tx_n_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_norm(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_norm(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  # Replace "input" with "live" in the names
  names_tx_b_vectors_live <- gsub("input", "live", names(tx_b_vectors_live))
  names_tx_g_vectors_live <- gsub("input", "live", names(tx_g_vectors_live))
  names_tx_n_vectors_live <- gsub("input", "live", names(tx_n_vectors_live))
  
  # Assign the new names back to the list
  names(tx_b_vectors_live) <- names_tx_b_vectors_live
  names(tx_g_vectors_live) <- names_tx_g_vectors_live
  names(tx_n_vectors_live) <- names_tx_n_vectors_live
  
  list2env(tx_b_vectors_live, .GlobalEnv)
  list2env(tx_g_vectors_live, .GlobalEnv)
  list2env(tx_n_vectors_live, .GlobalEnv)
  
  rm(names_tx_b_vectors_live,
     names_tx_g_vectors_live,
     names_tx_n_vectors_live)
  
  
  usage_sta_base <- 0.63
  
  # blood pressure
  SBP_reduct_input <- c("Mean" = 9.1, "SE" = 9.1 * 0.2)
  c_BP_input <- c("Mean" = 54.6, "SE" = 54.6 * 0.2)
  prop_unctrl <-0.33 # fixed
  RR_BP_hf_input <- c("Mean" = 0.72, "Lower" = 0.67, "Upper" = 0.78)
  RR_BP_mi_input <- c("Mean" = 0.83, "Lower" = 0.78, "Upper" = 0.88)
  RR_BP_str_input <- c("Mean" = 0.73, "Lower" = 0.68, "Upper" = 0.77)
  p_BP_itae_input <- c("Mean" = (0.005+0.009+0.016+0.023+0.03)/5, "SE" = (0.002+0.005+0.008+0.007+0.014)/5)
  p_BP_sae_input <- c("Mean" = (0.009+0.012)/2, "SE" = (0.002+0.005+0.008+0.007+0.014)/5)
  c_BP_itae_input <- c("Mean" = infl_to_latest(cost=110.28,year=2019), "SE" = infl_to_latest(cost=110.28,year=2019) * 0.1)
  c_BP_sae_input <- c("Mean" = infl_to_latest(cost=147.76,year=2019), "SE" = infl_to_latest(cost=147.76,year=2019) * 0.1)
  disu_BP_itae_input <- c("Mean" = 0.1*(2/365.25), "SE" = 0.1*(2/365.25) * 0.1)
  disu_BP_sae_input <- c("Mean" = 0.2*(2/52), "SE" = 0.2*(2/52) * 0.1)
  
  c_BP_live <- SA_gamma(vec.mean_se = c_BP_input, n = 6)[analysis] # gamma distribution, n = 6
  RR_BP_hf_live = SA_lognorm(vec.mean_l_u = RR_BP_hf_input, n = 3)[analysis] # log-normal distribution, n = 3
  RR_BP_mi_live = SA_lognorm(vec.mean_l_u = RR_BP_mi_input, n = 4)[analysis] # log-normal distribution, n = 4
  RR_BP_str_live = SA_lognorm(vec.mean_l_u = RR_BP_str_input, n = 5)[analysis] # log-normal distribution, n = 5
  
  tx_b_vectors_input <- list(p_BP_itae_input = p_BP_itae_input,
                             p_BP_sae_input = p_BP_sae_input) # beta distribution
  
  tx_g_vectors_input <- list(c_BP_itae_input = c_BP_itae_input,
                             c_BP_sae_input = c_BP_sae_input) # gamma distribution
  
  tx_n_vectors_input <- list(SBP_reduct_input = SBP_reduct_input,
                             disu_BP_itae_input = disu_BP_itae_input,
                             disu_BP_sae_input = disu_BP_sae_input) # normal distribution
  
  tx_b_vectors_live <- lapply(tx_b_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_beta(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_beta(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_g_vectors_live <- lapply(tx_g_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_gamma(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_gamma(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_n_vectors_live <- lapply(tx_n_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_norm(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_norm(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  # Replace "input" with "live" in the names
  names_tx_b_vectors_live <- gsub("input", "live", names(tx_b_vectors_live))
  names_tx_g_vectors_live <- gsub("input", "live", names(tx_g_vectors_live))
  names_tx_n_vectors_live <- gsub("input", "live", names(tx_n_vectors_live))
  
  # Assign the new names back to the list
  names(tx_b_vectors_live) <- names_tx_b_vectors_live
  names(tx_g_vectors_live) <- names_tx_g_vectors_live
  names(tx_n_vectors_live) <- names_tx_n_vectors_live
  
  list2env(tx_b_vectors_live, .GlobalEnv)
  list2env(tx_g_vectors_live, .GlobalEnv)
  list2env(tx_n_vectors_live, .GlobalEnv)
  
  rm(names_tx_b_vectors_live,
     names_tx_g_vectors_live,
     names_tx_n_vectors_live)
  
  usage_BP_base <- 0.98
  
  # aspirin
  RR_asp_hf_input <- c("Mean" = 1, "SE" = 0.01)
  RR_asp_mi_input <- c("Mean" = 0.86, "Lower" = 0.8, "Upper" = 1)
  RR_asp_str_input <- c("Mean" = 0.94, "Lower" = 0.84, "Upper" = 1.06)
  c_asp_input <- c("Mean" = 15.1, "SE" = 15.1 * 0.2)
  p_asp_mb_input <- c("Mean" = 1.39/1000, "Lower" = 0.7/1000, "Upper" = 2.28/1000)
  c_asp_mb_input <- c("Mean" = infl_to_latest(cost=13093,year=2018), "SE" = infl_to_latest(cost=13093,year=2018) * 0.1)
  disu_asp_mb_input <- c("Mean" = 0.0445, "Lower" = 0.016, "Upper" = 0.073)
  
  c_asp_live <- SA_gamma(vec.mean_se = c_asp_input, n = 6)[analysis] # gamma distribution, n = 6
  RR_asp_hf_live = SA_lognorm(vec.mean_se = RR_asp_hf_input, n = 3)[analysis] # log-normal distribution, n = 3
  RR_asp_mi_live = SA_lognorm(vec.mean_l_u = RR_asp_mi_input, n = 4)[analysis] # log-normal distribution, n = 4
  RR_asp_str_live = SA_lognorm(vec.mean_l_u = RR_asp_str_input, n = 5)[analysis] # log-normal distribution, n = 5
  
  tx_b_vectors_input <- list(p_asp_mb_input = p_asp_mb_input) # beta distribution
  
  tx_g_vectors_input <- list(c_asp_mb_input = c_asp_mb_input) # gamma distribution
  
  tx_n_vectors_input <- list(disu_asp_mb_input = disu_asp_mb_input) # normal distribution
  
  tx_b_vectors_live <- lapply(tx_b_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_beta(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_beta(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_g_vectors_live <- lapply(tx_g_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_gamma(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_gamma(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  tx_n_vectors_live <- lapply(tx_n_vectors_input, function(x) {
    if (length(x) == 2) {
      x <- SA_norm(vec.mean_se = x)[analysis]
    }else if (length(x) == 3){
      x <- SA_norm(vec.mean_l_u = x)[analysis]
    }
    return(x)
  })
  
  # Replace "input" with "live" in the names
  names_tx_b_vectors_live <- gsub("input", "live", names(tx_b_vectors_live))
  names_tx_g_vectors_live <- gsub("input", "live", names(tx_g_vectors_live))
  names_tx_n_vectors_live <- gsub("input", "live", names(tx_n_vectors_live))
  
  # Assign the new names back to the list
  names(tx_b_vectors_live) <- names_tx_b_vectors_live
  names(tx_g_vectors_live) <- names_tx_g_vectors_live
  names(tx_n_vectors_live) <- names_tx_n_vectors_live
  
  list2env(tx_b_vectors_live, .GlobalEnv)
  list2env(tx_g_vectors_live, .GlobalEnv)
  list2env(tx_n_vectors_live, .GlobalEnv)
  
  rm(names_tx_b_vectors_live,
     names_tx_g_vectors_live,
     names_tx_n_vectors_live)
  
  usage_asp_base <- 870/2578
  
  # SGLT2 inhibitor
  egfr_sglt_input <- c("Mean" = 1.34991272358995, "Lower" = 1.01175, "Upper" = 2.5965)
  effect_year_egfr <- 3
  UAR_sglt <- 0
  effect_year_UAR <- 3
  c_sglt_input <- c("Mean" = 3263, "SE" = 3263 * 0.2)
  RR_sglt_hf_input <- c("Mean" = 0.65, "Lower" = 0.55, "Upper" = 0.71)
  RR_sglt_mi_input <- c("Mean" = 0.84, "Lower" = 0.74, "Upper" = 0.95)
  RR_sglt_str_input <- c("Mean" = 0.84, "Lower" = 0.74, "Upper" = 0.95)
  
  egfr_sglt_live <- SA_gamma(vec.mean_l_u = egfr_sglt_input, n = 10)[analysis] # normal distribution, n = 10
  c_sglt_live <- SA_gamma(vec.mean_se = c_sglt_input, n = 11)[analysis] # gamma distribution, n = 11
  RR_sglt_hf_live = SA_lognorm(vec.mean_l_u = RR_sglt_hf_input, n = 7)[analysis] # log-normal distribution, n = 7
  RR_sglt_mi_live = SA_lognorm(vec.mean_l_u = RR_sglt_mi_input, n = 8)[analysis] # log-normal distribution, n = 8
  RR_sglt_str_live = SA_lognorm(vec.mean_l_u = RR_sglt_str_input, n = 9)[analysis] # log-normal distribution, n = 9
  
  usage_sglt_base <- 0
  
  # if the risk of these events is impacted by treated eGFR
  # set the RRs to be 1 to avoid double counting
  if(!effect_sglt_RR) { 
    RR_sglt_hf_live <- 1
    RR_sglt_mi_live <- 1
    RR_sglt_str_live <- 1
  }
  
  ### Scenario ----
  # Default: Usual care vs CV risk factor management
  usage_sta_target <- c(1, min(1.5 * usage_sta_base, 1))
  optm_sta <- c(FALSE, TRUE)
  
  usage_BP_target <- c(1, min(1.5 * usage_BP_base, 1))
  optm_BP <- c(FALSE, TRUE)
  
  usage_asp_target <- c(1, min(1.5 * usage_asp_base, 1))
  optm_asp <- c(FALSE, FALSE)
  
  usage_sglt_target <- c(0, 0)
  optm_sglt <- c(FALSE, FALSE)
  
  delayTreat <- c(1, 1) # time to treatment from diagnosis
  
  dict <- list(
    "Usual care" = list(
      usage_sta_target = usage_sta_base, 
      optm_sta = FALSE, 
      usage_BP_target = usage_BP_base, 
      optm_BP = FALSE,
      usage_asp_target = usage_asp_base,
      optm_asp = FALSE, 
      usage_sglt_target = 0, 
      optm_sglt = FALSE,
      delayTreat = 1
    ), 
    "CV risk factor management" = list(
      usage_sta_target = min(1.5 * usage_sta_base, 1), 
      optm_sta = TRUE, 
      usage_BP_target = min(1.5 * usage_BP_base, 1), 
      optm_BP = TRUE, 
      usage_asp_target = min(1.5 * usage_asp_base, 1),
      optm_asp = FALSE, 
      usage_sglt_target = 0, 
      optm_sglt = FALSE,
      delayTreat = 1
    ),
    "CV risk factor management for all eligible" = list(
      usage_sta_target = 1, 
      optm_sta = TRUE, 
      usage_BP_target = 1, 
      optm_BP = TRUE, 
      usage_asp_target = 1,
      optm_asp = FALSE, 
      usage_sglt_target = 0, 
      optm_sglt = FALSE,
      delayTreat = 1
    ),
    "CV risk factor management for all eligible and 50% CKD management" = list(
      usage_sta_target = 1, 
      optm_sta = TRUE, 
      usage_BP_target = 1, 
      optm_BP = TRUE, 
      usage_asp_target = 1,
      optm_asp = FALSE, 
      usage_sglt_target = 0.5, 
      optm_sglt = TRUE,
      delayTreat = 1
    ),
    "CV risk factor management for all eligible and full CKD management" = list(
      usage_sta_target = 1, 
      optm_sta = TRUE, 
      usage_BP_target = 1, 
      optm_BP = TRUE, 
      usage_asp_target = 1,
      optm_asp = FALSE, 
      usage_sglt_target = 1, 
      optm_sglt = TRUE,
      delayTreat = 1
    ),
    "CV risk factor management for all eligible and full CKD management, no delay" = list(
      usage_sta_target = 1, 
      optm_sta = TRUE, 
      usage_BP_target = 1, 
      optm_BP = TRUE, 
      usage_asp_target = 1,
      optm_asp = FALSE, 
      usage_sglt_target = 1, 
      optm_sglt = TRUE,
      delayTreat = 0
    )
  )
  
  # Assign the values based on the intervention and comparator
  if (intervention %in% names(dict)) {
    usage_sta_target[2] <- dict[[intervention]]$usage_sta_target
    optm_sta[2] <- dict[[intervention]]$optm_sta
    usage_BP_target[2] <- dict[[intervention]]$usage_BP_target
    optm_BP[2] <- dict[[intervention]]$optm_BP
    usage_asp_target[2] = dict[[intervention]]$usage_asp_target
    optm_asp[2] <- dict[[intervention]]$optm_asp
    usage_sglt_target[2] <- dict[[intervention]]$usage_sglt_target
    optm_sglt[2] <- dict[[intervention]]$optm_sglt
    delayTreat[2] <- dict[[intervention]]$delayTreat
  }
  
  if (comparator %in% names(dict)) {
    usage_sta_target[1] <- dict[[comparator]]$usage_sta_target
    optm_sta[1] <- dict[[comparator]]$optm_sta
    usage_BP_target[1] <- dict[[comparator]]$usage_BP_target
    optm_BP[1] <- dict[[comparator]]$optm_BP
    usage_asp_target[1] = dict[[comparator]]$usage_asp_target
    optm_asp[1] <- dict[[comparator]]$optm_asp
    usage_sglt_target[1] <- dict[[comparator]]$usage_sglt_target
    optm_sglt[1] <- dict[[comparator]]$optm_sglt
    delayTreat[1] <- dict[[comparator]]$delayTreat
  }
  
  ### treatment effects ----
  
  # statin
  inc_sta <- optm_sta * (usage_sta_target - usage_sta_base)
  
  RR_sta_hf_pop <- RR_sta_hf_live ^ inc_sta
  RR_sta_mi_pop <- RR_sta_mi_live ^ inc_sta
  RR_sta_str_pop <- RR_sta_str_live ^ inc_sta
  
  c_sta_pop <- c_sta_live * inc_sta
  c_sta_diabetes_pop <- c_sta_diabetes_live * p_sta_diabetes_live * inc_sta
  c_sta_mae_pop <- c_sta_mae_live * p_sta_mae_live * inc_sta
  c_sta_sae_pop <- c_sta_sae_live * p_sta_sae_live * inc_sta
  
  disu_sta_diabetes_pop <- disu_sta_diabetes_live * p_sta_diabetes_live * inc_sta
  disu_sta_mae_pop <- disu_sta_mae_live * p_sta_mae_live * inc_sta
  disu_sta_sae_pop <- disu_sta_sae_live * p_sta_sae_live * inc_sta
  
  # BP
  inc_BP <- optm_BP * input_hyp * ((usage_BP_target - usage_BP_base) + usage_BP_base * prop_unctrl)
  
  RR_BP_hf_pop <- (RR_BP_hf_live ^ (SBP_reduct_live/10)) ^ inc_BP
  RR_BP_mi_pop <- (RR_BP_mi_live ^ (SBP_reduct_live/10)) ^ inc_BP
  RR_BP_str_pop <- (RR_BP_str_live ^ (SBP_reduct_live/10)) ^ inc_BP
  
  c_BP_pop <- c_BP_live * inc_BP
  c_BP_itae_pop <- c_BP_itae_live * p_BP_itae_live * inc_BP
  c_BP_sae_pop <- c_BP_sae_live * p_BP_sae_live * inc_BP
  
  disu_BP_itae_pop <- disu_BP_itae_live * p_BP_itae_live * inc_BP
  disu_BP_sae_pop <- disu_BP_sae_live * p_BP_sae_live * inc_BP
  
  # aspirin
  inc_asp <- optm_asp * (usage_asp_target - usage_asp_base)
  
  RR_asp_hf_pop <- RR_asp_hf_live ^ inc_asp
  RR_asp_mi_pop <- RR_asp_mi_live ^ inc_asp
  RR_asp_str_pop <- RR_asp_str_live ^ inc_asp
  
  c_asp_pop <- c_asp_live * inc_asp
  c_asp_mb_pop <- c_asp_mb_live * p_asp_mb_live * inc_asp
  
  disu_asp_mb_pop <- disu_asp_mb_live * p_asp_mb_live * inc_asp
  
  # SGLT2i
  inc_sglt <- optm_sglt * (usage_sglt_target - usage_sglt_base)
  
  RR_sglt_hf_pop <- RR_sglt_hf_live ^ inc_sglt
  RR_sglt_mi_pop <- RR_sglt_mi_live ^ inc_sglt
  RR_sglt_str_pop <- RR_sglt_str_live ^ inc_sglt
  
  egfr_sglt_pop <- egfr_sglt_live * inc_sglt
  effect_year_egfr_pop <- effect_year_egfr + delayTreat
  UAR_sglt_pop <- UAR_sglt * inc_sglt
  effect_year_UAR_pop <- effect_year_UAR + delayTreat
  
  c_sglt_pop <- c_sglt_live * inc_sglt
  
  ### core model ----
  # "The annual risk of each event depends on the eGFR and albuminuria levels during each model cycle, as well as the baseline type 2 diabetes status"

  # "As eGFR declines and albuminuria increases over time, the risk of each clinical event increases"
  # "the rate of decline (of eGFR is) dependent on the current eGFR and current albuminuria in the cohort"
  # table 1
  changemat_egfr_plus_input <- array(data = c(-0.8, -0.3, -0.3, -0.1, 
                                              -2.2, -2.1, -1.5, -1.1, 
                                              -4.6, -4.6, -4.5, -3.6, 
                                              -0.1, -0.2, -0.2, -0.2, 
                                              -1.0, -1.5, -1.4, -1.2, 
                                              -3.1, -4.0, -3.2, -2.8),
                                     dim = c(n_cat_egfr, n_cat_alb, n_cat_diab),
                                     dimnames = list(eGFR = cat_egfr_names,
                                                     Alb = cat_alb_names,
                                                     Diab = cat_diab_names))
  apply_SA_norm <- function(x) {
    # Normal distribution, SE assumed to be 10% of mean 
    se_val <- 0.1 * abs(x)
    
    # Use the SA_norm function to generate the sensitivity analysis values
    sa_values <- SA_norm(vec.mean_se = c(x, se_val))
    
    return(sa_values[analysis])
  }
  
  # Apply the function to each element in changemat_egfr_plus_input
  changemat_egfr_plus_live <- apply(changemat_egfr_plus_input, 1:length(dim(changemat_egfr_plus_input)), apply_SA_norm)
  
  # weighted average over diabetes status
  changemat_egfr_plus <- changemat_egfr_plus_live[,,"With diabetes"]*input_dia + changemat_egfr_plus_live[,,"Without diabetes"]*(1-input_dia) 
  
  # "In each cycle, the mean albuminuria level in the cohort increases based on the previous level of albuminuria"
  # table 2
  changemat_alb_mult <- array(data = c(1.05, 1.10, 0.94),
                              dim = n_cat_alb,
                              dimnames = list(Alb = cat_alb_names))
  
  ## Mortality - Event ----
  # "the proportion of patients who die after their event"
  # table 7
  mt_hf_input <- c("Mean" = 0.043, "SE" = 0.043 * 0.1)
  mt_mi_input <- c("Mean" = 0.061, "SE" = 0.061 * 0.1)
  mt_str_input <- c("Mean" = 0.1, "SE" = 0.1 * 0.1)
  mt_esrd_input <- c("Mean" = 0.034, "SE" = 0.034 * 0.1)
  
  mt_hf_live <- SA_beta(vec.mean_se = mt_hf_input)[analysis]
  mt_mi_live <- SA_beta(vec.mean_se = mt_mi_input)[analysis]
  mt_str_live <- SA_beta(vec.mean_se = mt_str_input)[analysis]
  mt_esrd_live <- SA_beta(vec.mean_se = mt_esrd_input)[analysis]
  
  
  ## Heart failure ----
  # table 3
  # mean
  prob_mat_hf5y_input <- array(data = c(0.06, 0.07, 0.07, 0.24, 
                                        0.10, 0.13, 0.13, 0.19, 
                                        0.17, 0.14, 0.19, 0.24, 
                                        0.01, 0.02, 0.03, 0.11, 
                                        0.06, 0.03, 0.13, 0.09, 
                                        0.14, 0.11, 0.09, 0.15),
                               dim = c(n_cat_egfr, n_cat_alb, n_cat_diab),
                               dimnames = list(eGFR = cat_egfr_names,
                                               Alb = cat_alb_names,
                                               Diab = cat_diab_names))
  # lower bound
  prob_mat_hf5y_lb <- array(data = c(0.02, 0.04, 0.04, 0.14, 
                                     0.03, 0.08, 0.09, 0.12, 
                                     0.07, 0.09, 0.15, 0.18, 
                                     0.00, 0.01, 0.02, 0.05, 
                                     0.02, 0.01, 0.08, 0.04, 
                                     0.04, 0.05, 0.05, 0.09),
                            dim = dim(prob_mat_hf5y_input),
                            dimnames = dimnames(prob_mat_hf5y_input))
  # upper bound
  prob_mat_hf5y_ub <- array(data = c(0.12, 0.11, 0.11, 0.35, 
                                     0.22, 0.18, 0.18, 0.28, 
                                     0.29, 0.21, 0.23, 0.29, 
                                     0.03, 0.04, 0.06, 0.19, 
                                     0.14, 0.07, 0.18, 0.15, 
                                     0.28, 0.19, 0.14, 0.22),
                            dim = dim(prob_mat_hf5y_input),
                            dimnames = dimnames(prob_mat_hf5y_input))
  
  apply_SA_beta <- function(mean, lb, ub) {
    # Use the SA_beta function to generate the sensitivity analysis values
    sa_values <- SA_beta(num.mean = mean, num.lower = lb, num.upper = ub)
    
    return(sa_values[analysis])
  }
  
  # Apply the function to each corresponding set of elements in the arrays
  prob_mat_hf5y_live <- array(mapply(apply_SA_beta, prob_mat_hf5y_input, prob_mat_hf5y_lb, prob_mat_hf5y_ub), 
                              dim = dim(prob_mat_hf5y_input),
                              dimnames = dimnames(prob_mat_hf5y_input))
  
  # weighted average over diabetes status
  prob_mat_hf5y <- prob_mat_hf5y_live[,,"With diabetes"]*input_dia + prob_mat_hf5y_live[,,"Without diabetes"]*(1-input_dia) 
  
  prob_mat_nfatal_hf1y <- (1-mt_hf_live)*(1-exp(-(-log(1 - prob_mat_hf5y))/5)) #constant rate assumed
  prob_mat_fatal_hf1y <- mt_hf_live*(1-exp(-(-log(1 - prob_mat_hf5y))/5)) #constant rate assumed
  
  ## CVD ----
  # table 4
  # mean
  prob_mat_cvd5y_input <- array(data = c(0.06, 0.14, 0.16, 0.36, 
                                         0.15, 0.18, 0.22, 0.31, 
                                         0.31, 0.30, 0.35, 0.43, 
                                         0.04, 0.06, 0.10, 0.19, 
                                         0.08, 0.10, 0.22, 0.14, 
                                         0.20, 0.15, 0.14, 0.21),
                                dim = c(n_cat_egfr, n_cat_alb, n_cat_diab),
                                dimnames = list(eGFR = cat_egfr_names,
                                                Alb = cat_alb_names,
                                                Diab = cat_diab_names))
  # lower bound
  prob_mat_cvd5y_lb <- array(data = c(0.02, 0.09, 0.11, 0.24, 
                                      0.06, 0.13, 0.16, 0.21, 
                                      0.18, 0.23, 0.30, 0.36, 
                                      0.02, 0.04, 0.07, 0.11, 
                                      0.03, 0.06, 0.16, 0.08, 
                                      0.08, 0.08, 0.09, 0.14),
                             dim = dim(prob_mat_cvd5y_input),
                             dimnames = dimnames(prob_mat_cvd5y_input))
  # upper bound
  prob_mat_cvd5y_ub <- array(data = c(0.12, 0.19, 0.22, 0.47, 
                                      0.28, 0.25, 0.28, 0.41, 
                                      0.45, 0.38, 0.40, 0.49, 
                                      0.07, 0.09, 0.14, 0.28, 
                                      0.16, 0.16, 0.29, 0.21, 
                                      0.36, 0.23, 0.20, 0.29),
                             dim = dim(prob_mat_cvd5y_input),
                             dimnames = dimnames(prob_mat_cvd5y_input))
  
  
  # Apply the function to each corresponding set of elements in the arrays
  prob_mat_cvd5y_live <- array(mapply(apply_SA_beta, prob_mat_cvd5y_input, prob_mat_cvd5y_lb, prob_mat_cvd5y_ub), 
                               dim = dim(prob_mat_cvd5y_input),
                               dimnames = dimnames(prob_mat_cvd5y_input))
  
  prob_mat_cvd5y <- prob_mat_cvd5y_live[,,"With diabetes"]*input_dia + prob_mat_cvd5y_live[,,"Without diabetes"]*(1-input_dia) 
  
  prob_mat_cvd1y <- 1-exp(-(-log(1 - prob_mat_cvd5y))/5) #constant rate assumed
  
  # "For CVD, 65% of cardiovascular events are assumed to be MI and 35% are stroke"
  # table 6
  prop_mi_input <- c("Mean" = 389/(389+211), "SE" = sqrt(389*211/(389+211)^3))
  prop_mi_live <- SA_beta(vec.mean_se = prop_mi_input)[analysis]
  prop_str_live <- 1 - prop_mi_live
  
  prob_mat_nfatal_mi1y <- prob_mat_cvd1y*prop_mi_live*(1-mt_mi_live)
  prob_mat_fatal_mi1y <- prob_mat_cvd1y*prop_mi_live*mt_mi_live
  prob_mat_nfatal_str1y <- prob_mat_cvd1y*prop_str_live*(1-mt_str_live)
  prob_mat_fatal_str1y <- prob_mat_cvd1y*prop_str_live*mt_str_live
  
  ## ESRD ----
  # table 5
  # mean
  prob_mat_esrd5y_input <- array(data = c(0.01, 0.02, 0.03, 0.14, 
                                          0.08, 0.05, 0.13, 0.38, 
                                          0.09, 0.26, 0.46, 0.73, 
                                          0.01, 0.01, 0.02, 0.13, 
                                          0.01, 0.04, 0.13, 0.26, 
                                          0.06, 0.21, 0.34, 0.65),
                                 dim = c(n_cat_egfr, n_cat_alb, n_cat_diab),
                                 dimnames = list(eGFR = cat_egfr_names,
                                                 Alb = cat_alb_names,
                                                 Diab = cat_diab_names))
  # lower bound
  prob_mat_esrd5y_lb <- array(data = c(0.00, 0.00, 0.01, 0.07, 
                                       0.02, 0.02, 0.09, 0.28, 
                                       0.03, 0.19, 0.41, 0.67, 
                                       0.00, 0.00, 0.01, 0.07, 
                                       0.00, 0.02, 0.08, 0.19, 
                                       0.01, 0.13, 0.27, 0.55),
                              dim = dim(prob_mat_esrd5y_input),
                              dimnames = dimnames(prob_mat_esrd5y_input))
  # upper bound
  prob_mat_esrd5y_ub <- array(data = c(0.05, 0.04, 0.06, 0.24, 
                                       0.19, 0.10, 0.18, 0.47, 
                                       0.20, 0.33, 0.52, 0.78, 
                                       0.06, 0.02, 0.04, 0.22, 
                                       0.07, 0.08, 0.18, 0.35, 
                                       0.18, 0.30, 0.41, 0.73),
                              dim = dim(prob_mat_esrd5y_input),
                              dimnames = dimnames(prob_mat_esrd5y_input))
  
  
  # Apply the function to each corresponding set of elements in the arrays
  prob_mat_esrd5y_live <- array(mapply(apply_SA_beta, prob_mat_esrd5y_input, prob_mat_esrd5y_lb, prob_mat_esrd5y_ub), 
                                dim = dim(prob_mat_esrd5y_input),
                                dimnames = dimnames(prob_mat_esrd5y_input))
  
  prob_mat_esrd5y <- prob_mat_esrd5y_live[,,"With diabetes"]*input_dia + prob_mat_esrd5y_live[,,"Without diabetes"]*(1-input_dia) 
  
  prob_mat_esrd1y <- 1-exp(-(-log(1 - prob_mat_esrd5y))/5) #constant rate assumed
  
  # "For those reaching ESRD, it is assumed that 71% receive dialysis or supportive care (no treatment) and the remaining 29% receive a kidney graft"
  # table 6
  prop_transp_input <- c("Mean" = 0.29, "SE" = sqrt(0.29*(1-0.29)/(786000)))
  prop_transp_live <- SA_beta(vec.mean_se = prop_transp_input)[analysis]
  prop_dial_live <- 1 - prop_transp_live
  
  prob_mat_nfatal_esrd_d1y <- prob_mat_esrd1y*prop_dial_live*(1-mt_esrd_live)
  prob_mat_nfatal_esrd_t1y <- prob_mat_esrd1y*prop_transp_live*(1-mt_esrd_live)
  prob_mat_fatal_esrd1y <- prob_mat_esrd1y*mt_esrd_live
  
  ## Mortality - Other causes
  # US 2023 life tables, Verified by Svitlana Usachova on 09/03/26. https://www.cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-06.pdf
  LT <- read_csv("data/US_2023_Life_tables.csv",
                 col_types = cols(Age = col_integer(),
                                  Males = col_double(),
                                  Females = col_double(),
                                  Both = col_double()))
  LT$Overall <- LT$Females*input_sex+LT$Males*(1-input_sex)
  LT$age_range <- cut(LT$Age, breaks = c(15,20,seq(25, 85, by = 10),Inf), right = FALSE)
  
  # Cause specific mortality, US life tables 2023, Verified by Svitlana Usachova on 09/03/26. https://www.cdc.gov/nchs/data/nvsr/nvsr74/NVSR74-10.pdf
  CSM <- read_csv("data/Cause_specific_mortality_US_2023.csv", 
                  col_types = cols(age_range = col_character(), 
                                   female = col_factor(levels = c("m", "f")), 
                                   all_cause = col_double(), 
                                   disease_of_heart = col_double(), 
                                   diabetes = col_double(), 
                                   cerebrovascular_disease = col_double(), 
                                   nephritis_nephrosis = col_double()))
  
  CSM$prop_excl_dia <- (CSM$disease_of_heart+CSM$cerebrovascular_disease+CSM$nephritis_nephrosis) / CSM$all_cause
  
  CSM_bothgender <- data.frame(age_range = character(length(CSM$age_range)/2),
                               prop_excl_dia = numeric(length(CSM$age_range)/2))
  CSM_bothgender$age_range <- CSM$age_range[CSM$female=="f"]
  CSM_bothgender$prop_excl_dia <- CSM$prop_excl_dia[CSM$female=="f"]*input_sex+CSM$prop_excl_dia[CSM$female=="m"]*(1-input_sex)
  
  merged_LT <- merge(LT, CSM_bothgender, by.x = "age_range", by.y = "age_range")
  merged_LT$other_cause <- merged_LT$Overall*(1-merged_LT$prop_excl_dia)
  merged_LT$other_cause[merged_LT$Overall==1] <- 1
  
  event_names <- c("hf", "fatalhf", "mi", "fatalmi", "str", "fatalstr", "esrd_d" , "esrd_t", "fatalesrd", "otherdeath")
  prob_names <- c("prob_hf", "prob_fatalhf", "prob_mi", "prob_fatalmi", "prob_str", "prob_fatalstr", "prob_esrd_d" , "prob_esrd_t", "prob_fatalesrd", "prob_otherdeath")
  
  ## transition probabilities
  tp <- list()
  
  ## Markov trace
  trace <- list()
  
  ## loop for two arms
  for (i in 1:2) {
    tp[[i]] <- data.frame(year = numeric(n_cycles),
                          age = numeric(n_cycles),
                          eGFR0 = numeric(n_cycles), # eGFR if untreated 
                          eGFR = numeric(n_cycles),
                          alb = numeric(n_cycles),
                          prob_hf = numeric(n_cycles),
                          prob_fatalhf = numeric(n_cycles),
                          prob_mi = numeric(n_cycles),
                          prob_fatalmi = numeric(n_cycles),
                          prob_str = numeric(n_cycles),
                          prob_fatalstr = numeric(n_cycles),
                          prob_esrd_d = numeric(n_cycles),
                          prob_esrd_t = numeric(n_cycles),
                          prob_fatalesrd = numeric(n_cycles),
                          prob_otherdeath = numeric(n_cycles)
    )
    tp[[i]]$year[1] <- 0
    tp[[i]]$age[1] <- input_age
    tp[[i]]$eGFR0[1] <- input_egfr
    tp[[i]]$eGFR[1] <- input_egfr
    tp[[i]]$alb[1] <- input_alb
    
    for (j in 2:n_cycles) {
      tp[[i]]$year[j] <- j-1
      tp[[i]]$age[j] <- input_age+j-1
      
      egfr_grp <- cut(tp[[i]]$eGFR[j-1], breaks = c(0,30,45,60,Inf),right=FALSE,labels=c("G4/5", "G3b", "G3a", "G1/2"))
      egfr_grp0 <- cut(tp[[i]]$eGFR0[j-1], breaks = c(0,30,45,60,Inf),right=FALSE,labels=c("G4/5", "G3b", "G3a", "G1/2"))
      alb_grp <- cut(tp[[i]]$alb[j-1], breaks = c(0,30,300,Inf),right=FALSE,labels=c("A1", "A2", "A3"))
      
      egfr_value <- tp[[i]]$eGFR[j-1]
      egfr0_value <- tp[[i]]$eGFR0[j-1]
      alb_value <- tp[[i]]$alb[j-1]
      
      if (tp[[i]]$eGFR[j-1]<15) {
        tp[[i]]$eGFR[j] <- tp[[i]]$eGFR[j-1]
      }else{
        tp[[i]]$eGFR[j] <- tp[[i]]$eGFR[j-1] + 
          changemat_egfr_plus[as.character(egfr_grp),as.character(alb_grp)] + 
          egfr_sglt_pop[i] * (tp[[i]]$year[j]-delayTreat[i]) * (tp[[i]]$year[j]<=effect_year_egfr_pop[i]) * (tp[[i]]$year[j]>delayTreat[i])
      }
      
      if(tp[[i]]$eGFR0[j-1]<15) {
        tp[[i]]$eGFR0[j] <- tp[[i]]$eGFR0[j-1]
      }else{
        tp[[i]]$eGFR0[j] <- tp[[i]]$eGFR0[j-1] + 
          changemat_egfr_plus[as.character(egfr_grp0),as.character(alb_grp)] + 
          egfr_sglt_pop[i] * (tp[[i]]$year[j]-delayTreat[i]) * (tp[[i]]$year[j]<=effect_year_egfr_pop[i]) * (tp[[i]]$year[j]>delayTreat[i]) * !effect_sglt_RR
      }
      
      tp[[i]]$alb[j] <- tp[[i]]$alb[j-1] * 
        changemat_alb_mult[as.character(alb_grp)] + 
        UAR_sglt_pop[i] * tp[[i]]$year[j] * (tp[[i]]$year[j]<=effect_year_UAR_pop[i]) * (tp[[i]]$year[j]>delayTreat[i])
      
      # risks of HF MI stroke are determined by untreated eGFR if RRs are applied
      # egfr_grp0 used
      prob_hf1y <- fn_predict_prob(
        event = "HF", 
        eGFR = egfr0_value, 
        UACR = alb_value,
        Diabetes = 1
      ) * input_dia +
        fn_predict_prob(
          event = "HF", 
          eGFR = egfr0_value, 
          UACR = alb_value,
          Diabetes = 0
        ) * (1 - input_dia)
      
      prob_cvd1y <- fn_predict_prob(
        event = "CVD", 
        eGFR = egfr0_value, 
        UACR = alb_value,
        Diabetes = 1
      ) * input_dia +
        fn_predict_prob(
          event = "CVD", 
          eGFR = egfr0_value, 
          UACR = alb_value,
          Diabetes = 0
        ) * (1 - input_dia)
      
      tp[[i]]$prob_hf[j-1] <- prob_hf1y * (1-mt_hf_live) * 
        (RR_sta_hf_pop[i] * RR_BP_hf_pop[i] * RR_asp_hf_pop[i] * RR_sglt_hf_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      tp[[i]]$prob_fatalhf[j-1] <- prob_hf1y * mt_hf_live * 
        (RR_sta_hf_pop[i] * RR_BP_hf_pop[i] * RR_asp_hf_pop[i] * RR_sglt_hf_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      tp[[i]]$prob_mi[j-1] <- prob_cvd1y * prop_mi_live * (1-mt_mi_live) *
        (RR_sta_mi_pop[i] * RR_BP_mi_pop[i] * RR_asp_mi_pop[i] * RR_sglt_mi_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      tp[[i]]$prob_fatalmi[j-1] <- prob_cvd1y * prop_mi_live * mt_mi_live * 
        (RR_sta_mi_pop[i] * RR_BP_mi_pop[i] * RR_asp_mi_pop[i] * RR_sglt_mi_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      tp[[i]]$prob_str[j-1] <- prob_cvd1y * prop_str_live * (1-mt_str_live) * 
        (RR_sta_str_pop[i] * RR_BP_str_pop[i] * RR_asp_str_pop[i] * RR_sglt_str_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      tp[[i]]$prob_fatalstr[j-1] <- prob_cvd1y * prop_str_live * mt_str_live * 
        (RR_sta_str_pop[i] * RR_BP_str_pop[i] * RR_asp_str_pop[i] * RR_sglt_str_pop[i])^(tp[[i]]$year[j]>delayTreat[i])
      
      # risk of ESRD is still determined by treated eGFR
      # egfr_grp used
      prob_esrd1y <- fn_predict_prob(
        event = "ESRD", 
        eGFR = egfr_value, 
        UACR = alb_value,
        Diabetes = 1
      ) * input_dia +
        fn_predict_prob(
          event = "ESRD", 
          eGFR = egfr_value, 
          UACR = alb_value,
          Diabetes = 0
        ) * (1 - input_dia)
      
      tp[[i]]$prob_esrd_d[j-1] <- prob_esrd1y * prop_dial_live * (1-mt_esrd_live)
      tp[[i]]$prob_esrd_t[j-1] <- prob_esrd1y * prop_transp_live * (1-mt_esrd_live)
      tp[[i]]$prob_fatalesrd[j-1] <- prob_esrd1y * mt_esrd_live
      
      if(tp[[i]]$age[j-1]<=100){
        tp[[i]]$prob_otherdeath[j-1] <- merged_LT$other_cause[merged_LT$Age==tp[[i]]$age[j-1]]
      }else{
        tp[[i]]$prob_otherdeath[j-1] <- 1
      }
      tp[[i]]$prob_hf[j] <- NA
      tp[[i]]$prob_fatalhf[j] <- NA
      tp[[i]]$prob_mi[j] <- NA
      tp[[i]]$prob_fatalmi[j] <- NA
      tp[[i]]$prob_str[j] <- NA
      tp[[i]]$prob_fatalstr[j] <- NA
      tp[[i]]$prob_esrd_d[j] <- NA
      tp[[i]]$prob_esrd_t[j] <- NA
      tp[[i]]$prob_fatalesrd[j] <- NA
      tp[[i]]$prob_otherdeath[j] <- NA
    }
    
    # assumed 100% ESRD if (treated) eGFR < 15
    # adjust for sum of probabilities over 1
    for (j in 2:n_cycles) {
      if(tp[[i]]$eGFR[j-1]<15 & tp[[i]]$prob_otherdeath[j-1]<1){
        tp[[i]]$prob_hf[j-1] <- tp[[i]]$prob_fatalhf[j-1] <- tp[[i]]$prob_mi[j-1] <- tp[[i]]$prob_fatalmi[j-1] <- tp[[i]]$prob_str[j-1] <- tp[[i]]$prob_fatalstr[j-1] <- 0
        tp[[i]]$prob_esrd_d[j-1] <- prop_dial_live*(1-mt_esrd_live)
        tp[[i]]$prob_esrd_t[j-1] <- prop_transp_live*(1-mt_esrd_live)
        tp[[i]]$prob_fatalesrd[j-1] <- mt_esrd_live
        tp[[i]]$prob_otherdeath[j-1] <- 0
      }
      if(tp[[i]]$prob_otherdeath[j-1]==1){
        tp[[i]]$prob_hf[j-1] <- tp[[i]]$prob_fatalhf[j-1] <- tp[[i]]$prob_mi[j-1] <- tp[[i]]$prob_fatalmi[j-1] <- tp[[i]]$prob_str[j-1] <- tp[[i]]$prob_fatalstr[j-1] <- 0
        tp[[i]]$prob_esrd_d[j-1] <- tp[[i]]$prob_esrd_t[j-1] <- tp[[i]]$prob_fatalesrd[j-1] <- 0
      }
      
    }
    
    trace[[i]] <- data.frame(year = numeric(n_cycles),
                             age = numeric(n_cycles),
                             no_event = numeric(n_cycles),
                             hf = numeric(n_cycles),
                             fatalhf = numeric(n_cycles),
                             mi = numeric(n_cycles),
                             fatalmi = numeric(n_cycles),
                             str = numeric(n_cycles),
                             fatalstr = numeric(n_cycles),
                             esrd_d = numeric(n_cycles),
                             esrd_t = numeric(n_cycles),
                             fatalesrd = numeric(n_cycles),
                             otherdeath = numeric(n_cycles)
    )
    
    trace[[i]]$year <- tp[[i]]$year
    trace[[i]]$age <- tp[[i]]$age
    trace[[i]]$no_event[1] <- n_cohort
    
    
    for (j in 2:n_cycles) {
      
      for (k in 1:length(event_names)) {
        trace[[i]][j,event_names[k]] <- trace[[i]][j-1,event_names[k]] + trace[[i]]$no_event[j-1] * tp[[i]][j-1,prob_names[k]]
      }
      trace[[i]]$no_event[j] <- n_cohort - sum(trace[[i]][j,4:length(trace[[i]][j,])])
    }
    
    
  }
  
  # View(trace[[1]])
  # View(trace[[2]])
  # View(tp[[1]])
  # View(tp[[2]])
  
  ### sub-models ----
  sub_names <- c("hf", "mi", "str", "esrd_d", "esrd_t")
  n_sub <- length(sub_names)
  
  # cost set-up for sub-models
  c_mat <- array(dim = c(n_cat_egfr_6l, n_sub, n_age65, n_fy),
                 dimnames = list(eGFR = cat_egfr_6l_names,
                                 Event = sub_names,
                                 Age = l_age65,
                                 Year = l_first_year))
  
  for (cvd in levels(cost_ckdcvd$State)[-1]) {
    c_cvd_g3a_g5 <- cost_set_up_g3a_g5(data = cost_ckdcvd,
                                       cat1 = cat_egfr_6l_names[-(1:2)], 
                                       cat2 = l_age65, 
                                       cat3 = l_first_year,
                                       event = cvd)
    c_mat[-(1:2),cvd,,] <- c_cvd_g3a_g5[[1]] * (1 - input_dia) + c_cvd_g3a_g5[[2]] * input_dia
  }
  
  # G1 to G2, cvd
  c_mat[1:2,1:3,,"first"] <- c_g12_live # Assumption
  c_mat[1:2,1:3,,"subsequent"] <- c_mat[c(3,3),1:3,,"subsequent"] # Assumption
  
  # esrd_Dialysis
  c_mat[,"esrd_d",,"first"] <- c_esrd_d_4m_live + c_esrd_d_sub4m_live * 2
  c_mat[,"esrd_d",,"subsequent"] <- c_esrd_d_sub4m_live * 3
  
  # esrd_Transplant
  c_mat[,"esrd_t",,"first"] <- c_esrd_t_4m_live + c_esrd_t_sub4m_live * 2
  c_mat[,"esrd_t",,"subsequent"] <- c_esrd_t_sub4m_live * 3


  
  
  # c_mat <- array(data = rep(c(c(c_g12_hf,c_g12_hf,c_g3a_hf,c_g3b_hf,c_g4_hf,c_g5_hf),
  #                               rep(c(c_g12_cvd,c_g12_cvd,c_g3a_cvd,c_g3b_cvd,c_g4_cvd,c_g5_cvd),2),
  #                               rep(c_esrd_d_sub_live,n_cat_egfr_6l),
  #                               rep(c_esrd_t_sub_live,n_cat_egfr_6l)),
  #                             n_fy),
  #                dim = c(n_cat_egfr_6l, n_sub, n_fy),
  #                dimnames = list(eGFR = cat_egfr_6l_names,
  #                                Event = sub_names,
  #                                Year = l_first_year))
  
  
  # c_mat[,"hf",1] <- c_mat[,"hf",1] * (11/12) + c_hf_30d_live
  # c_mat[,"mi",1] <- c_mat[,"mi",1] * (11/12) + c_mi_30d_live
  # c_mat[,"str",1] <- c_mat[,"str",1] * (11/12) + c_str_30d_live
  
  
  # utility set-up
  u_mat <- array(data = c(rep(c(u_g12_live,u_g12_live,u_g3a_live,u_g3b_live,u_g4_live,u_g5_live), n_sub*n_fy)),
                 dim = c(n_cat_egfr_6l, n_sub, n_fy),
                 dimnames = list(eGFR = cat_egfr_6l_names,
                                 Event = sub_names,
                                 Year = l_first_year))
  
   
  u_mat[,"hf",1] <- u_mat[,"hf",1] - disu_hf_3m_live * 3/12 - disu_hf_sub_live * 9/12
  u_mat[,"mi",1] <- u_mat[,"mi",1] - disu_mi_3m_live * 3/12 - disu_mi_sub_live * 9/12
  u_mat[,"str",1] <- u_mat[,"str",1] - disu_str_3m_live * 3/12 - disu_str_sub_live * 9/12
  u_mat[,"hf",2] <- u_mat[,"hf",2] - disu_hf_sub_live
  u_mat[,"mi",2] <- u_mat[,"mi",2] - disu_mi_sub_live 
  u_mat[,"str",2] <- u_mat[,"str",2] - disu_str_sub_live
  
  u_mat[,"esrd_d",1] <- u_esrd_d_4m_live * 4/12 + u_esrd_d_sub_live * 8/12
  u_mat[,"esrd_t",1] <- u_esrd_t_4m_live * 4/12 + u_esrd_t_sub_live * 8/12
  u_mat[,"esrd_d",2] <- u_esrd_d_sub_live
  u_mat[,"esrd_t",2] <- u_esrd_t_sub_live
  
  # all-cause mortality risk function, NHLBI Pooled Cohorts Study
  
  # Re-calibration
  # "These risk functions were recalibrated for each sub-model, to replicate observed all-cause mortality for patients living with HF, MI, stroke, dialysis, and kidney transplantation in a range of published studies."
  cali_hf <- -1.2744402947477
  cali_mi <- -3.05705995765881
  cali_str <- -2.63721449871134
  cali_esrd_d <- -2.04315702469463
  cali_esrd_t <- -2.16337286506919
  
  alpha <- c(cali_hf, cali_mi, cali_str, cali_esrd_d, cali_esrd_t)
  
  Mortality_Risk_Function <- read_csv("data/Mortality_Risk_Function.csv", 
                                      col_types = cols(Variable = col_character(), 
                                                       Mean = col_double(), 
                                                       Coeff = col_double(), 
                                                       SE = col_double()))
  
  age_breaks <- c(18,seq(25, 95, by = 5),Inf)
  age_mid <- NA
  for (i in 1:(length(age_breaks)-2)){
    age_mid[i] <- age_breaks[i] + (age_breaks[i+1]-1 - age_breaks[i])/2
  }
  age_mid[length(age_breaks)-1] <- age_breaks[length(age_breaks)-1] + 2
  
  egfr_breaks <- c(seq(0, 60, by = 15),90,100)
  egfr_mid <- NA
  for (i in 1:(length(egfr_breaks)-1)){
    egfr_mid[i] <- mean(egfr_breaks[i:(i+1)])
  }
  
  payoffs <- function(age_grp, egfr_grp, first_event,forward.year=30){
    age_use <- NA
    egfr_use <- NA
    alive <- NA
    lp <- NA
    xb <- NA
    pr <- NA
    undisc_c <- NA
    disc_c <- NA
    accum_d_c <- NA
    undisc_u <- NA
    disc_u <- NA
    accum_d_u <- NA
    
    age_use[1] <- age_mid[age_grp]
    egfr_use[1] <- egfr_mid[egfr_grp]
    
    # survival
    alive[1] <- 1
    for (j in 1:forward.year){
      Mortality_Risk_Function[1,"Mean"] <- age_use[j]-55
      Mortality_Risk_Function[9,"Mean"] <- egfr_use[j]
      Mortality_Risk_Function[10,"Mean"] <- (age_use[j]-55) * Mortality_Risk_Function[2,"Mean"]
      Mortality_Risk_Function[11,"Mean"] <- (age_use[j]-55) * Mortality_Risk_Function[4,"Mean"]
      Mortality_Risk_Function[12,"Mean"] <- (age_use[j]-55) * Mortality_Risk_Function[8,"Mean"]
      lp[j] <- sum(Mortality_Risk_Function$Mean * Mortality_Risk_Function$Coeff)
      xb[j] <- lp[j] + alpha[match(first_event,sub_names)]
      pr[j] <- 1/(1+exp(-xb[j]))
      if (age_use[j]>=100) pr[j] <- 1
      alive[j+1] <- alive[j] * (1-pr[j])
      age_use[j+1] <- age_use[j] + 1
      egfr_grp4l <- cut(egfr_use[j], breaks = c(0,30,45,60,Inf),right=FALSE,labels=c("G4/5", "G3b", "G3a", "G1/2"))
      egfr_grp6l <- cut(egfr_use[j], breaks = egfr_breaks,right=FALSE,labels=c("G5", "G4","G3b", "G3a", "G2", "G1"))
      
      if(egfr_use[j]<15) {
        egfr_use[j+1] <- egfr_use[j]
      }else{
        egfr_use[j+1] <- egfr_use[j] + changemat_egfr_plus[as.character(egfr_grp4l),"A2"]
      }
      lp[j+1] <- NA
      xb[j+1] <- NA
      pr[j+1] <- NA
      
      # cost
      undisc_c[j] <- alive[j] * c_mat[as.character(egfr_grp6l),first_event,(age_use[j]>=65)+1,(j>1)+1]
      undisc_c[j+1] <-0
      
      disc_c[j] <- undisc_c[j] *(1/((1+disc.cost)^(j-1)))
      disc_c[j+1] <-0
      
      accum_d_c[j] <- sum(disc_c[1:j]) 
      
      # utility
      undisc_u[j] <- alive[j] * u_mat[as.character(egfr_grp6l),first_event,(j>1)+1]
      undisc_u[j+1] <-0
      
      disc_u[j] <- undisc_u[j] *(1/((1+disc.outcomes)^(j-1)))
      disc_u[j+1] <-0
      
      accum_d_u[j] <- sum(disc_u[1:j]) 
    }
    sub_survival <- cbind(age_use,egfr_use,lp,xb,pr,alive,undisc_c,disc_c,undisc_u,disc_u)
    
    # return accumulated cost and QALY by year

    return(matrix(data = c(accum_d_c[1:forward.year],accum_d_u[1:forward.year]),
                  ncol = 2,
                  dimnames = list(year = 1:forward.year,
                                  outcomes = c("Cost", "QALY"))))
  }
  
  agegrp_names <- c("[18,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", 
                    "[45,50)", "[50,55)", "[55,60)", "[60,65)", "[65,70)", 
                    "[70,75)", "[75,80)", "[80,85)", "[85,90)", "[90,95)", 
                    "[95,Inf)")
  egfrgrp_names <- c("[0,15)", "[15,30)", "[30,45)", "[45,60)", "[60,90)", "[90,100)")
  
  n_agegrp <- length(agegrp_names)
  n_egfrgrp <- length(egfrgrp_names)

  
  
  if (input_update_sub == F) {
    load("data/sub_model_po.RData")
  }
  
  if (input_update_sub == T) {
    sub_model_po <- array(
      dim = c(n_agegrp, n_egfrgrp, n_sub, 30, 2),
      dimnames = list(Age = agegrp_names,
                      eGFR = egfrgrp_names,
                      first_event = sub_names,
                      year = 1:30,
                      outcomes = c("Cost", "QALY")))
    
    for (i in 1:n_agegrp) { 
      #loop for age groups
      for (j in 1:n_egfrgrp) {
        #loop for eGFR groups
        for(name in sub_names) {
          #loop for first event: "hf", "mi", "str", "esrd_d", "esrd_t"
          sub_model_po[i,j,name,,] <- payoffs(age_grp = i, egfr_grp = j, first_event = name)
        }
      }
    }
  }

  
  remainyear <- c(pmin(pmax(n_cycles_short - 1:n_cycles_short, 0), 30), rep(0, n_cycles - n_cycles_short))

  payoffs_search <- function(age, eGFR, first_event, timehorz = 30) {
    if(timehorz==0) {
      return(c(0,0))
    }else{
      age_grp <- cut(age, breaks = age_breaks, right=FALSE)
      egfr_grp <- cut(eGFR, breaks = c(egfr_breaks[-length(egfr_breaks)], Inf), right=FALSE)
      return(sub_model_po[age_grp,egfr_grp,first_event,timehorz,])
    }
  }
  # payoffs_search(age=57,eGFR=52.5,"hf")
  
  cost_int <- matrix(NA,2,n_cycles)
  disc_cost_int <- matrix(NA,2,n_cycles)
  qaly_int <- matrix(NA,2,n_cycles)
  disc_qaly_int <- matrix(NA,2,n_cycles)
  
  cost_noevent <- matrix(NA,2,n_cycles)
  disc_cost_noevent <- matrix(NA,2,n_cycles)
  qaly_noevent <- matrix(NA,2,n_cycles)
  disc_qaly_noevent <- matrix(NA,2,n_cycles)
  
  cost_acute_hf <- matrix(NA,2,n_cycles)
  cost_acute_mi <- matrix(NA,2,n_cycles)
  cost_acute_str <- matrix(NA,2,n_cycles)
  cost_acute_esrd_d <- matrix(NA,2,n_cycles)
  cost_acute_esrd_t <- matrix(NA,2,n_cycles)
  
  disc_cost_acute_hf <- matrix(NA,2,n_cycles)
  disc_cost_acute_mi <- matrix(NA,2,n_cycles)
  disc_cost_acute_str <- matrix(NA,2,n_cycles)
  disc_cost_acute_esrd <- matrix(NA,2,n_cycles)
  
  cost_acute <- matrix(NA,2,n_cycles)
  disc_cost_acute <- matrix(NA,2,n_cycles)
  
  cq_sub_hf<-cq_sub_mi<-cq_sub_str<-cq_sub_esrd_d<-cq_sub_esrd_t<-
    array(NA, dim = c(2, n_cycles, 2), dimnames = list(arm = c("comparator","intervention"),
                                                       cycle = seq(to=n_cycles),
                                                       payoff = c("cost","QALY")))
  
  cost_sub_hf <- matrix(NA,2,n_cycles)
  cost_sub_mi <- matrix(NA,2,n_cycles)
  cost_sub_str <- matrix(NA,2,n_cycles)
  cost_sub_esrd_d <- matrix(NA,2,n_cycles)
  cost_sub_esrd_t <- matrix(NA,2,n_cycles)
  
  disc_cost_sub_hf <- matrix(NA,2,n_cycles)
  disc_cost_sub_mi <- matrix(NA,2,n_cycles)
  disc_cost_sub_str <- matrix(NA,2,n_cycles)
  disc_cost_sub_esrd <- matrix(NA,2,n_cycles)
  
  cost_sub <- matrix(NA,2,n_cycles)
  disc_cost_sub <- matrix(NA,2,n_cycles)
  
  qaly_sub_hf <- matrix(NA,2,n_cycles)
  qaly_sub_mi <- matrix(NA,2,n_cycles)
  qaly_sub_str <- matrix(NA,2,n_cycles)
  qaly_sub_esrd_d <- matrix(NA,2,n_cycles)
  qaly_sub_esrd_t <- matrix(NA,2,n_cycles)
  
  disc_qaly_sub_hf <- matrix(NA,2,n_cycles)
  disc_qaly_sub_mi <- matrix(NA,2,n_cycles)
  disc_qaly_sub_str <- matrix(NA,2,n_cycles)
  disc_qaly_sub_esrd <- matrix(NA,2,n_cycles)
  
  qaly_sub <- matrix(NA,2,n_cycles)
  disc_qaly_sub <- matrix(NA,2,n_cycles)
  
  total_cost <- matrix(NA,2,n_cycles)
  total_qaly <- matrix(NA,2,n_cycles)
  
  for (i in 1:2) {
    for (j in 2:n_cycles) {
      egfr_grp6l <- cut(tp[[i]]$eGFR[j], 
                        breaks = c(egfr_breaks[-length(egfr_breaks)], Inf),
                        right=FALSE,labels=c("G5", "G4","G3b", "G3a", "G2", "G1"))
      alb_grp <- cut(tp[[i]]$alb[j], breaks = c(0,30,300,Inf),right=FALSE,labels=c("A1", "A2", "A3"))
      egfr_value <- tp[[i]]$eGFR[j]
      alb_value <- tp[[i]]$alb[j]
      # intervention
      cost_int[i,j] <- (trace[[i]]$year[j]>delayTreat[i]) * (trace[[i]]$no_event[j]+trace[[i]]$hf[j]+trace[[i]]$mi[j]+trace[[i]]$str[j])*
        (c_sta_pop[i]+c_sta_diabetes_pop[i]+c_sta_mae_pop[i]+c_sta_sae_pop[i]+
           c_BP_pop[i]+c_BP_itae_pop[i]+c_BP_sae_pop[i]+
           c_asp_pop[i]+c_asp_mb_pop[i]+
           c_sglt_pop[i])
      
      disc_cost_int[i,j] <- cost_int[i,j] * (1/((1+disc.cost)^(j-1)))
      
      qaly_int[i,j] <- (trace[[i]]$year[j]>delayTreat[i]) * (trace[[i]]$no_event[j]+trace[[i]]$hf[j]+trace[[i]]$mi[j]+trace[[i]]$str[j])*
        (-disu_sta_diabetes_pop[i]-disu_sta_mae_pop[i]-disu_sta_sae_pop[i]-
           disu_BP_itae_pop[i]-disu_BP_sae_pop[i]-
           disu_asp_mb_pop[i])
      
      disc_qaly_int[i,j] <- qaly_int[i,j] * (1/((1+disc.outcomes)^(j-1)))
      
      # no event
      cost_noevent[i,j] <- trace[[i]]$no_event[j] * 
        (fn_predict_cost("CKD", 
                         eGFR = egfr_value,
                         UACR = alb_value,
                         Age_LT65 = trace[[i]]$age[j]<65,
                         Baseline_T2D = 0,
                         Year1 = j==2) * 
           (1 - input_dia) +
           fn_predict_cost("CKD", 
                           eGFR = egfr_value,
                           UACR = alb_value,
                           Age_LT65 = trace[[i]]$age[j]<65,
                           Baseline_T2D = 1,
                           Year1 = j==2) * 
           input_dia)
      
      # c_mat_noevent[as.character(egfr_grp6l),
      #               as.character(alb_grp),
      #               ifelse(trace[[i]]$age[j]<65, "lt 65", "ge 65"),
      #               ifelse(j==2, "first", "subsequent")]
     
      disc_cost_noevent[i,j] <- cost_noevent[i,j] * (1/((1+disc.cost)^(j-1)))
      
      qaly_noevent[i,j] <- trace[[i]]$no_event[j] * u_mat_noevent[as.character(egfr_grp6l),as.character(alb_grp)]
      disc_qaly_noevent[i,j] <- qaly_noevent[i,j] * (1/((1+disc.outcomes)^(j-1)))
      
      # acute cost of events not yet covered (for fatal events)
      cost_acute_hf[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_fatalhf[j-1] * c_hf_30d_live 
      cost_acute_mi[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_fatalmi[j-1] * c_mi_30d_live 
      cost_acute_str[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_fatalstr[j-1] * c_str_30d_live  
      cost_acute_esrd_d[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_fatalesrd[j-1] * prop_dial_live * c_esrd_d_4m_live
      cost_acute_esrd_t[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_fatalesrd[j-1] * prop_transp_live * c_esrd_t_4m_live
      
      cost_acute[i,j] <- cost_acute_hf[i,j] + cost_acute_mi[i,j] + cost_acute_str[i,j] + cost_acute_esrd_d[i,j] + cost_acute_esrd_t[i,j]
      disc_cost_acute_hf[i,j] <- cost_acute_hf[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_acute_mi[i,j] <- cost_acute_mi[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_acute_str[i,j] <- cost_acute_str[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_acute_esrd[i,j] <- cost_acute_esrd_d[i,j] *(1/((1+disc.cost)^(j-1))) + cost_acute_esrd_t[i,j] *(1/((1+disc.cost)^(j-1)))
      
      disc_cost_acute[i,j] <- cost_acute[i,j] *(1/((1+disc.cost)^(j-1)))
      
      # cost from sub-models
      cq_sub_hf[i,j,] <- payoffs_search(age=tp[[i]]$age[j],eGFR=tp[[i]]$eGFR[j],"hf",remainyear[j])
      cq_sub_mi[i,j,] <- payoffs_search(age=tp[[i]]$age[j],eGFR=tp[[i]]$eGFR[j],"mi",remainyear[j])
      cq_sub_str[i,j,] <- payoffs_search(age=tp[[i]]$age[j],eGFR=tp[[i]]$eGFR[j],"str",remainyear[j])
      cq_sub_esrd_d[i,j,] <- payoffs_search(age=tp[[i]]$age[j],eGFR=tp[[i]]$eGFR[j],"esrd_d",remainyear[j])
      cq_sub_esrd_t[i,j,] <- payoffs_search(age=tp[[i]]$age[j],eGFR=tp[[i]]$eGFR[j],"esrd_t",remainyear[j])
      
      cost_sub_hf[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_hf[j-1] * 
        cq_sub_hf[i,j,1]
      cost_sub_mi[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_mi[j-1] * 
        cq_sub_mi[i,j,1]
      cost_sub_str[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_str[j-1] * 
        cq_sub_str[i,j,1]
      cost_sub_esrd_d[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_esrd_d[j-1] * 
        cq_sub_esrd_d[i,j,1]
      cost_sub_esrd_t[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_esrd_t[j-1] * 
        cq_sub_esrd_t[i,j,1]
      
      cost_sub[i,j] <- cost_sub_hf[i,j] + cost_sub_mi[i,j] + cost_sub_str[i,j] + cost_sub_esrd_d[i,j] + cost_sub_esrd_t[i,j]
      disc_cost_sub_hf[i,j] <- cost_sub_hf[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_sub_mi[i,j] <- cost_sub_mi[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_sub_str[i,j] <- cost_sub_str[i,j] *(1/((1+disc.cost)^(j-1)))
      disc_cost_sub_esrd[i,j] <- cost_sub_esrd_d[i,j] *(1/((1+disc.cost)^(j-1))) + cost_sub_esrd_t[i,j] *(1/((1+disc.cost)^(j-1)))
      
      disc_cost_sub[i,j] <- cost_sub[i,j] *(1/((1+disc.cost)^(j-1)))
      
      # QALY from sub-models
      qaly_sub_hf[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_hf[j-1] * 
        cq_sub_hf[i,j,2]
      qaly_sub_mi[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_mi[j-1] * 
        cq_sub_mi[i,j,2]
      qaly_sub_str[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_str[j-1] * 
        cq_sub_str[i,j,2]
      qaly_sub_esrd_d[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_esrd_d[j-1] * 
        cq_sub_esrd_d[i,j,2]
      qaly_sub_esrd_t[i,j] <- trace[[i]]$no_event[j-1] * tp[[i]]$prob_esrd_t[j-1] * 
        cq_sub_esrd_t[i,j,2]
      
      qaly_sub[i,j] <- qaly_sub_hf[i,j] + qaly_sub_mi[i,j] + qaly_sub_str[i,j] + qaly_sub_esrd_d[i,j] + qaly_sub_esrd_t[i,j]
      disc_qaly_sub_hf[i,j] <- qaly_sub_hf[i,j] *(1/((1+disc.outcomes)^(j-1)))
      disc_qaly_sub_mi[i,j] <- qaly_sub_mi[i,j] *(1/((1+disc.outcomes)^(j-1)))
      disc_qaly_sub_str[i,j] <- qaly_sub_str[i,j] *(1/((1+disc.outcomes)^(j-1)))
      disc_qaly_sub_esrd[i,j] <- qaly_sub_esrd_d[i,j] *(1/((1+disc.outcomes)^(j-1))) + qaly_sub_esrd_t[i,j] *(1/((1+disc.outcomes)^(j-1)))
      
      disc_qaly_sub[i,j] <- qaly_sub[i,j] *(1/((1+disc.outcomes)^(j-1)))
      
      total_cost[i,j] <- disc_cost_int[i,j] + disc_cost_noevent[i,j] + disc_cost_acute[i,j] + disc_cost_sub[i,j]
      total_qaly[i,j] <- disc_qaly_int[i,j] + disc_qaly_noevent[i,j] + disc_qaly_sub[i,j]
    }
  }
  
  # short-term cost components
  ## no event
  disc_cost_s_noevent_comp <- sum(disc_cost_noevent[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_noevent_int <- sum(disc_cost_noevent[2,1:n_cycles_short],na.rm=T)
  
  ## ESRD
  disc_cost_s_esrd_comp <- sum(disc_cost_sub_esrd[1,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_esrd[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_esrd_int <- sum(disc_cost_sub_esrd[2,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_esrd[2,1:n_cycles_short],na.rm=T)
  
  ## Stroke
  disc_cost_s_str_comp <- sum(disc_cost_sub_str[1,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_str[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_str_int <- sum(disc_cost_sub_str[2,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_str[2,1:n_cycles_short],na.rm=T)
  
  ## MI
  disc_cost_s_mi_comp <- sum(disc_cost_sub_mi[1,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_mi[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_mi_int <- sum(disc_cost_sub_mi[2,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_mi[2,1:n_cycles_short],na.rm=T)
  
  ## HF
  disc_cost_s_hf_comp <- sum(disc_cost_sub_hf[1,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_hf[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_hf_int <- sum(disc_cost_sub_hf[2,1:n_cycles_short],na.rm=T) + sum(disc_cost_acute_hf[2,1:n_cycles_short],na.rm=T)
  
  ## Medication
  disc_cost_s_int_comp <- sum(disc_cost_int[1,1:n_cycles_short],na.rm=T)
  disc_cost_s_int_int <- sum(disc_cost_int[2,1:n_cycles_short],na.rm=T)
  
  
  # short-term QALY components
  ## no event
  disc_qaly_s_noevent_comp <- sum(disc_qaly_noevent[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_noevent_int <- sum(disc_qaly_noevent[2,1:n_cycles_short],na.rm=T)
  
  ## ESRD
  disc_qaly_s_esrd_comp <- sum(disc_qaly_sub_esrd[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_esrd_int <- sum(disc_qaly_sub_esrd[2,1:n_cycles_short],na.rm=T)
  
  ## Stroke
  disc_qaly_s_str_comp <- sum(disc_qaly_sub_str[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_str_int <- sum(disc_qaly_sub_str[2,1:n_cycles_short],na.rm=T)
  
  ## MI
  disc_qaly_s_mi_comp <- sum(disc_qaly_sub_mi[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_mi_int <- sum(disc_qaly_sub_mi[2,1:n_cycles_short],na.rm=T)
  
  ## HF
  disc_qaly_s_hf_comp <- sum(disc_qaly_sub_hf[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_hf_int <- sum(disc_qaly_sub_hf[2,1:n_cycles_short],na.rm=T)
  
  # total cost and QALY
  disc_cost_s_comp <- sum(total_cost[1,1:n_cycles_short],na.rm=T)
  disc_qaly_s_comp <- sum(total_qaly[1,1:n_cycles_short],na.rm=T)
  
  disc_cost_s_int <- sum(total_cost[2,1:n_cycles_short],na.rm=T)
  disc_qaly_s_int <- sum(total_qaly[2,1:n_cycles_short],na.rm=T)
  
  disc_incr.cost_s <- disc_cost_s_int - disc_cost_s_comp
  disc_incr.qaly_s <- disc_qaly_s_int - disc_qaly_s_comp
  
  ICER_s <- disc_incr.cost_s / disc_incr.qaly_s
  
  # life time cost components
  ## no event
  disc_cost_noevent_comp <- sum(disc_cost_noevent[1,],na.rm=T)
  disc_cost_noevent_int <- sum(disc_cost_noevent[2,],na.rm=T)
  
  ## ESRD
  disc_cost_esrd_comp <- sum(disc_cost_sub_esrd[1,],na.rm=T) + sum(disc_cost_acute_esrd[1,],na.rm=T)
  disc_cost_esrd_int <- sum(disc_cost_sub_esrd[2,],na.rm=T) + sum(disc_cost_acute_esrd[2,],na.rm=T)
  
  ## Stroke
  disc_cost_str_comp <- sum(disc_cost_sub_str[1,],na.rm=T) + sum(disc_cost_acute_str[1,],na.rm=T)
  disc_cost_str_int <- sum(disc_cost_sub_str[2,],na.rm=T) + sum(disc_cost_acute_str[2,],na.rm=T)
  
  ## MI
  disc_cost_mi_comp <- sum(disc_cost_sub_mi[1,],na.rm=T) + sum(disc_cost_acute_mi[1,],na.rm=T)
  disc_cost_mi_int <- sum(disc_cost_sub_mi[2,],na.rm=T) + sum(disc_cost_acute_mi[2,],na.rm=T)
  
  ## HF
  disc_cost_hf_comp <- sum(disc_cost_sub_hf[1,],na.rm=T) + sum(disc_cost_acute_hf[1,],na.rm=T)
  disc_cost_hf_int <- sum(disc_cost_sub_hf[2,],na.rm=T) + sum(disc_cost_acute_hf[2,],na.rm=T)
  
  ## Medication
  disc_cost_int_comp <- sum(disc_cost_int[1,],na.rm=T)
  disc_cost_int_int <- sum(disc_cost_int[2,],na.rm=T)
  
  
  # life time QALY components
  ## no event
  disc_qaly_noevent_comp <- sum(disc_qaly_noevent[1,],na.rm=T)
  disc_qaly_noevent_int <- sum(disc_qaly_noevent[2,],na.rm=T)
  
  ## ESRD
  disc_qaly_esrd_comp <- sum(disc_qaly_sub_esrd[1,],na.rm=T)
  disc_qaly_esrd_int <- sum(disc_qaly_sub_esrd[2,],na.rm=T)
  
  ## Stroke
  disc_qaly_str_comp <- sum(disc_qaly_sub_str[1,],na.rm=T)
  disc_qaly_str_int <- sum(disc_qaly_sub_str[2,],na.rm=T)
  
  ## MI
  disc_qaly_mi_comp <- sum(disc_qaly_sub_mi[1,],na.rm=T)
  disc_qaly_mi_int <- sum(disc_qaly_sub_mi[2,],na.rm=T)
  
  ## HF
  disc_qaly_hf_comp <- sum(disc_qaly_sub_hf[1,],na.rm=T)
  disc_qaly_hf_int <- sum(disc_qaly_sub_hf[2,],na.rm=T)
  
  
  # total cost and QALY
  disc_cost_comp <- sum(total_cost[1,],na.rm=T)
  disc_qaly_comp <- sum(total_qaly[1,],na.rm=T)
  
  disc_cost_int <- sum(total_cost[2,],na.rm=T)
  disc_qaly_int <- sum(total_qaly[2,],na.rm=T)
  
  disc_incr.cost <- disc_cost_int - disc_cost_comp
  disc_incr.qaly <- disc_qaly_int - disc_qaly_comp
  
  ICER <- disc_incr.cost / disc_incr.qaly
  
  # Save the sub model payoffs for fast running
  save(sub_model_po, file = "data/sub_model_po.RData")
  
  # Create the HTML table
  CEAtable_s <- htmlTable(matrix(c(paste0("$",formatC(disc_cost_s_comp, format = "f", big.mark = ",", digits = 0)),
                                 formatC(disc_qaly_s_comp, format = "f", digits = 2), 
                                 "-",
                                 paste0("$",formatC(disc_cost_s_int, format = "f", big.mark = ",", digits = 0)), 
                                 formatC(disc_qaly_s_int, format = "f", digits = 2), 
                                 "-",
                                 paste0("$",formatC(disc_incr.cost_s, format = "f", big.mark = ",", digits = 0)), 
                                 formatC(disc_incr.qaly_s, format = "f", digits = 2), 
                                 ifelse(ICER_s<0, 
                                        ifelse(disc_incr.qaly_s<0, 
                                               "Dominated", 
                                               "Dominant"), 
                                        paste0("$",formatC(ICER_s, format = "f", big.mark = ",", digits = 0),"/QALY"))), 
                               nrow = 3, byrow = TRUE),
                        header = c("Costs", "QALYs", "ICER"),
                        rnames = c(comparator, intervention, "Incremental"),
                        css.cell = "padding: 10px;",
                        css.table = "width: 700px;")
  
  CEAtable <- htmlTable(matrix(c(paste0("$",formatC(disc_cost_comp, format = "f", big.mark = ",", digits = 0)),
                                 formatC(disc_qaly_comp, format = "f", digits = 2), 
                                 "-",
                                 paste0("$",formatC(disc_cost_int, format = "f", big.mark = ",", digits = 0)), 
                                 formatC(disc_qaly_int, format = "f", digits = 2), 
                                 "-",
                                 paste0("$",formatC(disc_incr.cost, format = "f", big.mark = ",", digits = 0)), 
                                 formatC(disc_incr.qaly, format = "f", digits = 2), 
                                 ifelse(ICER<0, 
                                        ifelse(disc_incr.qaly<0, 
                                               "Dominated", 
                                               "Dominant"),
                                        paste0("$",formatC(ICER, format = "f", big.mark = ",", digits = 0),"/QALY"))), 
                               nrow = 3, byrow = TRUE),
                        header = c("Costs", "QALYs", "ICER"),
                        rnames = c(comparator, intervention, "Incremental"),
                        css.cell = "padding: 10px;",
                        css.table = "width: 700px;")
  
  # CEAtable
  for (i in 1:2) {
    trace[[i]] <- transform(trace[[i]], 
                            esrd_all = trace[[i]]$esrd_d + trace[[i]]$esrd_t,
                            death_all = trace[[i]]$otherdeath + trace[[i]]$fatalhf + trace[[i]]$fatalmi + trace[[i]]$fatalstr + trace[[i]]$fatalesrd)
  }
  
  plot1 <- ggplot(trace[[1]][trace[[1]]$year <= 40, ], aes(x = year)) +
    geom_line(aes(y = no_event, color = "No event"), linewidth = 1.3) +
    geom_line(aes(y = hf, color = "HF"), linewidth = 1.3) +
    geom_line(aes(y = mi, color = "MI"), linewidth = 1.3) +
    geom_line(aes(y = str, color = "Stroke"), linewidth = 1.3) +
    geom_line(aes(y = esrd_all, color = "ESRD"), linewidth = 1.3) +
    geom_line(aes(y = death_all, color = "Death"), linewidth = 1.3) +
    scale_color_manual("", values = c("No event" = "lightblue", "HF" = "darkred", "MI" = "gold", "Stroke" = "orange", "ESRD" = "purple", "Death" = "black")
                       , breaks = c("No event", "HF", "MI", "Stroke", "ESRD", "Death")
                       ) +
    labs(title = comparator,
         x = "Year",
         y = "Cumulative Probability of Event") +
    theme_minimal()
  
  plot2 <- ggplot(trace[[2]][trace[[2]]$year <= 40, ], aes(x = year)) +
    geom_line(aes(y = no_event, color = "No event"), linewidth = 1.3) +
    geom_line(aes(y = hf, color = "HF"), linewidth = 1.3) +
    geom_line(aes(y = mi, color = "MI"), linewidth = 1.3) +
    geom_line(aes(y = str, color = "Stroke"), linewidth = 1.3) +
    geom_line(aes(y = esrd_all, color = "ESRD"), linewidth = 1.3) +
    geom_line(aes(y = death_all, color = "Death"), linewidth = 1.3) +
    scale_color_manual("", values = c("No event" = "lightblue", "HF" = "darkred", "MI" = "gold", "Stroke" = "orange", "ESRD" = "purple", "Death" = "black")
                       , breaks = c("No event", "HF", "MI", "Stroke", "ESRD", "Death")
    ) +
    labs(title = intervention,
         x = "Year",
         y = "Cumulative Probability of Event") +
    theme_minimal()
  
  
  # For short term events
  for (i in 1:2) {
    trace[[i]]$cr_death <- trace[[i]]$fatalesrd + trace[[i]]$fatalhf + trace[[i]]$fatalmi + trace[[i]]$fatalstr
  }
  
  event_short_com <- unlist(trace[[1]][trace[[1]]$year == n_cycles_short - 1, c("no_event", "esrd_all", "str", "mi", "hf", "cr_death", "otherdeath")])
  event_short_int <- unlist(trace[[2]][trace[[2]]$year == n_cycles_short - 1, c("no_event", "esrd_all", "str", "mi", "hf", "cr_death", "otherdeath")])
  
  event_lt_com <- unlist(trace[[1]][n_cycles, c("esrd_all", "str", "mi", "hf", "cr_death", "otherdeath")])
  event_lt_int <- unlist(trace[[2]][n_cycles, c("esrd_all", "str", "mi", "hf", "cr_death", "otherdeath")])
  
  # sum of the time to event: sumproduct(trace_n - trace_n-1, n)
  sum_tte_esrd <- t(diff(cbind(trace[[1]][,"esrd_all"], trace[[2]][,"esrd_all"]))) %*% 1:(n_cycles - 1)
  sum_tte_str <- t(diff(cbind(trace[[1]][,"str"], trace[[2]][,"str"]))) %*% 1:(n_cycles - 1)
  sum_tte_mi <- t(diff(cbind(trace[[1]][,"mi"], trace[[2]][,"mi"]))) %*% 1:(n_cycles - 1)
  sum_tte_hf <- t(diff(cbind(trace[[1]][,"hf"], trace[[2]][,"hf"]))) %*% 1:(n_cycles - 1)
  sum_tte_cr_death <- t(diff(cbind(trace[[1]][,"cr_death"], trace[[2]][,"cr_death"]))) %*% 1:(n_cycles - 1)
  sum_tte_otherdeath <- t(diff(cbind(trace[[1]][,"otherdeath"], trace[[2]][,"otherdeath"]))) %*% 1:(n_cycles - 1)
  
  tte_esrd_com <- unname(sum_tte_esrd[1]/event_lt_com["esrd_all"])
  tte_str_com <- unname(sum_tte_str[1]/event_lt_com["str"])
  tte_mi_com <- unname(sum_tte_mi[1]/event_lt_com["mi"])
  tte_hf_com <- unname(sum_tte_hf[1]/event_lt_com["hf"])
  tte_cr_death_com <- unname(sum_tte_cr_death[1]/event_lt_com["cr_death"])
  tte_otherdeath_com <- unname(sum_tte_otherdeath[1]/event_lt_com["otherdeath"])
  
  tte_esrd_int <- unname(sum_tte_esrd[2]/event_lt_int["esrd_all"])
  tte_str_int <- unname(sum_tte_str[2]/event_lt_int["str"])
  tte_mi_int <- unname(sum_tte_mi[2]/event_lt_int["mi"])
  tte_hf_int <- unname(sum_tte_hf[2]/event_lt_int["hf"])
  tte_cr_death_int <- unname(sum_tte_cr_death[2]/event_lt_int["cr_death"])
  tte_otherdeath_int <- unname(sum_tte_otherdeath[2]/event_lt_int["otherdeath"])
  
  tte_com <- data.frame(Event = c("Non-fatal ESRD", "Non-fatal stroke", "Non-fatal MI", "Non-fatal HF", "CR death", "Other death"), 
                        TimeToEvent = c(tte_esrd_com, tte_str_com, tte_mi_com, tte_hf_com, tte_cr_death_com, tte_otherdeath_com))
  tte_int <- data.frame(Event = c("Non-fatal ESRD", "Non-fatal stroke", "Non-fatal MI", "Non-fatal HF", "CR death", "Other death"), 
                        TimeToEvent = c(tte_esrd_int, tte_str_int, tte_mi_int, tte_hf_int, tte_cr_death_int, tte_otherdeath_int))
  
  components <- c("No event", "ESRD", "Stroke", "MI", "HF", "Medication")
  components_q <- c("No event", "ESRD", "Stroke", "MI", "HF")
  # For short term disagregated cost 
  cost_short <- data.frame(
    category = c(rep("SOC", 6),rep("Intervention", 6)),
    cost = c(disc_cost_s_noevent_comp,
             disc_cost_s_esrd_comp,
             disc_cost_s_str_comp,
             disc_cost_s_mi_comp,
             disc_cost_s_hf_comp,
             disc_cost_s_int_comp,
             
             disc_cost_s_noevent_int,
             disc_cost_s_esrd_int,
             disc_cost_s_str_int,
             disc_cost_s_mi_int,
             disc_cost_s_hf_int,
             disc_cost_s_int_int),
    component = rep(components, 2)
  )
  
  # For short term disagregated QALY
  qaly_short <- data.frame(
    category = c(rep("SOC", 5),rep("Intervention", 5)),
    qaly = c(disc_qaly_s_noevent_comp,
             disc_qaly_s_esrd_comp,
             disc_qaly_s_str_comp,
             disc_qaly_s_mi_comp,
             disc_qaly_s_hf_comp,
             
             disc_qaly_s_noevent_int,
             disc_qaly_s_esrd_int,
             disc_qaly_s_str_int,
             disc_qaly_s_mi_int,
             disc_qaly_s_hf_int),
    component = rep(components_q, 2)
  )
  
  
  # For disagregated cost lifetime
  cost_lt <- data.frame(
    category = c(rep("SOC", 6),rep("Intervention", 6)),
    cost = c(disc_cost_noevent_comp,
             disc_cost_esrd_comp,
             disc_cost_str_comp,
             disc_cost_mi_comp,
             disc_cost_hf_comp,
             disc_cost_int_comp,
             
             disc_cost_noevent_int,
             disc_cost_esrd_int,
             disc_cost_str_int,
             disc_cost_mi_int,
             disc_cost_hf_int,
             disc_cost_int_int),
    component = rep(components, 2)
  )
  
  # For disagregated QALY lifetime
  qaly_lt <- data.frame(
    category = c(rep("SOC", 5),rep("Intervention", 5)),
    qaly = c(disc_qaly_noevent_comp,
             disc_qaly_esrd_comp,
             disc_qaly_str_comp,
             disc_qaly_mi_comp,
             disc_qaly_hf_comp,
             
             disc_qaly_noevent_int,
             disc_qaly_esrd_int,
             disc_qaly_str_int,
             disc_qaly_mi_int,
             disc_qaly_hf_int),
    component = rep(components_q, 2)
  )
  
  
  if (analysis == 1){
    list(plot1 = plot1, plot2 = plot2, 
         CEAtable_s = CEAtable_s,
         event_short_com = event_short_com, event_short_int = event_short_int,
         event_lt_com = event_lt_com, event_lt_int = event_lt_int,
         tte_com = tte_com, tte_int = tte_int,
         cost_short = cost_short, qaly_short = qaly_short,
         cost_lt = cost_lt, qaly_lt = qaly_lt,
         CEAresults_s = c("Incremental cost" = disc_incr.cost_s,
                          "Incremental QALY" = disc_incr.qaly_s,
                          "ICER" = ICER_s),
         CCQQresults_s = c("Cost_int" = disc_cost_s_int,
                           "Cost_comp" = disc_cost_s_comp,
                           "QALY_int" = disc_qaly_s_int,
                           "QALY_comp" = disc_qaly_s_comp))
  }else if (analysis == 2){
    return(c("Cost_int" = disc_cost_s_int,
             "Cost_comp" = disc_cost_s_comp,
             "QALY_int" = disc_qaly_s_int,
             "QALY_comp" = disc_qaly_s_comp))
  }else{
    return(c("Incremental cost" = disc_incr.cost_s,
             "Incremental QALY" = disc_incr.qaly_s,
             "ICER" = ICER_s))
  }

}

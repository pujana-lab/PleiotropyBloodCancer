# Figure 1. Forestplot depicting multivariate Cox proportional hazards model for the risk of incident cancer diagnosis within one year from recruitment.

library(dplyr) 
library(tidyr) 
library(ggplot2) 
library(survival)
library(forestmodel)

#UKK Biobank Dataset
df <-"meta/ukbiobankFile.txt" 

# Cox model
outcome <- paste("Surv(tsurv, evsurv) ~ ")

demographics_cont_predictors <- paste(c("pspline(age_at_recruitment)","pspline(body_mass_index)"), collapse= "+")

demographics_cat <- c(
  "strata(sex)",
  "smoking",
  "strata(drinking)",
  "qualifications",
  "strata(nComorbiditiesRec)",
  "townsend",
  "strata(assessment_centre)"
)

genetics <- 

blood_sample_log_linear <- c(
  "erythrocyte_count_log", 
  "platelet_crit_log",  
  "erythrocyte_distr_log", 
  "lymphocyte_count_log", 
  "monocyte_count_log", 
  "neutrophill_count_log",  
  "basophill_count_log", 
  "nucleated_rbd_pct_rec", 
  "reticulocyte_pct_log", 
  "platelet_distr_log",
  "cor_haemoglobin_log", 
  "eosinophill_count_log"  )

demographics_cat_predictors <- paste(demographics_cat, collapse= "+")
genetics_predictors <- paste(genetics, collapse= "+")
blood_sample_log_linearpredictors <- paste(blood_sample_log_linear, collapse= "+")

linearPredictor <- paste(c(demographics_cont_predictors, demographics_cat_predictors, genetics_predictors,"C_reactive_protein_log",blood_sample_log_linearpredictors), collapse= "+")

coxModel <- as.formula( paste(outcome, linearPredictor) )

coxres <-  coxph(coxModel , data = df,  x = TRUE)

# Forestplot
forestmodel::forest_model(coxres, 
                          covariates = c("C_reactive_protein_log", blood_sample_log_linear), 
                          exclude_infinite_cis = TRUE,
                          format_options = forestmodel::forest_model_format_options( 
                          colour = "black",
                          color = NULL,
                          shape = 16, #forma del punt del HR
                          text_size = 3,
                          point_size = 3,
                          banded = TRUE
                          )
)

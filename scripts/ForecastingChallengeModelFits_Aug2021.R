#####
## NSF C-Accel Forecasting Challenge Analysis -- Updated Aug 2021
## Key updates:
  # Using new NOAA data, along with covariates calculated through the sdm_workflow project
  # Based on previous runs, I think we should be able to run all of the groundfish with the random walk model
#####


# Libraries and preliminaries ---------------------------------------------
# _targets.R
# Libraries
# install.packages("devtools")
library(devtools)
# devtools::install_version("Matrix", version = "1.2.8")
library(Matrix)
# devtools::install_version("TMB", "1.7.18")
# devtools::install_github("James-Thorson-NOAA/FishStatsUtils", force = TRUE, upgrade = FALSE)
# devtools::install_github("James-Thorson-NOAA/VAST", ref = "main", force = TRUE, upgrade = FALSE)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(here)
library(splines)
library(parallel)
library(doFuture)
library(tools)
library(gmRi)
os_use<- "unix"

rw_only<- TRUE
docker<- FALSE
if(docker){
  #source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VASTfunctions_AJAedits.R")
  #source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VAST_wrapper_func.R")
  source(here::here("scripts", "VASTfunctions_AJAedits.R"))
  source(here::here("scripts", "VAST_wrapper_func.R"))
  source(here::here("scripts", "fit_model_eff.R"))
} else {
  source(here::here("scripts", "VASTfunctions_AJAedits.R"))
  source(here::here("scripts", "VAST_wrapper_func.R"))
  source("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/R/vast_functions.R")
}

## Key paths
res_data_path<- "~/Box/RES_Data/"
out_folder<- shared.path(os = os_use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results/")

## Key settings
testing<- FALSE

# Loading data ---------------------------------------------
all_dat<- readRDS("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/combined/tidy_mod_data.rds")

# We don't want the DFO data...
nmfs_dat<- all_dat %>%
  dplyr::filter(., SURVEY == "NMFS")

# We do want the nice species names...
nice_spp_names<- read_csv("~/Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/Assesmentfishspecies.csv")

nmfs_dat<- nmfs_dat %>%
  left_join(., nice_spp_names, by = c("NMFS_SVSPP" = "SVSPP"))

# Only 12 groundfish
groundfish<- c("AMERICAN PLAICE", "ATLANTIC COD", "ATLANTIC HALIBUT", "POLLOCK", "ATLANTIC WOLFFISH", "HADDOCK", "OCEAN POUT", "ACADIAN REDFISH", "WHITE HAKE", "WINDOWPANE", "WINTER FLOUNDER", "WITCH FLOUNDER", "YELLOWTAIL FLOUNDER")
nmfs_dat<- nmfs_dat %>%
  dplyr::filter(., COMNAME %in% groundfish)

# Finally year as a factor
nmfs_dat$Year_Cov<- factor(nmfs_dat$EST_YEAR)

# Sample data columns
samp_dat_cols<- c("ID", "EST_YEAR", "Year_Cov", "DECDEG_BEGLAT", "DECDEG_BEGLON", "BIOMASS")
area_swept<- 0.0384

# Covariate data columns
nmfs_dat$Year_Cov<- factor(nmfs_dat$EST_YEAR)
covs<- c("Year_Cov", "Depth", "BT_seasonal", "SST_seasonal")
cov_dat_cols<- c("ID", "EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", covs)

# Model fitting settings ---------------------------------------------
## Extrapolation Grid
# VAST
n_x_use<- 400
region_code<- "northwest_atlantic"
strat_limits<- list('All_areas' = 1:1e5)

## Spatial, Spatio-temporal factors, and then temporal correlation in intercepts and spatio-temporal factors.
fieldconfig_forebase<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_forebase<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2)

## Observation model
obsmodel_use<- c(2, 1)

## Options -- want SD_observation_density
options_use<- c("SD_site_logdensity" = FALSE, "Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "Calculate_Cov_SE" = FALSE, "Calculate_Synchrony" = FALSE, "Calculate_proportion" = FALSE, "SD_observation_density" = TRUE)

## Settings for base model
settings_forebase<- make_settings(n_x = n_x_use, Region = region_code, strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_forebase, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)

## Model formula/fit stuff
gam_degree<- 3
formula_use<- ~ bs(Depth, degree = 3, intercept = FALSE) + bs(SST_seasonal, degree = 3, intercept = FALSE) + bs(BT_seasonal, degree = 3, intercept = FALSE) 
hab_env_coeffs_n<- length(attributes(terms.formula(formula_use))$term.labels)


# Model fitting loop ---------------------------------------------
# Create a nested data frame
dat_nest<- nmfs_dat %>%
  group_by(., SEASON, COMNAME) %>%
  nest()

# Model fitting loop -- outside, overforecast "challenges"...
if(testing){
  fore_challenges<- seq(from = 2017, to = 2019, by = 1)
} else {
  fore_challenges<- seq(from = 2004, to = 2019, by =  1)
}

# Add in the forechallenges piece...
fore_challenges_df<- data.frame("COMNAME" = rep(unique(dat_nest$COMNAME), each = length(fore_challenges)), "fore_challenge" = rep(fore_challenges, length(unique(dat_nest$COMNAME))))

dat_nest<- dat_nest %>%
  left_join(., fore_challenges_df, by = c("COMNAME" = "COMNAME"))

# Detecting cores
cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail-2)

for(i in 2:nrow(dat_nest)){
  
  # Get description for forecast challenge, based on years model will predict to
  fore_dates<- paste(dat_nest$fore_challenge[i], "to2019", sep = "")
  
  # Some descriptors about the specific run...
  season_run<- as.character(dat_nest$SEASON)[i]
  spp_run<- as.character(dat_nest$COMNAME)[i]
  dat_run<- dat_nest$data[[i]]
  
  # Create output folder
  outfolder<- paste(out_folder, paste("VAST", season_run, spp_run, sep = "_"), sep = "")
  if(!file.exists(outfolder)){
    dir.create(outfolder)
  }
  
  # Create sub-folder
  outfile<- paste(outfolder, "/", fore_dates, sep = "")

  # Text file to print progress
  progress_out<- "Starting"
  write(progress_out, file = paste(outfile, "_progress.txt", sep = ""), append = FALSE)
  
  # Create sample data frame and manipulating PRED_TF such that when PRED_TF = 0, the observation is used in the calculation of the likelihood. In contrast, observations with PRED_TF = 1 will not be used in the likelihood, but we will get a model prediction for those observations.
  samp_dat_all<- dat_run %>%
    dplyr::select(., one_of(samp_dat_cols))
  samp_dat_all$PRED_TF<- ifelse(samp_dat_all$EST_YEAR < dat_nest$fore_challenge[i], 0, 1)
  samp_dat_all$Swept<- area_swept
  
  # Finally, need to rename columns...
  samp_dat_all<- samp_dat_all %>%
    rename(., c("Year" = "EST_YEAR", "Year_Cov" = "Year_Cov", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON", "Biomass" = "BIOMASS", "Pred_TF" = "PRED_TF")) %>%
    dplyr::select(., Year, Lat, Lon, Biomass, Swept, Pred_TF) %>%
    data.frame()
  
  # Base, set all years to be the same -- this could be needed if we struggle to fit any of the more complex models
  # samp_dat_simp<- samp_dat_all
  # samp_dat_simp$Year<- rep(min(samp_dat_base$Year), nrow(samp_dat_base))
  
  # Now, covariate data frame. 
  cov_dat<- dat_run %>%
    dplyr::select(., one_of(cov_dat_cols)) %>%
    drop_na(., {{covs}})
  
  # Finally, need to rename columns...
  cov_dat_all<- cov_dat %>%
    rename(., c("Year" = "EST_YEAR", "Year_Cov" = "Year_Cov", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON")) %>%
    dplyr::select(., one_of(c("Year", "Lat", "Lon", covs))) %>%
    data.frame()
 
  # Model fit 
  fit_base<- try(fit_model("settings" = settings_forebase, "Lat_i" = samp_dat_all[, 'Lat'], "Lon_i" = samp_dat_all[, 'Lon'], "t_i" = samp_dat_all[, 'Year'], "c_i" = rep(0, nrow(samp_dat_all)), "b_i" = samp_dat_all[, 'Biomass'], "a_i" = samp_dat_all[, 'Swept'], "PredTF_i" = samp_dat_all[, 'Pred_TF'], "covariate_data" = cov_dat_all, "X1_formula" = formula_use, "X2_formula" = formula_use, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "working_dir" = "~/Box/Mills Lab/Projects/ForecastingChallenge", CompileDir = "~/Box/Mills Lab/Projects/ForecastingChallenge"), silent = TRUE)
 
  if(class(fit_base) == "fit_model" && (max(abs(fit_base$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
    fit_coveff<- get_vast_covariate_effects(vast_fit = fit_base, params_plot = c("Depth", "SST_seasonal", "BT_seasonal"), params_plot_levels = 100, effects_pad_values = c(), nice_category_names = spp_run, out_dir = outfolder)
    
    plot_coveff<- plot_vast_covariate_effects(vast_covariate_effects = fit_coveff, vast_fit = fit_base, nice_category_names = spp_run, out_dir = outfolder)
    
    progress_new<- "Excellent, model with RW on beta and Epsilon worked"
    write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
    fit_base$parameter_estimates$fitted_settings <- "BetaRWEpsRW"
    saveRDS(fit_base, file = paste(outfile, "_betaRWepsRW_modelfit.rds", sep = ""))
  } else {
    progress_new<- "Model with RW on beta and Epsilon failed"
    write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
  }
  # End for each loop over species, seasons, fore_challenge
}


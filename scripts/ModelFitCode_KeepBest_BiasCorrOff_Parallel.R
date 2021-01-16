#####
## NSF C-Accel Forecasting Challenge (12/16/2020)
#####
# Preliminaries ----------------------------------------------------------
# True/False switch as there might be a few packages that need to be installed IF running on DO server
# Libraries and sourcing functions
library(TMB)
library(VAST)
library(tidyverse)
library(here)
library(sf)
library(googledrive)
library(snakecase)
library(rnaturalearth)
library(rnaturalearthdata)
library(conflicted)
library(gmRi)
os.use<- "unix"
library(parallel)
library(doParallel)
library(doFuture)
filter <- dplyr::filter
select <- dplyr::select

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
  source(here::here("scripts", "fit_model_eff.R"))
}

# Load in data -----------------------------------------------------------
if(docker){
  # Fisheries data
  dat<- read.csv(here::here("data", "ForecastingChallengeModelData.csv"))
  # Land shapefile for mapping
  land<- st_read(here::here("data", "NELME_sf.shp"))
} else {
  # Fisheries data
  dat<- read.csv(here::here("data", "ForecastingChallengeModelData.csv"))
  # Land shapefile for mapping
  land<- st_read(here::here("data", "NELME_sf.shp"))
}

# For testing, reduce data set to speed up model fits
testing<- FALSE
if(testing){
  dat<- dat %>%
    filter(., EST_YEAR >= 2012)
}

# Model fitting SEASONALLY INDEPENDENT  --------------------------------------------------------------
### Sample and covariate data columns
# Sample data columns
samp_dat_cols<- c("ID", "EST_YEAR", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", "SUM_BIOMASS")
area_swept<- 0.01

# Covariate data columns
covs<- c("AVGDEPTH", "SODA_BT", "SODA_SST")
cov_dat_cols<- c("ID", "EST_YEAR", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", covs)

### Model fit settings
# Extrapolation Grid
# VAST
n_x_use<- 400
max_dist_from_sample<- 25
grid_dim_km<- c(25, 25)
region_code<- "northwest_atlantic"
strat_limits<- data.frame("STRATA" = unique(dat$STRATUM)[order(unique(dat$STRATUM))])

# User supplied
# Make extrapolation grid from shapefile
if(docker){
  # NELME grid
  nelme_grid<- convert_shapefile(here::here("data", "NELME_sf.shp"), projargs = "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", projargs_for_shapefile = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", grid_dim_km = grid_dim_km, make_plots = FALSE, area_tolerance = 2)
} else {
  # NELME grid
  nelme_grid<- convert_shapefile(here::here("data", "NELME_sf.shp"), projargs = "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", projargs_for_shapefile = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", grid_dim_km = grid_dim_km, make_plots = FALSE, area_tolerance = 2)
}

# Visualize
xlim<- c(-76, -65) 
ylim<- c(35, 45)

extrap_plot<- ggplot() +
  geom_point(data = nelme_grid$extrapolation_grid, aes(x = Lon, y = Lat), size = 0.05) +
  geom_sf(data = land, fill = "#f6f6f6", color = "light gray") +
  coord_sf(xlim, ylim, expand = FALSE) +    
  xlab("") +
  ylab("") +
  theme_bw()
extrap_plot

# Spatial, Spatio-temporal factors, and then temporal correlation in intercepts and spatio-temporal factors.
# For running this in the "keep best" situation, going to want to start somewhere close to the ideal model, rather than building up from the pure environment-only SDM to hopefully reduce the number of models that need to be fit and increase speed and efficiency. With an eye towards using these models to forecast species distributions/abundance, the ideal forecasting model would include both spatial and spatio-temporal correlation to account for unmeasured environmental/biological processes that are either persistent throughout the entire time series or ephemeral. In addition, we are going to want to have some temporal autoregressive structure for both the intercepts and the spatio-temporal correlation as this will facilitate predicting these model components in the future.  
fieldconfig_forebase<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_forebase<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 1, "Epsilon2" = 1)

# Observation model ("Poisson" link function)
obsmodel_use<- c("PosDist" = 4, "EncProbForm" = 1)

# Options -- want SD_observation_density
options_use<- c("SD_site_logdensity" = FALSE, "Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "Calculate_Cov_SE" = FALSE, "Calculate_Synchrony" = FALSE, "Calculate_proportion" = FALSE, "SD_observation_density" = TRUE)

settings_forebase<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_forebase, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)

### Model formula/fit stuff
formula_use<- ~ splines::bs(AVGDEPTH, df = 3, intercept = FALSE) +
  splines::bs(SODA_BT, df = 3, intercept = FALSE) +
  splines::bs(SODA_SST, df = 3, intercept = FALSE)
# formula_use<- "~ AVGDEPTH + SODA_BT"

# Create a nested data frame
dat_nest<- dat %>%
  group_by(., SEASON, COMNAME) %>%
  nest()

# Model fitting loop -- outside, overforecast "challenges"...
if(testing){
  fore_challenges<- seq(from = 2017, to = 2019, by = 1)
} else {
  fore_challenges<- seq(from = 2004, to = 2019, by =  1)
}

# Output folder
if(docker){
  out_folder<- "/home/andrew.allyn@gmail.com/Temp\ Results/"
  if(!file.exists(out_folder)){
    dir.create(out_folder)
  }
} else {
  out_folder<- here::here("temp results/")
  if(!file.exists(out_folder)){
    dir.create(out_folder)
  }
}

# dat_nest<- dat_nest %>% 
#   filter(., COMNAME == "HADDOCK")
dat_nest<- dat_nest %>% 
  filter(., COMNAME == "ATLANTIC COD")

# Add in the forechallenges piece...
fore_challenges_df<- data.frame("COMNAME" = rep(unique(dat_nest$COMNAME), each = length(fore_challenges)), "fore_challenge" = rep(fore_challenges, length(unique(dat_nest$COMNAME))))

dat_nest<- dat_nest %>%
  left_join(., fore_challenges_df, by = c("COMNAME" = "COMNAME"))

# Detecting cores
cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail-1)

# I kept getting errors with the compilation here when trying to run in parallel. Not sure why, but trying something by running a simple example to get the VAST_cpp and VAST_so files. Then copying those into the folders as needed.
## Setup really simple example for testing
example <- load_example( data_set="EBS_pollock" )
dat <- subset(example$sampling_data, Year==2013)
settings <- make_settings(n_x=100, Region=example$Region,
                          purpose="index2", bias.correct=FALSE )
settings$FieldConfig[1:2, 1:2] <- 0

## Default settings to compile the latest CPP in the R package
## directory 'executables'
fit = fit_model(settings=settings,
                working_dir=here::here(),
                Lat_i=dat$Lat, Lon_i=dat$Lon, t_i=dat$Year,
                b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2, run_model = FALSE)
future:::ClusterRegistry("stop")

### Parallel loop!
all<- dat_nest
dat_nest<- all

# Cluster stuff
cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail-1)
start_time<- Sys.time()

foreach(i = 1:nrow(dat_nest)) %dopar% {
  
  library(VAST)
  
  #for(i in 1:nrow(dat_nest)){
  
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
    drop_na(., {{covs}}) %>%
    select(., one_of(samp_dat_cols))
  samp_dat_all$PRED_TF<- ifelse(samp_dat_all$EST_YEAR < dat_nest$fore_challenge[i], 0, 1)
  
  # Finally, need to rename columns...
  samp_dat_all<- samp_dat_all %>%
    rename(., c("Year" = "EST_YEAR", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON", "Catch_KG" = "SUM_BIOMASS")) %>%
    data.frame()
  
  # Base, set all years to be the same -- this could be needed if we struggle to fit any of the more complex models
  samp_dat_simp<- samp_dat_all
  samp_dat_simp$Year<- rep(min(samp_dat_simp$Year), nrow(samp_dat_simp))
  
  # Now, covariate data frame. 
  cov_dat<- dat_run %>%
    select(., one_of(cov_dat_cols)) %>%
    drop_na(., {{covs}})
  
  # Rescale covariates to SD of 1
  cov_dat<- cov_dat %>%
    mutate_at(., {{covs}}, vast_scale_func, type = "AJA")
  
  # Finally, need to rename columns...
  cov_dat_all<- cov_dat %>%
    rename(., c("Year" = "EST_YEAR", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON")) %>%
    select(., one_of(c("Year", "Lat", "Lon", covs))) %>%
    data.frame()
  
  if(explore){
    summary(cov_dat_all)
    summary(samp_dat_all)
    hist(samp_dat_all$Catch_KG)
    
    # Extreme values?
    samp_dat_all<- samp_dat_all %>%
      filter(., Catch_KG <= 200)
    cov_dat_all<- cov_dat_all
  }
  
  # Base, set all years to be the same -- this could be needed if we struggle to fit any of the more complex models
  cov_dat_simp<- cov_dat_all
  cov_dat_simp$Year<- rep(min(cov_dat_simp$Year), nrow(cov_dat_simp))
  
  # Model fit wrapper function
  fit_base<- try(fit_model("settings" = settings_forebase,
                               #bias.correct.control = list(nsplit = 5),
                               # Spatial info
                               observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                               # Model info
                               "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
  
  # If our basic forecasting model worked, let's add in some complexity...
  if(class(fit_base) == "fit_model" && (max(abs(fit_base$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
    # Append progress file
    progress_new<- "Base model has passed fit checks, trying to add complexity"
    write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
   
    # Try adding on autoregressive structure to the intercept and remake settings
    rhoconfig_betaAR1<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 1, "Epsilon2" = 1)
    settings_betaAR1<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_betaAR1, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
    
    # Refit the model
    fit_betaAR1<- try(fit_model("settings" = settings_betaAR1,
                                    #bias.correct.control = list(nsplit = 5),
                                    # Spatial info
                                    observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                    # Model info
                                    "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
    
    # If that worked, let's try AR1 on spatio-temporal structure too...
    if(class(fit_betaAR1) == "fit_model" && (max(abs(fit_betaAR1$parameter_estimates$diagnostics$final_gradient)) <= 0.01)) {
      # Append progress file
      progress_new<- "Adding AR1 to beta converged, now trying to add temporal structure to spatio-temporal variability too"
      write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
      
      rhoconfig_bothAR1<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 4, "Epsilon2" = 4)
      settings_bothAR1<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_bothAR1, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
      
      # Refit the model
      fit_bothAR1<- try(fit_model("settings" = settings_bothAR1,
                                      #bias.correct.control = list(nsplit = 5),
                                      # Spatial info
                                      observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                      # Model info
                                      "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
      
      # If that worked, we are all done. Save fit_bothAR1 and make note of the model settings
      if(class(fit_bothAR1) == "fit_model" && max(abs(fit_bothAR1$parameter_estimates$diagnostics$final_gradient)) <= 0.01){
        # Append progress file
        progress_new<- "Excellent, model with AR1 structure on both beta and spatio-temporal variation converged"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        fit_bothAR1$parameter_estimates$fitted_settings <- "BothAR1"
        saveRDS(fit_bothAR1, file = paste(outfile, "_bothAR1_modelfit.rds", sep = ""))
      } 
      
      # If that didn't work, next option is to try a RW model for spatio-temporal variation...
      if(class(fit_bothAR1) == "try-error" || (!class(fit_bothAR1) == "try-error" & max(abs(fit_bothAR1$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        # Append progress file
        progress_new<- "Unable to fit model with AR1 on spatio-temporal variation, trying RW instead"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        
        # Adjusting settings
        rhoconfig_betaAR1stRW<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 2, "Epsilon2" = 2)
        settings_betaAR1stRW<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_betaAR1stRW, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
        
        # Refit the model
        fit_betaAR1stRW<- try(fit_model("settings" = settings_betaAR1stRW,
                                            #bias.correct.control = list(nsplit = 5),
                                            # Spatial info
                                            observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                            # Model info
                                            "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
        
        # Alright, after running that (even if it failed) we are at the end of this subset. If it worked, we save the betaAR1stRW. 
        if(class(fit_betaAR1stRW) == "fit_model" && max(abs(fit_betaAR1stRW$parameter_estimates$diagnostics$final_gradient)) <= 0.01){
          # Append progress file
          progress_new<- "Excellent, model with AR1 on beta and RW on spatio-temporal variability converged"
          write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
          fit_betaAR1stRW$parameter_estimates$fitted_settings <- "BetaAR1stRW"
          saveRDS(fit_betaAR1stRW, file = paste(outfile, "_betaAR1stRW_modelfit.rds", sep = "")) 
        }
        
        # If not, we go back and save the one with just temporal structure on beta
        if(class(fit_betaAR1stRW) == "try-error" || (!class(fit_betaAR1stRW) == "try-error" & max(abs(fit_betaAR1stRW$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
          # Append progress file
          progress_new<- "Despite trying to add temporal structure to spatio-temporal variability, none of the models (AR1 or RW converged), so just using model with temporal structure on beta"
          write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
          fit_betaAR1$parameter_estimates$fitted_settings <- "BetaAR1"
          saveRDS(fit_betaAR1, file = paste(outfile, "_betaAR1_modelfit.rds", sep = ""))
        } 
      }
      # End trying AR1 options for beta
    } 
    
    # If that didn't work, then we can try the RW on beta and either save that model OR if that one doesn't work, then just keep the basic forecast model
    if(class(fit_betaAR1) == "try-error" || (!class(fit_betaAR1) == "try-error" & max(abs(fit_betaAR1$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
      
      # Append progress file
      progress_new<- "Unable to fit AR1 to beta, trying simpler RW instead"
      write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
      
      # Try RW
      rhoconfig_betaRW<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 1, "Epsilon2" = 1)
      settings_betaRW<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_betaRW, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
      
      # Refit the model
      fit_betaRW<- try(fit_model("settings" = settings_betaRW, 
                                     #bias.correct.control = list(nsplit = 5),
                                     # Spatial info
                                     observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                     # Model info
                                     "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
      
      # If that worked, save it
      if(class(fit_betaRW) == "fit_model" && (max(abs(fit_betaRW$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
        # Append progress file
        progress_new<- "Excellent, model with RW on beta worked"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        fit_betaRW$parameter_estimates$fitted_settings <- "BetaRW"
        saveRDS(fit_betaRW, file = paste(outfile, "_betaRW_modelfit.rds", sep = ""))
      }
      
      # If it didn't work, back to the basic model
      if(class(fit_betaRW) == "try-error" || (!class(fit_betaRW) == "try-error" & max(abs(fit_betaRW$parameter_estimates$diagnostics$final_gradient)) > 0.01)) {
        # Append progress file
        progress_new<- "Unable to fit a model with temporal structure on beta, so just going to use the basic forecast model without temporal structure on beta or on spatio-temporal variability"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        fit_base$parameter_estimates$fitted_settings <- "Base"
        saveRDS(fit_base, file = paste(outfile, "_base_modelfit.rds", sep = ""))
      }
      # End trying RW on beta
    }
    # End trying to add complexity to the base model
  }
  
  # If an issue with the base model
  if(class(fit_base) == "try-error" || (!class(fit_base) == "try-error" & max(abs(fit_base$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
    # Append progress file
    progress_new<- "Convergence issues with just the basic forecast model, trying now to turn on/off spatial or spatio-temporal variability as needed"
    write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
    
    # We were unable to even fit the base forecasting model. This (most likely) means either a problem with estimating the spatial or the spatio-temporal variability. 
    # Turn off spatio-temporal variability
    fieldconfig_nost<- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 0)
    rhoconfig_nost_betarw<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 0, "Epsilon2" = 0)
    
    settings_nost_betarw<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_nost, RhoConfig = rhoconfig_nost_betarw, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
    
    # Model fit wrapper function
    fit_nost_betarw<- try(fit_model("settings" = settings_nost_betarw,
                                        #bias.correct.control = list(nsplit = 5),
                                        # Spatial info
                                        observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                        # Model info
                                        "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
    
    # If that worked, save result
    if(class(fit_nost_betarw) == "fit_model" && (max(abs(fit_nost_betarw$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
      # Append progress file
      progress_new<- "Model with just spatial variability and betaRW converged"
      write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
      fit_nost_betarw$parameter_estimates$fitted_settings <- "Nost_betarw"
      saveRDS(fit_nost_betarw, file = paste(outfile, "_nost_betarw_modelfit.rds", sep = ""))
    } 
    
    # If it didn't, try turning off betaRW
    if(class(fit_nost_betarw) == "try-error" || (!class(fit_nost_betarw) == "try-error" & max(abs(fit_nost_betarw$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
      # Append progress file
      progress_new<- "Convergence issues with spatial and betaRW, turning off betaRW"
      write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
      
      # We were unable to even fit the base forecasting model. This (most likely) means either a problem with estimating the spatial or the spatio-temporal variability. 
      
      # Turn off spatio-temporal variability
      fieldconfig_nost<- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 0)
      rhoconfig_nost<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 0, "Epsilon2" = 0)
      
      settings_nost<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_nost, RhoConfig = rhoconfig_nost, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
      
      # Model fit wrapper function
      fit_nost<- try(fit_model("settings" = settings_nost,
                                   #bias.correct.control = list(nsplit = 5),
                                   # Spatial info
                                   observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                   # Model info
                                   "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
      
      # If that worked, save result
      if(class(fit_nost) == "fit_model" && (max(abs(fit_nost$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
        # Append progress file
        progress_new<- "Model with just spatial variability and betaRW converged"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        fit_nost$parameter_estimates$fitted_settings <- "Nost"
        saveRDS(fit_nost, file = paste(outfile, "_nost_modelfit.rds", sep = ""))
      }
      
      # If that still didn't work, turn off spatial variability
      if(class(fit_nost) == "try-error" || (!class(fit_nost) == "try-error" & max(abs(fit_nost$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        # Append progress file
        progress_new<- "Turning off spatio-temporal variability did not solve the convergence issues, now trying model with spatio-temporal and without spatial"
        write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
        
        # Turn off spatial variability
        fieldconfig_nosp<- c("Omega1" = 0, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 1)
        rhoconfig_nosp<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 1, "Epsilon2" = 1)
        
        settings_nosp<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_nosp, RhoConfig = rhoconfig_nosp, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
        
        # Model fit wrapper function
        fit_nosp<- try(fit_model("settings" = settings_nosp,
                                     #bias.correct.control = list(nsplit = 5),
                                     # Spatial info
                                     observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                     # Model info
                                     "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
        
        # If that worked, all done...
        if(class(fit_nosp) == "fit_model" && (max(abs(fit_nosp$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
          # Append progress file
          progress_new<- "Model without spatial variability converged"
          write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
          fit_nosp$parameter_estimates$fitted_settings <- "Nosp"
          saveRDS(fit_nosp, file = paste(outfile, "_nosp_modelfit.rds", sep = ""))
        } 
        
        # If not...turn everything off
        if(class(fit_nosp) == "try-error" || (!class(fit_nosp) == "try-error" & max(abs(fit_nosp$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
          # Append progress file
          progress_new<- "Model still not converging, going to have to turn off both spatial and spatio-temporal variability"
          write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
          
          # Turn off everything
          fieldconfig_simp<- c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
          rhoconfig_simp<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 0, "Epsilon2" = 0)
          
          settings_simp<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_simp, RhoConfig = rhoconfig_simp, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)
          
          # Model fit wrapper function
          fit_simp<- try(fit_model("settings" = settings_simp,
                                       #bias.correct.control = list(nsplit = 2),
                                       # Spatial info
                                       observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                       # Model info
                                       "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = FALSE), silent = TRUE)
          
          # Did that work?
          if(class(fit_simp) == "try-error" || (!class(fit_simp) == "try-error" & max(abs(fit_simp$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
            # Append progress file noting failure
            progress_new<- "Simple model did not converge"
            write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
          }
          
          if(class(fit_simp) == "fit_model" && (max(abs(fit_simp$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
            # Append progress file
            progress_new<- "Simple model converged"
            write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
            fit_simp$parameter_estimates$fitted_settings <- "Simp"
            saveRDS(fit_simp, file = paste(outfile, "_simp_modelfit.rds", sep = ""))
          }
        }
      }
    }
    # End trying even simpler models than the basic forecast model
  }
  # End for each loop over species, seasons, fore_challenge
}

# Cluster clean up
future:::ClusterRegistry("stop")
gc()

elapsed_time<- Sys.time() - start_time

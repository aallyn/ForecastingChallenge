#####
## NSF C-Accel Forecasting Challenge (09/15/2020)
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
library(future.apply)
library(conflicted)
library(gmRi)
filter <- dplyr::filter
select <- dplyr::select

docker<- FALSE
if(docker){
  source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VASTfunctions_AJAedits.R")
  source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VAST_wrapper_func.R")
} else {
  source(here::here("scripts", "VASTfunctions_AJAedits.R"))
  source(here::here("scripts", "VAST_wrapper_func.R"))
}

# Load in data -----------------------------------------------------------
if(docker){
  # Fisheries data
  dat<- read.csv(paste("/home/andrew.allyn@gmail.com/ForecastingChallenge/Data/", 
                         "ForecastingChallengeModelData.csv", sep = ""))
  # Land shapefile for mapping
  land<- st_read("/home/andrew.allyn@gmail.com/Shapefiles/ne_50m_land/ne_50m_land.shp")
} else {
  # Fisheries data
  fish_dat_path<- shared.path(os = os.use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/data/")
  dat<- read.csv(paste(fish_dat_path, "ForecastingChallengeModelData.csv", sep = ""))
  # Land shapefile for mapping
  res_dat_path<- shared.path(os = os.use, group = "RES Data")
  land<- st_read(paste(res_dat_path, "Shapefiles/ne_50m_land/ne_50m_land.shp", sep = "")) 
}

# For testing, reduce data set to speed up model fits
testing<- TRUE
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
n_x_use<- 100
max_dist_from_sample<- 25
grid_dim_km<- c(50, 50)
region_code<- "northwest_atlantic"
strat_limits<- data.frame("STRATA" = unique(dat$STRATUM)[order(unique(dat$STRATUM))])

# User supplied
# Make extrapolation grid from shapefile
if(docker){
  # NELME grid
  nelme_grid<- convert_shapefile("/home/andrew.allyn@gmail.com/Shapefiles/NELME_regions/NELME_sf.shp", projargs = NULL, projargs_for_shapefile = NULL, grid_dim_km = grid_dim_km, make_plots = FALSE, area_tolerance = 0.05)
} else {
  # NELME grid
  nelme_grid<- convert_shapefile(paste(res_dat_path, "Shapefiles/NELME_regions/NELME_sf.shp", sep = ""), projargs = "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", projargs_for_shapefile = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", grid_dim_km = grid_dim_km, make_plots = FALSE, area_tolerance = 2)
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
# First try -- the basic model. Ideally, want this model to as closely as possible resemble what we would get using mgcv::gam(covariates). I have gone back and forth on this ended up treating the intercept as a fixed effect and for the baseline data setting all of the years to be equal. I also thought about estimating the average occurrence across time steps (e.g., the intercept) as a random, CONSTANT, effect, where BetaN = 3, but that still seems different than the "traditional" GAM environment-only model with DEPTH, BT and SST. 
fieldconfig_base<- c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
rhoconfig_base<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 0, "Epsilon2" = 0)

# Observation model ("Poisson" link function)
obsmodel_use<- c("PosDist" = 4, "EncProbForm" = 1)

# Options -- want SD_observation_density
options_use<- c("SD_site_logdensity" = FALSE, "Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "Calculate_Cov_SE" = FALSE, "Calculate_Synchrony" = FALSE, "Calculate_proportion" = FALSE, "SD_observation_density" = TRUE)

settings_base<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_base, RhoConfig = rhoconfig_base, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = FALSE, bias.correct = FALSE)

# Persistent spatial correlation settings
fieldconfig_sp<- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 0)
settings_sp<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_sp, RhoConfig = rhoconfig_base, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

# Persistent spatial and spatio-temporal correlation settings. Does this one make sense? Still have one global intercept for each linear predictor (average catch rate/probability of presence across the entire domain and time series)?
fieldconfig_spst<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
settings_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_spst, RhoConfig = rhoconfig_base, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

# Fit intercept as random effect with RW structure or AR1 structure, persistent spatial and spatio-temporal correlation settings. Does this one make sense? Still have one global intercept for each linear predictor (average catch rate/probability of presence across the entire domain and time series)?
fieldconfig_rwint_spst<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_rwint_spst<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 0, "Epsilon2" = 0)
settings_rwint_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_rwint_spst, RhoConfig = rhoconfig_rwint_spst, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

fieldconfig_AR1int_spst<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_AR1int_spst<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 0, "Epsilon2" = 0)
settings_AR1int_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_AR1int_spst, RhoConfig = rhoconfig_AR1int_spst, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

# Now also add temporal structure to spatio-temporal variability...
rhoconfig_rwintst_spst<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2)
settings_rwintst_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_rwint_spst, RhoConfig = rhoconfig_rwintst_spst, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

rhoconfig_AR1intst_spst<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 4, "Epsilon2" = 4)
settings_AR1intst_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_rwint_spst, RhoConfig = rhoconfig_AR1intst_spst, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

rhoconfig_RWintAR1st_spst<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 4, "Epsilon2" = 4)
settings_RWintAR1st_spst<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_rwint_spst, RhoConfig = rhoconfig_RWintAR1st_spst, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = TRUE)

### Model formula/fit stuff
formula_use<- ~ splines::bs(AVGDEPTH, knots = NULL, intercept = FALSE) +
  splines::bs(SODA_BT, knots = NULL, intercept = FALSE) +
  splines::bs(SODA_SST, knots = NULL, intercept = FALSE)

# Create a nested data frame
dat_nest<- dat %>%
  group_by(., SEASON, COMNAME) %>%
  nest()

# Model fitting loop -- outside, overforecast "challenges"...
if(testing){
  fore_challenges<- seq(from = 2015, to = 2019, by = 1)
} else {
  fore_challenges<- seq(from = 2004, to = 2019, by =  1)
}

# Output folder
out_folder<- shared.path(os = os.use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results/")

for(h in 1:forecast_challenges){
  
  # Get description for forecast challenge, based on years model will predict to
  fore_dates<- paste(fore_challenges[h], "to2019", sep = "")
  
  ## Inner loop -- each species-season, in parellel
  foreach(i = nrow(dat_nest)) %dopar% {
  #for(i in 1:nrow(dat_nest)){
    
    # Some descriptors about the specific run...
    season_run<- as.character(dat_nest$SEASON)[i]
    spp_run<- as.character(dat_nest$COMNAME)[i]
    dat_run<- dat_nest$data[[i]]
    
    # Create output file
    outfile<- paste(out_folder, paste("VAST", season_run, spp_run, sep = "_"), sep = "")
    if(!file.exists(outfile)){
      dir.create(outfile)
    }
    
    # Create sample data frame and manipulating PRED_TF such that when PRED_TF = 0, the observation is used in the calculation of the likelihood. In contrast, observations with PRED_TF = 1 will not be used in the likelihood, but we will get a model prediction for those observations.
    samp_dat_all<- dat_run %>%
      drop_na(., {{covs}}) %>%
      select(., one_of(samp_dat_cols))
    samp_dat_all$PRED_TF<- ifelse(samp_dat_all$EST_YEAR < fore_challenges[h], 0, 1)
   
    # Finally, need to rename columns...
    samp_dat_all<- samp_dat_all %>%
      rename(., c("Year" = "EST_YEAR", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON", "Catch_KG" = "SUM_BIOMASS")) %>%
      data.frame()
    
    # Base, set all years to be the same
    samp_dat_base<- samp_dat_all
    samp_dat_base$Year<- rep(min(samp_dat_base$Year), nrow(samp_dat_base))
    
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
    
    # Base, set all years to be the same
    cov_dat_base<- cov_dat_all
    cov_dat_base$Year<- rep(min(cov_dat_base$Year), nrow(cov_dat_base))
    
    # Model fit wrapper function
    fit_base<- fit_model("settings" = settings_base,
                          # Spatial info
                          observations_LL = cbind("Lat" = samp_dat_base[,'Lat'], "Lon" = samp_dat_base[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                          # Model info
                          "Lat_i" = samp_dat_base[,'Lat'], "Lon_i" = samp_dat_base[,'Lon'], "t_i" = as.vector(samp_dat_base[,'Year']), "b_i" = samp_dat_base[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_base)), "v_i" = rep(0, nrow(samp_dat_base)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_base)), covariate_data = cov_dat_base, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_base[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_base, file = paste(outfile, "/", fore_dates, "_base_modelfit.rds", sep = ""))
    
    # Turn on persistent spatial correlation
    fit_sp<- fit_model("settings" = settings_sp,
                         # Spatial info
                         observations_LL = cbind("Lat" = samp_dat_base[,'Lat'], "Lon" = samp_dat_base[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE,
                         # Model info
                         "Lat_i" = samp_dat_base[,'Lat'], "Lon_i" = samp_dat_base[,'Lon'], "t_i" = as.vector(samp_dat_base[,'Year']), "b_i" = samp_dat_base[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_base)), "v_i" = rep(0, nrow(samp_dat_base)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_base)), covariate_data = cov_dat_base, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_base[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_sp, file = paste(outfile, "/", fore_dates, "_sp_modelfit.rds", sep = ""))
    
    
    # Turn on persistent spatial and ephemeral spatio-temporal correlation
    fit_spst<- fit_model("settings" = settings_spst,
                       # Spatial info
                       observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                       # Model info
                       "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_spst, file = paste(outfile, "/", fore_dates, "_spst_modelfit.rds", sep = ""))
    
    # Turn on RW intercept and persistent spatial and ephemeral spatio-temporal correlation
    fit_rwint_spst<- fit_model("settings" = settings_rwint_spst,
                               # Spatial info
                               observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                               # Model info
                               "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_rwint_spst, file = paste(outfile, "/", fore_dates, "_rwint_spst_modelfit.rds", sep = ""))

    # Turn on RE AR1 intercept and persistent spatial and ephemeral spatio-temporal correlation
    fit_AR1int_spst<- fit_model("settings" = settings_AR1int_spst,
                         # Spatial info
                         observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                         # Model info
                         "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_AR1int_spst, file = paste(outfile, "/", fore_dates, "_AR1int_spst_modelfit.rds", sep = ""))
    
    # Turn on RE RW intercept AND spatio-temporal correlation, and persistent spatial and ephemeral spatio-temporal correlation
    fit_rwintst_spst<- fit_model("settings" = settings_rwintst_spst,
                                # Spatial info
                                observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                                # Model info
                                "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_all[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = paste(outfile, "/", sep = ""), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
    save(fit_rwintst_spst, file = paste(outfile, "/", fore_dates, "_rwintst_spst_modelfit.rds", sep = ""))
  }
}

# Model fitting INTRA ANNUAL  --------------------------------------------------------------
# Ideally, we have this fitted with the intra-annual approach (e.g., spring 2000 influences fall 2000, then spring 2001, etc..) That said, a lot of work to do that. While I have some code that I *think* works (replicates JT paper results), going to try to keep this all with JT's code/implementation.
glimpse(dat)

# Sample data columns
samp.dat.cols<- c("ID", "EST_YEAR", "SEASON", "SEASON_YEAR", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", "SUM_BIOMASS")
area.swept<- 0.01

# Covariate data columns
covs.use<- c("AVGDEPTH")
cov.dat.cols<- c("ID", "EST_YEAR", "SEASON", "SEASON_YEAR", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SEASON_NUM")

### Model fit settings
# Extrapolation Grid
knots.use<- 400
maximum.dist.from.sample<- 10
grid.dim.km<- c(5, 5)
region.code<- "northwest_atlantic"
strat.limits<- unique(dat$STRATUM)[order(unique(dat$STRATUM))]

### Setting up seasons and years
# Set of seasons and years
season.set<- c("SPRING", "FALL")
year.set<- sort(unique(dat$EST_YEAR))

# Time levels and adding a "season_year" factor
# First, make sure we have all unique combinations of season and year. This doesn't matter for us, but would if we wanted to predict to seasons NOT sampled in a given year. 
seasonyear.grid<- expand.grid("season" = season.set, "year" = year.set)

# Combine "season" and "year" into a "season_year"
t.levels<- apply(seasonyear.grid, MARGIN = 1, FUN = paste, collapse = "_")

# Make that a numeric value with year and decimal for the "season" (.0 or 0.5 for "fall").
t.labels = round(seasonyear.grid[, 'year'] + (as.numeric(factor(seasonyear.grid[, 'season'], levels = season.set))-1)/length(season.set), digits = 1)

# Now, do the same but for our actual data
seasonyear.i<- apply(dat[, c("SEASON", "EST_YEAR")], MARGIN = 1, FUN = paste, collapse = "_")
seasonyear.i<- factor(seasonyear.i, levels = t.levels)
dat<- cbind(dat, "SEASON_YEAR" = seasonyear.i)
dat$SEASON = factor(dat$SEASON, levels = season.set)

### Now, just need the sampling and covariate data...
# Create a nested data frame
dat.nest<- dat %>%
  group_by(., COMNAME) %>%
  nest()

# Some descriptors about the specific run...
spp.run<- as.character(dat.nest$COMNAME)[1]
dat.run<- dat.nest$data[[1]]

# Create output file
OutFile<- here::here(paste("VAST_IntraAnnual", spp.run, sep = "_"))
if(!file.exists(OutFile)){
  dir.create(OutFile)
}

# Create sample data frame. 
samp.dat<- dat.run %>%
  select(., one_of(samp.dat.cols))

# Define the PRED_TF vector, if 0, then observation is used in likelihood, and if 1, then only in predictive probability...
samp.dat$PRED_TF<- ifelse(samp.dat$EST_YEAR < fore.challenges[h], 0, 1)

# Finally, need to rename columns...
samp.dat<- samp.dat %>%
  rename(., c("Year" = "EST_YEAR", "Season" = "SEASON", "Season_Year" = "SEASON_YEAR", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON", "Catch_KG" = "SUM_BIOMASS")) %>%
  data.frame()

# Now, covariate data frame. 
cov.dat<- dat.run %>%
  select(., one_of(cov.dat.cols)) %>%
  drop_na(., {{covs.use}})

# Rescale covariates to SD of 1
cov.dat<- cov.dat %>%
  mutate_at(., {{covs.use}}, vast_scale_func, type = "AJA")

# Finally, need to rename columns...
cov.dat<- cov.dat %>%
  rename(., c("Year" = "EST_YEAR", "Season" = "SEASON", "Season_Year" = "SEASON_YEAR", "Lat" = "DECDEG_BEGLAT", "Lon" = "DECDEG_BEGLON")) %>%
  select(., one_of(c("Year", "Season", "Season_Year", "Lat", "Lon", covs.use))) %>%
  data.frame()

# Error: Year 0 not found in 'covariate_data' popped up. I think this is because of how time is specified. With t_i as.numeric(samp.dat[,"Season_Year"])-1.
cov.dat$Year<- as.numeric(cov.dat[, "Season_Year"])-1

### Field/Rho configurations. Field shouldn't change from the seasonal-independent model. Rho, though, may want to go to having both as AR1 processes? 
FieldConfig.JT<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
RhoConfig.JT<- c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4)
XiConfig.JT<- c("Xi1_season" = 3, "Xi1_year" = 1, "Xi2_season" = 1, "Xi2_year" = 1) # 0: Off; 1 = Linear; 2 = Spatially-varying, zero centered; 3 = Spatially-varying, linear effect. Why are these set up the way they are, would think we'd want these all to be set at 2 for spatially varying zero-centered linear effect?

FieldConfig.all1<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
RhoConfig.allAR1<- c("Beta1" = 4, "Beta2" = 4, "Epsilon1" = 4, "Epsilon2" = 4)
XiConfig.AJA<- c("Xi1_season" = 2, "Xi1_year" = 1, "Xi2_season" = 2, "Xi2_year" = 1)

# Capture settings
input.settings<- list("FieldConfig" = FieldConfig.all1, "RhoConfig" = RhoConfig.allAR1, "XiConfig" = XiConfig.AJA, "Options" = options.use)
capture.output(input.settings, file = paste0(OutFile, "Input_settings.txt"))

# Make settings
settings.use<- make_settings(n_x = knots.use, Region = region.code, strata.limits = strat.limits, purpose = "index2", FieldConfig = FieldConfig.all1, RhoConfig = RhoConfig.allAR1, ObsModel = ObsModel.use, Options = options.use, use_anisotropy = TRUE, fine_scale = TRUE, bias.correct = FALSE)
season.names = season.labels = season.set

# Build initial model
Design<- "Both"
Use_REML<- FALSE
fit.init<- fit_model("settings" = settings.use, strata.limits = strat.limits, "Lat_i" = samp.dat[,'Lat'], "Lon_i" = samp.dat[,'Lon'], "t_i" = as.numeric(samp.dat[,"Season_Year"])-1, "b_i" = samp.dat[,'Catch_KG'], "c_iz" = rep(0, nrow(samp.dat)), "v_i" = rep(0, nrow(samp.dat)), "Q_ik" = NULL, "a_i" = rep(area.swept, nrow(samp.dat)), "PredTF_i" = samp.dat[,'PRED_TF'], covariate_data = cov.dat, formula = formula.use, "observations_LL" = cbind("Lat" = samp.dat[,'Lat'], "Lon" = samp.dat[, 'Lon']), "maximum_distance_from_sample" = maximum.dist.from.sample, "grid_dim_km" = grid.dim.km, "working_dir" = OutFile, "Use_REML" = Use_REML, "run_model" = TRUE)

# Check out initial settings of key objects
str(fit$data_list$Xconfig_zcp)
str(fit$data_list$X_itp)
summary(fit$data_list$X_itp[,,1])

### Make covariate matrix -- SEASON
# For observations
season.i<- match(samp.dat[, 'Season'], season.set)
Xseason.ip<- ThorsonUtilities::vector_to_design_matrix(season.i)
Xseason.itp<- aperm(Xseason.ip %o% rep(1, fit.init$data_list$n_t), c(1,3,2))

# For grid
season.t<- match(seasonyear.grid[, 'Season'], season.set)
Xseason.tp<- ThorsonUtilities::vector_to_design_matrix(season.t)
Xseason.gtp = rep(1, fit.init$data_list$n_g) %o% Xseason.tp

#### Make covariate matrix -- YEAR
# For observations
year.i<- match(samp.dat[, 'Year'], year.set)
Xyear.ip<- ThorsonUtilities::vector_to_design_matrix(year.i)
Xyear.itp<- aperm(Xyear.ip %o% rep(1, fit.init$data_list$n_t), c(1, 3, 2))

# For grid
year.t<- match(seasonyear.grid[, 'Year'], year.set)
Xyear.tp<- ThorsonUtilities::vector_to_design_matrix(year.t)
Xyear.gtp<- rep(1, fit.init$data_list$n_g) %o% Xyear.tp

# Define season-year main effects as zero-centered. What were they before??
season.Xconfig.zcp<- array(NA, dim = c(2, 1, ncol(Xseason.ip)), dimnames = list(NULL, NULL, colnames(Xseason.ip)))
# Impose identifiability constraints on year effects
year.Xconfig.zcp = array(NA, dim = c(2, 1, ncol(Xyear.ip)), dimnames = list(NULL, NULL, colnames(Xyear.ip)))

# Year and season effects
XiConfig
season.Xconfig.zcp[1, ,]<- XiConfig["Xi1_season"]
year.Xconfig.zcp[1, ,]<- XiConfig["Xi1_year"]
season.Xconfig.zcp[2, ,]<- XiConfig["Xi2_season"]
year_Xconfig_zcp[2, ,]<- XiConfig["Xi2_year"]

# If Design="Both" and no RhoConfig structure, then drop two DF (corner constraint for season and year main effects, one for corner constraint of main effects + interactions)
# If using RhoConfig structure, drop one less DF. What is this doing? Ensuring that the season and year effect for the first time step is 0?
# If Design!="Both, drop one less DF
# SOLUTION:  Drop one level from each Season and Year effect
FUN = function(num){
  if(num == 1) return(0)
  if(num == 3) return(2)
  return(num)
}

# Adjusting our original settings.
season.Xconfig.zcp[1, , 1]<- FUN(season.Xconfig.zcp[1, , 1])
season.Xconfig.zcp[2, , 1]<- FUN(season.Xconfig.zcp[2, , 1])
year.Xconfig.zcp[1, , 1]<- FUN(year.Xconfig.zcp[1, , 1])
year.Xconfig.zcp[2, , 1]<- FUN(year.Xconfig.zcp[2, , 1])

# The end result here makes some sense for year, where the effect of year 1 is 0 for both linear predictors, which is the boundary constraint piece. I don't understand the result for season, where now we have...first season spatially-varying zero-centered effect for lin pred 1 and no effect of season for lin pred 2, and then second season we have a spatially-varying linear effect for lin pred 1 and a linear effect for lin pred 2?

# Combine season and year matrices
Xconfig.zcp = X.itp = X.gtp = NULL
if(Design %in% c("Season", "Both")){
  X.itp<- abind::abind(X.itp, Xseason.itp, along = 3)
  X.gtp<- abind::abind(X.gtp, Xseason.gtp, along = 3)
  Xconfig.zcp<- abind::abind(Xconfig.zcp, season.Xconfig_zcp)
}
if(Design %in% c("Year", "Both")){
  X.itp<- abind::abind(X.itp, Xyear.itp, along = 3)
  X.gtp<- abind::abind(X.gtp, Xyear.gtp, along = 3)
  Xconfig.zcp<- abind::abind(Xconfig.zcp, year.Xconfig.zcp)
}

# A lot of changes there. For X_itp, we now have Nobs * Nseasonyears * Nyears + season levels (i.e., covariates) and then for Xconfig_zcp we have Nlinpreds * Ncategories * Nyears+season levels (22).
fit.intra<- fit_model("settings" = settings.use, strata.limits = strat.limits, "Lat_i" = samp.dat[,'Lat'], "Lon_i" = samp.dat[,'Lon'], "t_i" = as.numeric(samp.dat[,"Season_Year"])-1, "b_i" = samp.dat[,'Catch_KG'], "c_iz" = rep(0, nrow(samp.dat)), "v_i" = rep(0, nrow(samp.dat)), "Q_ik" = NULL, "a_i" = rep(area.swept, nrow(samp.dat)), "PredTF_i" = samp.dat[,'PRED_TF'], covariate_data = cov.dat, formula = formula.use, "observations_LL" = cbind("Lat" = samp.dat[,'Lat'], "Lon" = samp.dat[, 'Lon']), "maximum_distance_from_sample" = maximum.dist.from.sample, "grid_dim_km" = grid.dim.km, "working_dir" = OutFile, "Use_REML" = Use_REML, "X_itp" = X.itp, "X_gtp" = X.gtp, "Xconfig_zcp" = Xconfig.zcp, "test_fit" = FALSE, "run_model" = FALSE)

# Customize Map
Map.init<- fit.init$tmb_list$Map
Map.seas<- fit.intra$tmb_list$Map

# I'm a bit unclear what is happening here. The logsigmaXi piece is related to the variability in spatially varying GMRF, with the 1 and the 2 index corresponding to the linear predictor. If we had Xconfig be equal to 2 instead of 3 for season, would this stll be relevant or would that set log_sigma to unit variance?
if(Design %in% c("Both")){
  Map<- fit.intra$tmb_list$Map
  Map.logsigma.orig<- Map$log_sigmaXi1_cp
  str(Map.logsigma.orig) # Vector of length 22, [1] = 1, [2] = 2 and rest = NA
  Map$log_sigmaXi1_cp = factor(c(rep(Map$log_sigmaXi1_cp[1], dim(Xseason.itp)[3]), rep(Map$log_sigmaXi1_cp[dim(Xseason.itp)[3]+1], dim(Xyear.itp)[3])))
  Map$log_sigmaXi2_cp = factor(c(rep(Map$log_sigmaXi2_cp[1], dim(Xseason.itp)[3]), rep(Map$log_sigmaXi2_cp[dim(Xseason.itp)[3]+1], dim(Xyear.itp)[3])) )
  
  str(Map$log_sigmaXi1_cp) # Vector of length 22, [1] = 1, [2] = 1 and rest  = NA.
  
  # Again, still a bit unclear here what is going on here. Is this manually setting the variance to 1?
}else{
  stop("Not implemented")
}

# Refit with new mapping argument
# Check
fit.intra<- fit_model("settings" = settings.use, "Map" = Map, strata.limits = strat.limits, "Lat_i" = samp.dat[,'Lat'], "Lon_i" = samp.dat[,'Lon'], "t_i" = as.numeric(samp.dat[,"Season_Year"])-1, "b_i" = samp.dat[,'Catch_KG'], "c_iz" = rep(0, nrow(samp.dat)), "v_i" = rep(0, nrow(samp.dat)), "Q_ik" = NULL, "a_i" = rep(area.swept, nrow(samp.dat)), "PredTF_i" = samp.dat[,'PRED_TF'], covariate_data = cov.dat, formula = formula.use, "observations_LL" = cbind("Lat" = samp.dat[,'Lat'], "Lon" = samp.dat[, 'Lon']), "maximum_distance_from_sample" = maximum.dist.from.sample, "grid_dim_km" = grid.dim.km, "working_dir" = OutFile, "Use_REML" = Use_REML, "X_itp" = X.itp, "X_gtp" = X.gtp, "Xconfig_zcp" = Xconfig.zcp, "test_fit" = FALSE, "run_model" = FALSE)

# Run
fit.intra<- fit_model("settings" = settings.use, "Map" = Map, strata.limits = strat.limits, "Lat_i" = samp.dat[,'Lat'], "Lon_i" = samp.dat[,'Lon'], "t_i" = as.numeric(samp.dat[,"Season_Year"])-1, "b_i" = samp.dat[,'Catch_KG'], "c_iz" = rep(0, nrow(samp.dat)), "v_i" = rep(0, nrow(samp.dat)), "Q_ik" = NULL, "a_i" = rep(area.swept, nrow(samp.dat)), "PredTF_i" = samp.dat[,'PRED_TF'], covariate_data = cov.dat, formula = formula.use, "observations_LL" = cbind("Lat" = samp.dat[,'Lat'], "Lon" = samp.dat[, 'Lon']), "maximum_distance_from_sample" = maximum.dist.from.sample, "grid_dim_km" = grid.dim.km, "working_dir" = OutFile, "Use_REML" = Use_REML, "X_itp" = X.itp, "X_gtp" = X.gtp, "Xconfig_zcp" = Xconfig.zcp, "test_fit" = FALSE, "run_model" = TRUE)
  

#####
## NSF C-Accel Forecasting Challenge (12/16/2020)
#####
# Preliminaries ----------------------------------------------------------
# True/False switch as there might be a few packages that need to be installed IF running on DO server
# Libraries and sourcing functions
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
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

docker<- TRUE
if(docker){
  #source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VASTfunctions_AJAedits.R")
  #source("/home/andrew.allyn@gmail.com/GitHub/ForecastingChallenge/scripts/VAST_wrapper_func.R")
  source(here::here("scripts", "VASTfunctions_AJAedits.R"))
  source(here::here("scripts", "VAST_wrapper_func.R"))
  source(here::here("scripts", "fit_model_eff.R"))
} else {
  source(here::here("scripts", "VASTfunctions_AJAedits.R"))
  source(here::here("scripts", "VAST_wrapper_func.R"))
}

# Load in data -----------------------------------------------------------
if(docker){
  # Fisheries data
  dat<- read.csv(here::here("data", "ForecastingChallengeModelData.csv"))
  # Land shapefile for mapping
  #land<- st_read("/home/andrew.allyn@gmail.com/Shapefiles/ne_50m_land/ne_50m_land.shp")
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
max_dist_from_sample<- 50
grid_dim_km<- c(50, 50)
region_code<- "northwest_atlantic"
strat_limits<- data.frame("STRATA" = unique(dat$STRATUM)[order(unique(dat$STRATUM))])

# User supplied
# Make extrapolation grid from shapefile
if(docker){
  # NELME grid
  nelme_grid<- convert_shapefile(here::here("data", "NELME_sf.shp"), projargs = "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", projargs_for_shapefile = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", grid_dim_km = grid_dim_km, make_plots = FALSE, area_tolerance = 2)
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
# For running this in the "keep best" situation, going to want to start somewhere close to the ideal model, rather than building up from the pure environment-only SDM to hopefully reduce the number of models that need to be fit and increase speed and efficiency. With an eye towards using these models to forecast species distributions/abundance, the ideal forecasting model would include both spatial and spatio-temporal correlation to account for unmeasured environmental/biological processes that are either persistent throughout the entire time series or ephemeral. In addition, we are going to want to have some temporal autoregressive structure for both the intercepts and the spatio-temporal correlation as this will facilitate predicting these model components in the future.  
fieldconfig_forebase<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_forebase<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 1, "Epsilon2" = 1)

# Observation model ("Poisson" link function)
obsmodel_use<- c("PosDist" = 4, "EncProbForm" = 1)

# Options -- want SD_observation_density
options_use<- c("SD_site_logdensity" = FALSE, "Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "Calculate_Cov_SE" = FALSE, "Calculate_Synchrony" = FALSE, "Calculate_proportion" = FALSE, "SD_observation_density" = FALSE)

settings_forebase<- make_settings(n_x = n_x_use, Region = "other", strata.limits = strat_limits, purpose = "index2", FieldConfig = fieldconfig_forebase, RhoConfig = rhoconfig_forebase, ObsModel = obsmodel_use, Options = options_use, use_anisotropy = TRUE, bias.correct = FALSE)

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
  out_folder<- shared.path(os = os.use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results/")
}

dat_nest<- dat_nest %>% 
  filter(., COMNAME == "ATLANTIC COD")

# Add in the forechallenges piece...
fore_challenges_df<- data.frame("COMNAME" = rep(unique(dat_nest$COMNAME), each = length(fore_challenges)), "fore_challenge" = rep(fore_challenges, length(unique(dat_nest$COMNAME))))

dat_nest<- dat_nest %>%
  left_join(., fore_challenges_df, by = c("COMNAME" = "COMNAME"))

# Detecting cores
cores_avail<- detectCores()
library(doFuture)
registerDoFuture()
plan(multisession, workers = cores_avail-2)

# I kept getting errors with the compilation here when trying to run in parallel. Not sure why, but trying something by running a simple example to get the VAST_cpp and VAST_so files. Then copying those into the folders as needed.
## Setup really simple example for testing
library(VAST)
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

vast_files<- c(paste(here::here(), "VAST_v12_0_0.cpp", sep = "/"), paste(here::here(), "VAST_v12_0_0.so", sep = "/"), paste(here::here(), "VAST_v12_0_0.o", sep = "/"))

# Seemed really slow...
#all<- dat_nest
all<- dat_nest[1:2,]

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
  
  #file.copy(vast_files, outfolder)
  
  # Create sub-folder
  # outfile<- paste(outfolder, "/", paste(fore_dates, season_run, spp_run, sep = "_"), sep = "")
  outfile<- paste(outfolder, "/", fore_dates, sep = "")
  # if(!file.exists(outfile)){
  #   dir.create(outfile)
  # }
  
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
  
  # Base, set all years to be the same -- this could be needed if we struggle to fit any of the more complex models
  cov_dat_base<- cov_dat_all
  cov_dat_base$Year<- rep(min(cov_dat_base$Year), nrow(cov_dat_base))
  
  # Model fit wrapper function
  fit_base<- fit_model_eff("settings" = settings_forebase,
                               # Spatial info
                               observations_LL = cbind("Lat" = samp_dat_all[,'Lat'], "Lon" = samp_dat_all[, 'Lon']), grid_dim_km = grid_dim_km, make_plots = TRUE, 
                               # Model info
                               "Lat_i" = samp_dat_all[,'Lat'], "Lon_i" = samp_dat_all[,'Lon'], "t_i" = as.vector(samp_dat_all[,'Year']), "b_i" = samp_dat_base[,'Catch_KG'], "c_iz" = rep(0, nrow(samp_dat_all)), "v_i" = rep(0, nrow(samp_dat_all)), "Q_ik" = NULL, "a_i" = rep(area_swept, nrow(samp_dat_all)), covariate_data = cov_dat_all, X1_formula = formula_use, X2_formula = formula_use, "PredTF_i" = samp_dat_all[,'PRED_TF'], "working_dir" = outfolder, "CompileDir" = here::here(), "run_model" = TRUE, "test_fit" = FALSE, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, Use_REML = TRUE)
  
  # If our basic forecasting model worked, let's add in some complexity...
  if(class(fit_base) == "fit_model" && (max(abs(fit_base$parameter_estimates$diagnostics$final_gradient)) <= 0.01)){
    # Append progress file
    progress_new<- "Base model has passed fit checks, trying to add complexity"
    write(progress_new, file = paste(outfile, "_progress.txt", sep = ""), append = TRUE)
    saveRDS(fit_base, file = paste(outfile, "_base_modelfit.rds", sep = "")) 
  }
}

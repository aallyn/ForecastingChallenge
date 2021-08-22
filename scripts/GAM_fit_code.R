#####
## Single species delta GAM model
#####

# Preliminaries -----------------------------------------------------------
# To run this code, we are going to need a few different functions in existing libraries. If for some reason these lines throw a "library does not exist" error, then you will need to install them and then load them.
library(mgcv)
library(tidyverse)
library(snakecase)
library(here)
library(conflicted)
library(gmRi)
library(sf)
source(here::here("scripts", "GAM_fit_functions.R"))
source(here::here("scripts", "VASTfunctions_AJAedits.R"))
filter <- dplyr::filter
select <- dplyr::select

# Load in data -----------------------------------------------------------
# Fisheries data
dat<- read.csv(here::here("data", "ForecastingChallengeModelData.csv"))
# Land shapefile for mapping
land<- st_read(here::here("data", "NELME_sf.shp"))

# Some processing
dat<- dat %>%
  mutate(., "PRESENCE" = ifelse(SUM_BIOMASS > 0, 1, 0),
         "LOGBIO" = log(SUM_BIOMASS))

# Create a nested data frame
dat_nest<- dat %>%
  group_by(., SEASON, COMNAME) %>%
  nest()

# Add in the forechallenges piece...
fore_challenges<- seq(from = 2004, to = 2019, by =  1)
fore_challenges_df<- data.frame("COMNAME" = rep(unique(dat_nest$COMNAME), each = length(fore_challenges)), "fore_challenge" = rep(fore_challenges, length(unique(dat_nest$COMNAME))))

dat_nest<- dat_nest %>%
  left_join(., fore_challenges_df, by = c("COMNAME" = "COMNAME"))

# Model fitting ----------------------------------------------------------
# Set path to main root for model output
out_folder<- "~/Box/Mills Lab/Projects/ForecastingChallenge/Temp Results/"

# Map the model fitting function to the nested data frame
dat_nest<- dat_nest %>%
  mutate(., "GAM_Fits_Naive" = pmap(list(season = SEASON, comname = COMNAME, data = data, fore_challenge = fore_challenge, covs = list(c("AVGDEPTH", "SODA_BT", "SODA_SST")), out_folder = out_folder), possibly(deltalog_mod_fit_naive, NA)),
         "GAM_Fits_Fore_M1" = pmap(list(season = SEASON, comname = COMNAME, data = data, fore_challenge = fore_challenge, covs = list(c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")), covs_scale = list(covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")), year_m = list(1), out_folder = out_folder), possibly(deltalog_mod_fit_fore, NA)),
         "GAM_Fits_Fore_M2" = pmap(list(season = SEASON, comname = COMNAME, data = data, fore_challenge = fore_challenge, covs = list(c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")), covs_scale = list(covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")), year_m = list(2), out_folder = out_folder), possibly(deltalog_mod_fit_fore, NA)))
    
# Model predictions  -------------------------------------------------------
dat_nest<- dat_nest %>%
  mutate(., "GAM_Preds_Naive" = pmap(list(season = SEASON, comname = COMNAME, data = data, mod_list = GAM_Fits_Naive, fore_challenge = fore_challenge, covs = list(c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")), covs_scale = list(covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")), mod_form = "Naive", out_folder = out_folder), possibly(deltalog_mod_predict, NA)),
         "GAM_Preds_Fore_M1" = pmap(list(season = SEASON, comname = COMNAME, data = data, mod_list = GAM_Fits_Fore_M1, fore_challenge = fore_challenge, covs = list(c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")), covs_scale = list(covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")), mod_form = "Fore_YearM_1", out_folder = out_folder), possibly(deltalog_mod_predict, NA)),
         "GAM_Preds_Fore_M2" = pmap(list(season = SEASON, comname = COMNAME, data = data, mod_list = GAM_Fits_Fore_M2, fore_challenge = fore_challenge, covs = list(c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")), covs_scale = list(covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")), mod_form = "Fore_YearM_2", out_folder = out_folder), possibly(deltalog_mod_predict, NA)))















# One residual plot -- biomass
pred.df<- pred.df %>%
  mutate(., "Sample" = sample(1:nrow(.)))

# Randomly shuffle the data to remove time to match preliminary plot from Cornell
pred.df.plot<- pred.df[pred.df$Sample, ]
pred.df.plot$xlab<- seq(from = 1, to = nrow(pred.df.plot))
pred.vs.obs<- ggplot() +
  geom_point(data = pred.df.plot, aes(x = xlab, y = wtcpue), shape = 1, color = "blue") +
  geom_point(data = pred.df.plot, aes(x = xlab, y = predicted.bio), shape = 4, color = "green") +
  ylab("Weight Catch Per Unit Effort (kg)") +
  xlab("Sample") +
  theme_bw()
ggsave(here::here("temp results", paste(spp.mod, "_GAM_PredObs_", dat.run, ".jpg", sep = "")), pred.vs.obs)

# Resid summary
resids<- ggplot() +
  geom_point(data = pred.df.plot, aes(x = xlab, y = bio.resid, color = bio.resid)) +
  scale_color_gradient2(name = "Residual", midpoint = 0) +
  ylab("Residaul weight Catch Per Unit Effort (kg)") +
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(here::here("temp results", paste(spp.mod, "_GAM_Resids_", dat.run, ".jpg", sep = "")), resids)

# Save
write.csv(pred.df, here::here("temp results", paste(spp.mod, "Predictions_", dat.run, "Data.csv", sep = "")))

# Things to do moving forward:
# More formal prediction assessment (RMSE?)
# Spatial/spatio-temporal patterns in residuals
#####
## Single species delta GAM model
#####

# Preliminaries -----------------------------------------------------------
# To run this code, we are going to need a few different functions in existing libraries. If for some reason these lines throw a "library does not exist" error, then you will need to install them and then load them.
library(mgcv)
library(tidyverse)
library(snakecase)
library(here)
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")
source(here::here("scripts", "DeltaLogNormal_Functions.R"))

# Load data, model settings -----------------------------------------------
dat.all<- read_csv(here::here("data", "AllDataForModels_03192020.csv"))

# OLD!!! Data depending on if we are using Michelle's (groundfish_trianing/groundfish_testing) or Shufeng's (Morely_biomass/Morley_environmental). To switch, change the dat.run to dat.run<- "Michelle" or dat.run<- "Shufeng"
if(FALSE){
  dat.run<- "Michelle"
  if(dat.run == "Michelle"){
    # Load in the two RData files
    dat.train<- readRDS(here::here("data", "groundfish_training.rds"))
    dat.test<- readRDS(here("data", "groundfish_testing.rds"))
    
    # Combining them to one data set so we have more flexibility defining when the training/testing time period cut is.
    dat.all<- bind_rows(dat.train, dat.test)
    
    # Table with common names
    name.table<- read_csv(here::here("data", "sppocean_and_commonnames_table.csv"))
    
    # Add common names to the data
    dat.all<- dat.all %>%
      left_join(., name.table)
    
    # Presence/Absence column
    dat.all<- dat.all %>%
      mutate(., "presenceabsence" = ifelse(wtcpue > 0, 1, 0))
    
  } else if(dat.run == "Shufeng"){
    
    # Load in the two csv files
    dat.features<- read.csv(here::here("data", "features_processed.csv"))
    dat.labels<- read.csv(here::here("data", "labels_processed.csv"))
    
    # Going to want one dataframe...
    # First, let's gather the dat.labels so that we have haulid, year, species, rather than each species as it's own column. Also, add logwtcpue.
    dat.labels.l<- dat.labels %>%
      gather(., commonname, wtcpue, -haulid, -year) %>%
      mutate(., "logwtcpue" = log(wtcpue))
    dat.labels.l$logwtcpue[is.infinite(dat.labels.l$logwtcpue)]<- 0 # This seems....weird. But, log(0) in other dataset is 0, so  doing the same here.
    
    # Arrange by year and haulid
    dat.labels.l<- dat.labels.l %>%
      arrange(., year, haulid, commonname)
    
    # Now, join the environmental data to the species data
    dat.all<- dat.labels.l %>%
      left_join(., dat.features, by = c("haulid"))
    
    # Adjust common name
    dat.all$commonname<- to_any_case(dat.all$commonname, case = "sentence")
    
    # No surveyfact variable here...need to add it from Michelle's data
    dat.supp<- bind_rows(readRDS(here::here("data", "groundfish_training.rds")), readRDS(here::here("data", "groundfish_testing.rds"))) %>%
      dplyr::select(., haulid, lat, lon, year, month, regionfact, surveyfact) %>%
      distinct()
    
    dat.all<- dat.all %>%
      left_join(., dat.supp) %>%
      drop_na(surveyfact)
    
    # Some naming differences of variables/features. Adjusting those so that they match with what we have below for the model formula
    names(dat.all)[c(18, 19, 20, 21, 22)]<- c("SST.min", "SST.max", "SBT.seasonal", "SBT.min", "SBT.max")
    
    # Presence/absence column
    dat.all<- dat.all %>%
      mutate(., "presenceabsence" = ifelse(wtcpue > 0, 1, 0))
  }
}

## Model settings
# Setting up the forecasting experiment...
year.cuts<- c(1999, 2011, 2012, 2013)

# Only going to use the NEFSC Spring and Fall survey
surveyfact.keep<- c("NEFSC_NEUSFall", "NEFSC_NEUSSpring")

# For each of these, going to fit models to each species, save the model, make predictions and save the predictions. To make that a process a bit easier, going to nest the data and then map model fitting/prediction functions to the nested dataframe.
dat.nest<- dat.all %>%
  filter(., surveyfact %in% surveyfact.keep) %>%
  group_by(commonname, surveyfact) %>%
  nest()

# Model fitting ----------------------------------------------------------
# Map the model fitting function to the nested dataframe, once for each time period
dat.nest<- dat.nest %>%
  mutate(., "Mod.1999" = map2(data, year.cuts[1], deltalog_mod_fit),
         "Mod.2011" = map2(data, year.cuts[2], deltalog_mod_fit),
         "Mod.2012" = map2(data, year.cuts[3], deltalog_mod_fit),
         "Mod.2013" = map2(data, year.cuts[4], deltalog_mod_fit))

# Model predictions  -------------------------------------------------------
# Map prediction function to the nested dataframe, once for each time period
dat.nest<- dat.nest %>%
  mutate(., "Preds.2000to2014" = pmap(list(data = data, mod.fit = Mod.1999, split.year = year.cuts[1]), deltalog_mod_predict),
         "Preds.2012to2014" = pmap(list(data = data, mod.fit = Mod.2011, split.year = year.cuts[2]), deltalog_mod_predict),
         "Preds.2013to2014" = pmap(list(data = data, mod.fit = Mod.2012, split.year = year.cuts[3]), deltalog_mod_predict),
         "Preds.2014" = pmap(list(data = data, mod.fit = Mod.2013, split.year = year.cuts[4]), deltalog_mod_predict))
save(dat.nest, file = "~/Box/Mills Lab/Projects/SDM-convergence/temp results/DeltaLogNormalGAMs.rds")

# Residuals
pred.df<- pred.df %>%
  mutate(., "presenceabsence.resid" = presenceabsence - predicted.prob.presence,
         "logbio.resid" = logwtcpue - predicted.logbio,
         "bio.resid" = wtcpue - predicted.bio)

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
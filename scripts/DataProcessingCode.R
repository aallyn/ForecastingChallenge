#####
## Data prep for NSF C-Accel forecasting data
#####
# Preliminaries ----------------------------------------------------------
# Libraries and sourcing functions
library(tidyverse)
library(here)
library(sf)
library(raster)
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")

# Load raw trawl data ----------------------------------------------------------
load(paste(res.data.path, "NMFS_trawl/Survdat_Nye_Aug 2020.RData", sep = ""))

# Create a "null" dataset FIRST as this will be used to impute absences and we want to make sure it includes all tows within area of interest and seasons, regardless of what was actually caught there.
# Start and end years
year.start<- 1980
year.end<- 2017

null.dat.temp<- survdat %>%
  filter(STRATUM >= 01010 & STRATUM <= 01760) %>%
  filter(STRATUM!=1310 & STRATUM!=1320 & STRATUM!=1330 & STRATUM!=1350 &
           STRATUM!=1410 & STRATUM!=1420 & STRATUM!=1490) %>%
  filter(SEASON == "SPRING" | SEASON == "FALL") %>%
  filter(EST_YEAR >= year.start & EST_YEAR <= year.end) 

# Expand tow/species grid
null.dat<- expand.grid("ID" = unique(null.dat.temp$ID), "COMNAME" = unique(null.dat.temp$COMNAME)) %>%
  mutate(., "COMNAME" = as.character(COMNAME))

# Processing -- following LG's processing code:
  # Select columns of interest 
  # Select stratum for full offshore survey
  # Eliminate strata not sampled for full time span
  # Select fall and spring surveys only
  # Select years based on year.start and year.end
  # If biomass is NA but abundance is >0, make biomass be 0.01
  # If biomass is 0 and abundance is >0, make biomass be 0.01
  # If abundance is NA and biomass >0, make abundance 1
  # If abundance is 0 and biomass >0, make abundance 1
  # Remove values where biomass AND abundance are NAs
  # Keep unique tow-species-sex rows
  # Calculate sum of biomass and abundance per tow-species

# Create trawl data set
dat<- survdat %>%
  dplyr::select(ID, EST_YEAR, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, CATCHSEX, BIOMASS, AVGDEPTH, ABUNDANCE, LENGTH, NUMLEN) %>%
  filter(STRATUM >= 01010 & STRATUM <= 01760) %>%
  filter(STRATUM!=1310 & STRATUM!=1320 & STRATUM!=1330 & STRATUM!=1350 &
           STRATUM!=1410 & STRATUM!=1420 & STRATUM!=1490) %>%
  filter(SEASON == "SPRING" | SEASON == "FALL") %>%
  filter(EST_YEAR >= year.start & EST_YEAR <= year.end) %>%
  mutate(BIOMASS = ifelse(is.na(BIOMASS) == TRUE & ABUNDANCE > 0,0.01, BIOMASS)) %>% 
  mutate(BIOMASS = ifelse(BIOMASS == 0 & ABUNDANCE > 0,0.01, BIOMASS)) %>% 
  mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == TRUE & BIOMASS > 0, 1, ABUNDANCE)) %>% 
  mutate(ABUNDANCE = ifelse(ABUNDANCE == 0 & BIOMASS > 0, 1, ABUNDANCE)) %>%   
  filter(!is.na(BIOMASS),
         !is.na(ABUNDANCE)) %>%
  distinct(ID, EST_YEAR, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, CATCHSEX, .keep_all = TRUE) %>%
  group_by(., ID, EST_YEAR, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, AVGDEPTH) %>%
  summarize(., "SUM_BIOMASS" = sum(BIOMASS, na.rm = TRUE),
            "SUM_ABUNDANCE" = sum(ABUNDANCE, na.rm = TRUE))

glimpse(dat)

# Now, create groundfish dataset...
groundfish<- c("AMERICAN PLAICE", "ATLANTIC COD", "ATLANTIC HALIBUT", "POLLOCK", "ATLANTIC WOLFFISH", "HADDOCK", "OCEAN POUT", "ACADIAN REDFISH", "WHITE HAKE", "WINDOWPANE", "WINTER FLOUNDER", "WITCH FLOUNDER", "YELLOWTAIL FLOUNDER")

dat.ground<- dat %>%
  ungroup() %>%
  mutate(., "COMNAME" = as.character(COMNAME)) %>%
  filter(., COMNAME %in% groundfish) 

null.dat.ground<- null.dat %>%
  filter(., COMNAME %in% groundfish)

# Add in "absences" from null.dat. Just focus on ID, COMNAME, SUM_BIOMASS, SUM_ABUNDANCE.  
occu.ground<- dat.ground %>%
  dplyr::select(., ID, COMNAME, SUM_BIOMASS, SUM_ABUNDANCE)

# Impute missing ID/COMNAME combos
occu.ground<- occu.ground %>%
  full_join(., null.dat.ground) %>%
  arrange(., ID, COMNAME)

# Update SUM_BIOMASS, SUM_ABUNDANCE to zero if NA 
occu.ground<- occu.ground %>%
  mutate_at(., c("SUM_BIOMASS", "SUM_ABUNDANCE"), ~replace_na(., 0))
  
# Add back in trawl info...
env.ground<- dat.ground %>%
  select(., c(ID, EST_YEAR, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, AVGDEPTH))

# Final model data set
dat.mod<- env.ground %>%
  left_join(., occu.ground, by = c("ID" = "ID"))

# Check -- tows per species species per tow
dat.mod %>%
  group_by(., ID) %>%
  summarize(., "SppPerTow" = n_distinct(COMNAME))

dat.mod %>%
  group_by(., COMNAME) %>%
  summarize(., "TowPerSpp" = n_distinct(ID))

# Extracting environmental covariates at trawl locations ----------------------------------------------------------
## SODA environmental data
soda<- raster::stack("/home/andrew.allyn@gmail.com/ForecastingChallenge/Data/soda3.4.2_mn_ocean_reg_1980.nc")

# Test save
writeRaster(soda, filename = "/home/rstudio/ForecastingChallenge/Temp Results/soda_temp.grd")
# Save it....
saveRDS(dat.mod, file = here::here("data", "ForecastingChallengeModelData.rds"))
write.csv(dat.mod, file = here::here("data", "ForecastingChallengeModelData.csv"))

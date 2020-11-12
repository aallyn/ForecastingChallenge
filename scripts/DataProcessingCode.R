#####
## Data prep for NSF C-Accel forecasting data
#####
# Preliminaries ----------------------------------------------------------
# Libraries and sourcing functions
library(tidyverse)
library(here)
library(sf)
library(raster)
library(ncdf4)
library(patchwork)
library(conflicted)
library(gmRi)
filter<- dplyr::filter
select<- dplyr::select
extract<- raster::extract
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")

# Land shapefile for visualization
res_dat_path<- shared.path(os.use = os.use, group = "RES Data")
land<- st_read(paste(res_dat_path, "Shapefiles/ne_50m_land/ne_50m_land.shp", sep = "")) 
xlim_use<- c(-76, -65) 
ylim_use<- c(35, 45)

# Depth bathymetry
depth<- raster(paste(res_dat_path, "Shapefiles/NEShelf_Etopo1_bathy.tiff", sep = "")) 

# Load raw trawl data ----------------------------------------------------------
load(paste(res_dat_path, "NMFS_trawl/Survdat_Nye_Aug 2020.RData", sep = ""))

# Create a "null" dataset FIRST as this will be used to impute absences and we want to make sure it includes all tows within area of interest and seasons, regardless of what was actually caught there.
# Start and end years
year_start<- 1980
year_end<- 2019

null_dat_temp<- survdat %>%
  filter(STRATUM >= 01010 & STRATUM <= 01760) %>%
  filter(STRATUM!=1310 & STRATUM!=1320 & STRATUM!=1330 & STRATUM!=1350 &
           STRATUM!=1410 & STRATUM!=1420 & STRATUM!=1490) %>%
  filter(SEASON == "SPRING" | SEASON == "FALL") %>%
  filter(EST_YEAR >= year_start & EST_YEAR <= year_end) 

# Expand tow/species grid
null_dat<- expand.grid("ID" = unique(null_dat_temp$ID), "COMNAME" = unique(null_dat_temp$COMNAME)) %>%
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
  dplyr::select(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, CATCHSEX, BIOMASS, AVGDEPTH, ABUNDANCE, LENGTH, NUMLEN) %>%
  filter(STRATUM >= 01010 & STRATUM <= 01760) %>%
  filter(STRATUM!=1310 & STRATUM!=1320 & STRATUM!=1330 & STRATUM!=1350 &
           STRATUM!=1410 & STRATUM!=1420 & STRATUM!=1490) %>%
  filter(SEASON == "SPRING" | SEASON == "FALL") %>%
  filter(EST_YEAR >= year_start & EST_YEAR <= year_end) %>%
  mutate(BIOMASS = ifelse(is.na(BIOMASS) == TRUE & ABUNDANCE > 0,0.01, BIOMASS)) %>% 
  mutate(BIOMASS = ifelse(BIOMASS == 0 & ABUNDANCE > 0,0.01, BIOMASS)) %>% 
  mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == TRUE & BIOMASS > 0, 1, ABUNDANCE)) %>% 
  mutate(ABUNDANCE = ifelse(ABUNDANCE == 0 & BIOMASS > 0, 1, ABUNDANCE)) %>%   
  filter(!is.na(BIOMASS),
         !is.na(ABUNDANCE)) %>%
  distinct(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, CATCHSEX, .keep_all = TRUE) %>%
  group_by(., ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, SVSPP, COMNAME, AVGDEPTH) %>%
  summarize(., "SUM_BIOMASS" = sum(BIOMASS, na.rm = TRUE),
            "SUM_ABUNDANCE" = sum(ABUNDANCE, na.rm = TRUE))

glimpse(dat)

# Now, create groundfish dataset...
groundfish<- c("AMERICAN PLAICE", "ATLANTIC COD", "ATLANTIC HALIBUT", "POLLOCK", "ATLANTIC WOLFFISH", "HADDOCK", "OCEAN POUT", "ACADIAN REDFISH", "WHITE HAKE", "WINDOWPANE", "WINTER FLOUNDER", "WITCH FLOUNDER", "YELLOWTAIL FLOUNDER")

dat_ground<- dat %>%
  ungroup() %>%
  mutate(., "COMNAME" = as.character(COMNAME)) %>%
  filter(., COMNAME %in% groundfish) 

null_dat_ground<- null_dat %>%
  filter(., COMNAME %in% groundfish)

# Add in "absences" from null.dat. Just focus on ID, COMNAME, SUM_BIOMASS, SUM_ABUNDANCE.  
occu_ground<- dat_ground %>%
  dplyr::select(., ID, COMNAME, SUM_BIOMASS, SUM_ABUNDANCE)

# Impute missing ID/COMNAME combos
occu_ground<- occu_ground %>%
  full_join(., null_dat_ground) %>%
  arrange(., ID, COMNAME)

# Update SUM_BIOMASS, SUM_ABUNDANCE to zero if NA 
occu_ground<- occu_ground %>%
  mutate_at(., c("SUM_BIOMASS", "SUM_ABUNDANCE"), ~replace_na(., 0))

# Add back in trawl info...
env_ground<- dat_ground %>%
  dplyr::select(., c(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, AVGDEPTH)) %>%
  distinct()

# Final model data set
dat_mod<- env_ground %>%
  left_join(., occu_ground, by = c("ID" = "ID")) 
dat_mod$ROW_ID<- seq(from = 1, to = nrow(dat_mod), by = 1)

# Check -- tows per species species per tow
dat_mod %>%
  group_by(., ID) %>%
  summarize(., "SppPerTow" = n_distinct(COMNAME))

dat_mod %>%
  group_by(., COMNAME) %>%
  summarize(., "TowPerSpp" = n_distinct(ID))

# Extracting environmental covariates at trawl locations ----------------------------------------------------------
## SODA environmental data
# Have we done this already?
soda_done<- file.exists("~/SODA/sst.grd")
if(!soda_done){
  soda_files<- list.files("~/SODA/", full.names = TRUE)
  soda_files<- soda_files[grepl("[[:digit:]]", soda_files)]
  
  # Combine them
  soda_sst_out<- raster::stack()
  
  for(i in seq_along(soda_files)){
    # SST 
    sst_rast_temp<- raster::rotate(raster::brick(soda_files[[i]], varname = "temp", "level" = 1))
    if(i == 1){
      soda_sst_out<- raster::stack(sst_rast_temp, soda_sst_out)
    } else {
      soda_sst_out<- raster::stack(soda_sst_out, sst_rast_temp)
    }
  }
  # Save em
  writeRaster(soda_sst_out, filename = "~/SODA/sst.grd")
} else {
  soda_files<- list.files("~/SODA/", full.names = TRUE)
  soda_files<- soda_files[grepl("[[:digit:]]", soda_files)]
  soda_sst_out<- raster::stack("~/SODA/sst.grd")
}

# Bottom temperature is going to be more complicated, mostly because the model depth (z level) is going to be a discrete bin and might not match the "true" depth. I went about this one way before where I basically got the "non-NA" temperature value from the model depth that was closest to the true depth. This worked well on occasion, but also ended up with some very problematic results (for example, Wilkinson Basin where the true depth is ~200 m and there isn't a temperature value from the model until depth level 5 (45 meters)). So...for any given time step, want to find the "bottom" depth at each of the locations. To make this a bit easier, we can try to reduce this a bit as we don't need to go through ALL of the potential depth levels.
soda_temp<- nc_open("~/SODA/soda3.4.2_mn_ocean_reg_1980.nc")
soda_depths<- ncvar_get(soda_temp, var = "st_ocean")

closest_depth_func<- function(x, vec = soda_depths){
  temp_depths<- soda_depths[which(soda_depths <= x)]
  if(length(temp_depths) == 0){
    closest_depth_out<- min(soda_depths)
  } else {
    temp_depths2<- temp_depths[!temp_depths == x]
    closest_depth_out<- temp_depths2[which.min(abs(temp_depths2 - x))]
  }
  return(closest_depth_out)
}

dat_mod<- dat_mod %>%
  mutate(., "SODA_DEPTH" = map_dbl(AVGDEPTH, closest_depth_func))
dat_mod$SODA_LEVEL<- match(dat_mod$SODA_DEPTH, soda_depths)

soda_get_levels<- sort(unique(dat_mod$SODA_LEVEL))
soda_out<- vector("list", length(soda_get_levels))

rast_bbox<- extent(c(-80, -65, 34, 46))

# Mask depth raster, which is the resolution we will interpolate the bottom temperatures too...
bt_interp<- crop(depth, rast_bbox)
bt_interp[]<- NA

for(i in seq_along(soda_files)){
  
  for(j in seq_along(soda_get_levels)){
    # Get temp  
    lt_rast_temp<- raster::rotate(raster::brick(soda_files[[i]], varname = "temp", "level" = soda_get_levels[j]))
    lt_rast_crop<- crop(lt_rast_temp, rast_bbox)
    
    # Store it in the right list area
    if(i == 1){
      soda_out[[j]]<- lt_rast_crop
    } else {
      soda_out[[j]]<- raster::stack(soda_out[[j]], lt_rast_crop)
    }  
  }
  print(paste(soda_files[[i]], " is done!", sep = ""))
}

# Okay...so now we have a list. There are 23 elements in the list (1:level 23 of the model) and then within each level of the model (list element) there are...480 rasters, which correspond to year-months. Ultimately, want a stack by DAY, where the layers are depths (down the list)...
bt_stack_out<- raster::stack()

for(i in 1:dim(soda_out[[1]])[3]){
  level_use<- i
  yrmonth_stack<- raster::stack(sapply(soda_out, "[[", level_use))
  
  # Convert to dataframe
  yrmonth_df<- as.data.frame(yrmonth_stack, xy = TRUE)
  
  # Add a row ID to make things easier
  yrmonth_df$row_id<- seq(from = 1, to = nrow(yrmonth_df))
  
  # Alright, some rows are going to have all NAs...no sense dealing with those.
  vars_ignore<- c("x", "y", "row_id")
  filter_vars<- colnames(yrmonth_df)[!colnames(yrmonth_df) %in% vars_ignore]
  yrmonth_good<- yrmonth_df %>%
    filter_at(all_of(filter_vars), any_vars(!is.na(.))) %>%
    select(., all_of(c(vars_ignore, filter_vars)))
  
  # That's good...now, for each location or row, we want to keep the last value BEFORE an NA...
  yrmonth_good$BT_keep<- apply(yrmonth_good, 1, function(x) x[max(which(!is.na(x)))])
  
  # New dataframe...
  keep_first<- colnames(yrmonth_df)[3]
  bt_out<- yrmonth_df %>%
    select(., all_of(c(vars_ignore, keep_first)))
  
  # Add in "bottom" temps...
  bt_out[,4][match(yrmonth_good$row_id, bt_out$row_id)]<- yrmonth_good$BT_keep
  names(bt_out)[4]<- gsub(".[^.]+$", "", keep_first)
  
  # Convert back to a raster...
  bt_rast<- rasterFromXYZ(bt_out[,c(1,2,4)])
  bt_rast_out<- resample(bt_rast, bt_interp, method = "bilinear")
  
  if(i == 1){
    bt_stack_out<- stack(bt_rast_out, bt_stack_out)
  } else {
    bt_stack_out<- stack(bt_stack_out, bt_rast_out)
  }

}

# I think that worked. What about checking the difference between this and sst?
soda_sst_fine<- crop(soda_sst_out, rast_bbox)
soda_sst_fine<- resample(soda_sst_fine, bt_interp, method = "bilinear")

# Difference
sst_bt_diff<- soda_sst_fine - bt_stack_out

# Now, extract BT and SST. Only need the spatial information here and don't want to extract from duplicated tows because of multiple species records per tow.
dat_env<- dat_mod %>%
  distinct(EST_YEAR, EST_MONTH, DECDEG_BEGLON, DECDEG_BEGLAT) %>%
  st_as_sf(., coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
dat_extract_bt<- extract(bt_stack_out, dat_env)

# Get the right year-month
i_temp<- match(paste("X", dat_env$EST_YEAR, ".", str_pad(dat_env$EST_MONTH, 2, "left", pad = 0), sep = ""), gsub(".[^.]+$", "", names(bt_stack_out))) 
i_mat<- cbind(1:nrow(dat_env), i_temp) 

# Add it...
dat_env$SODA_BT<- rep(NA, nrow(dat_env))
dat_env$SODA_BT<- dat_extract_bt[i_mat]

# SST
dat_extract_sst<- extract(soda_sst_out, dat_env)

# Get the right year-month
i_temp<- match(paste("X", dat_env$EST_YEAR, ".", str_pad(dat_env$EST_MONTH, 2, "left", pad = 0), sep = ""), gsub(".[^.]+$", "", names(soda_sst_out))) 
i_mat <- cbind(1:nrow(dat_env), i_temp) 

# Bind result to dat env
dat_env$SODA_SST<- rep(NA, nrow(dat_env))
dat_env$SODA_SST<- dat_extract_sst[i_mat]

# Now, bring over BT and SST
dat_env_join<- dat_env %>%
  select(., EST_YEAR, EST_MONTH, DECDEG_BEGLON, DECDEG_BEGLAT, SODA_BT, SODA_SST) %>%
  st_drop_geometry()

dat_mod<- dat_mod %>%
  left_join(., dat_env_join, by = c("EST_YEAR", "EST_MONTH", "DECDEG_BEGLON", "DECDEG_BEGLAT"))

# Save it....
saveRDS(dat_mod, file = "~/Box/Mills Lab/Projects/ForecastingChallenge/Data/ForecastingChallengeModelData.rds")
write.csv(dat_mod, file = "~/Box/Mills Lab/Projects/ForecastingChallenge/Data/ForecastingChallengeModelData.csv")

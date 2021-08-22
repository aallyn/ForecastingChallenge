#####
## NSF C-Accel Forecasting Challenge (01/19/2021): Visualizing results
#####
# Preliminaries ----------------------------------------------------------
library(tidyverse)
library(gmRi)
library(PresenceAbsence)
library(MLmetrics)
library(forecast)
library(raster)
library(sf)
library(conflicted)
library(stringi)
library(rnaturalearth)
library(rnaturalearthdata)
select<- dplyr::select
filter<- dplyr::filter
os_use<- "unix"
source(here::here("scripts", "SDM_PredValidation_Functions.R"))
proj_box_path<- shared.path(os.use = os_use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/")

# Getting haul information ------------------------------------------------
dat_all<-  read.csv(here::here("data", "ForecastingChallengeModelData.csv"))

# Going to want some other information from the dat_all as these just have haulid
hauls<- dat_all %>%
  dplyr::select(., ROW_ID, EST_YEAR, EST_MONTH, DECDEG_BEGLON, DECDEG_BEGLAT) %>%
  distinct() 
names(hauls)<- c("haulid", "year", "month", "lon", "lat")

# Visualizing study area ------------------------------------------------
land<- st_read("~/Box/RES_Data/Shapefiles/NELME_Regions/NELME_sf.shp")

# Visualizing forecasting challenge options
# Need SST data 
oisst_path<- shared.path(os.use = os_use, group = "root", folder = "RES_Data/OISST/")
oisst_ts<- raster::stack(paste(oisst_path, "DailyAnomsThroughFeb2020.grd", sep = ""))

# Convert to data frame
oisst_df<- as.data.frame(oisst_ts, xy = TRUE) %>%
  gather(., "Date", "SST", -x, -y)

oisst_df2<- oisst_df %>%
  drop_na(SST) %>%
  mutate(., "Date" = as.Date(gsub("[.]", "-", gsub("X", "", Date))))

oisst_yrmonth<- oisst_df2 %>%
  mutate(., "Year_Month" = format(Date, "%Y-%m")) %>%
  group_by(., Year_Month) %>%
  summarize(., "Mean_SST" = mean(SST, na.rm = TRUE))

oisst_yr<- oisst_df2 %>%
  mutate(., "Year" = format(Date, "%Y")) %>%
  group_by(., Year) %>%
  summarize(., "Mean_SST" = mean(SST, na.rm = TRUE)) %>%
  filter(., Year < 2020) %>%
  mutate(., Year = as.numeric(Year),
         "Year_Model" = Year - min(Year, na.rm = TRUE))
# Plots
# Full trend line plot
lm_full<- lm(Mean_SST ~ Year_Model, data = oisst_yr)
my_formula_full<- paste("SST Anom = ", as.numeric(round(coef(lm_full)[1], 2)), " + ", as.numeric(round(coef(lm_full)[2], 2)), "*Year", sep = "")
adj_r2_full<- paste("Adj R2 = ", round(summary(lm_full)$adj.r.squared, 3), sep = "")
adj_r2_full

# Plot
ts<- ggplot(data = oisst_yr, aes(x = Year, y = Mean_SST, group = 1)) +
  geom_point(col = gmri_cols("gmri blue"), size = 3) +
  geom_path(col = gmri_cols("gmri blue"), lwd = 0.5) +
  geom_smooth(data = oisst_yr, aes(x = Year, y = Mean_SST), method = "lm", formula = y ~ x, col = gmri_cols("gmri blue"), size = 1, lty = "dashed", se = FALSE) +
  scale_x_continuous(breaks = pretty(oisst_yr$Year, n = 15)) +
  ylim(c(-2, 2)) +
  ylab("Annual SST Anomaly") +
  xlab("Year") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_text(aes(x = 1987, y = -1.75, label = my_formula_full), col = gmri_cols("gmri blue")) +
  geom_text(aes(x = 1997, y = -1.75, label = adj_r2_full), col = gmri_cols("gmri blue"))
ggsave(paste(proj_box_path, "Temp Results/TempTS.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")

# Add in 2004 warming trend
oisst_yr_2004<- oisst_yr %>%
  filter(., Year >= 2004)
lm_2004<- lm(Mean_SST ~ Year_Model, data = oisst_yr_2004)
my_formula_2004<- paste("SST Anom = ", as.numeric(round(coef(lm_2004)[1], 2)), " + ", as.numeric(round(coef(lm_2004)[2], 2)), "*Year", sep = "")
adj_r2_2004<- paste("Adj R2 = ", round(summary(lm_2004)$adj.r.squared, 3), sep = "")
adj_r2_2004

ts_with2004<- ggplot(data = oisst_yr, aes(x = Year, y = Mean_SST, group = 1)) +
  geom_point(col = gmri_cols("gmri blue"), size = 3) +
  geom_path(col = gmri_cols("gmri blue"), lwd = 0.5) +
  geom_smooth(data = oisst_yr, aes(x = Year, y = Mean_SST), method = "lm", formula = y ~ x, col = gmri_cols("gmri blue"), size = 1, lty = "dashed", se = FALSE) +
  geom_smooth(data = oisst_yr_2004, aes(x = Year, y = Mean_SST), method = "lm", formula = y ~ x, col = gmri_cols("orange"), size = 1, lty = "dashed", se = FALSE) +
  scale_x_continuous(breaks = pretty(oisst_yr$Year, n = 15)) +
  ylim(c(-2, 2)) +
  ylab("Annual SST Anomaly") +
  xlab("Year") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_text(aes(x = 1987, y = -1.85, label = my_formula_full), col = gmri_cols("gmri blue")) +
  geom_text(aes(x = 1997, y = -1.85, label = adj_r2_full), col = gmri_cols("gmri blue")) +
  geom_text(aes(x = 1987, y = -1.55, label = my_formula_2004), col = gmri_cols("orange")) +
  geom_text(aes(x = 1997, y = -1.55, label = adj_r2_2004), col = gmri_cols("orange"))
ggsave(paste(proj_box_path, "Temp Results/TempTS_2004.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")

# Bubble plot displaying the forecasting scenarios
row_ids<- seq(from = 2003, to = 2019, by = 1)
complete_rows<- length(seq(from = 2003, to = 2019, by = 1))
complete_cols<- length(seq(from = 1980, to = 2019, by = 1))
stop_start<- 2003

for(i in 1:complete_rows){
  row_temp<- rep("Training", length(seq(from = 1980, to = stop_start, by = 1)))
  
  # How much lead times are there for this row?
  letters_add<- 2019-stop_start
  if(letters_add > 0){
    row_temp<- c(row_temp, letters[1:min(letters_add, 10)])
    # Fill rest with NA
    if(length(row_temp) < complete_cols){
      row_temp<- c(row_temp, rep("k", complete_cols - length(row_temp)))
    }
  }
  
  if(i == 1){
    flower_plot_df<- row_temp
    stop_start<- stop_start+1
  } else {
    flower_plot_df<- rbind(flower_plot_df, row_temp)
    stop_start<- stop_start+1
  }
}

t<- data.frame(flower_plot_df)
rownames(t)<- c(paste("Training data through", seq(from = 2003, to = 2019, by = 1)))
colnames(t)<- seq(from = 1980, to = 2019, by = 1)
t$ForecastingChallenge<- rownames(t)

flower_plot_df_long<- t %>%
  gather(., "Year", "PlotGroup", -ForecastingChallenge) %>%
  filter(., ForecastingChallenge != "Training data through 2019")
flower_plot_df_long$Year<- as.numeric(flower_plot_df_long$Year)
flower_plot_df_long$PlotGroup<- factor(flower_plot_df_long$PlotGroup, levels = c(letters[1:10], "k", "Training"), labels = c("1y", "2y", "3y", "4y", "5y", "6y", "7y", "8y", "9y", "10y", "Not used", "Training"))
flower_plot_df_long$ForecastingChallenge<- factor(flower_plot_df_long$ForecastingChallenge, levels = rev(c(paste("Training data through", seq(from = 2003, to = 2018, by = 1)))), labels = rev(c(paste("Training data through", seq(from = 2003, to = 2018, by = 1)))))

fore_challenge_plot_out<- ggplot() + 
  geom_path(data = flower_plot_df_long, aes(x = Year, y = ForecastingChallenge, group = ForecastingChallenge), lwd = 0.5, alpha = 0.5) +
  geom_point(data = flower_plot_df_long, aes(x = Year, y = ForecastingChallenge, fill = PlotGroup, group = Year), pch = 21, size = 4) + 
  scale_fill_manual(name = "Forecast lead time", breaks = c("1y", "2y", "3y", "4y", "5y", "6y", "7y", "8y", "9y", "10y"), values = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026','#6a3d9a', "white", "#969696")) +
  scale_x_continuous(breaks = pretty(as.numeric(flower_plot_df_long$Year), n = 10)) +
  ylab("Forecasting Challenge") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste(proj_box_path, "Temp Results/ForeChallengePlot.jpg", sep = ""), dpi = 300, width = 12, height = 8, units = "in")

# Getting fitted results ----------------------------------------------------------
fore_challenges<- seq(from = 2004, to = 2019, by =  1)

res_root<- paste(shared.path(os.use = os_use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results"))

species_get<- c("HADDOCK")

#####
## VAST results
#####
res_folders<- list.files(res_root, pattern = "VAST", full.names = TRUE)
res_folders<- res_folders[which(grepl(species_get, res_folders))]

for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelfit.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Get fit settings
    t1_field<- t1$settings$FieldConfig
    t1_rho<- t1$settings$RhoConfig
    
    # Sample data...
    t1_sampdat<- t1$data_frame %>%
      select(., Lat_i, Lon_i, b_i, t_i)
    names(t1_sampdat)<- c("lat", "lon", "wtcpue", "year")
    t1_sampdat$surveyfact<- ifelse(grepl("FALL", res_folder_use[i]), "NEFSC_NEUSFall", "NEFSC_NEUSSpring")
    t1_sampdat$ForecastScenario<- paste(fore_challenges[j], " to 2019", sep = "")
    t1_sampdat$commonname<- sub('.*_', '', res_folder_use)
    t1_sampdat$presenceabsence<- ifelse(t1_sampdat$wtcpue > 0, 1, 0)
    
    # Model fit predictions
    t1_preds<- t1$Report
    
    # Combine em
    t1_out<- data.frame(t1_sampdat, "predicted.prob.presence" = t1_preds$R1_i, "predicted.bio" = t1_preds$D_i)
    
    # Subset
    t1_out$PredTF<- t1$data_list$PredTF_i
    t1_out_sub<- t1_out[t1_out$PredTF == 1, ]
    
    # Arrange things
    t1_out_nest<- t1_out_sub %>%
      select(., commonname, surveyfact, ForecastScenario, lat, lon, year, wtcpue, presenceabsence, predicted.prob.presence, predicted.bio) %>%
      group_by(commonname, ForecastScenario, surveyfact) %>%
      nest()
    
    if(i == 1 && j == 1){
      vast_fits_out<- t1_out_nest
      print(paste(res_use, " is done!", sep = ""))
    } else {
      vast_fits_out<- bind_rows(vast_fits_out, t1_out_nest)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}

# Name clean up
names(vast_fits_out)[4]<- "PredictionDF"
vast_fits_out$Model<- rep("VAST", nrow(vast_fits_out))

# Now getting the results we want. For each row, we want to go in and get the Prediction statistics at intervals from 1 to 10 years where possible...
forecast_lead_time<- data.frame(expand.grid("surveyfact" = c("NEFSC_NEUSFall", "NEFSC_NEUSSpring"), "ForecastScenario" = unique(vast_fits_out$ForecastScenario), "LeadTime" = seq(from = 0, to = 9, by = 1)))
vast_fits_out<- vast_fits_out %>%
  left_join(., forecast_lead_time)

#####
## GAM -- Naive
#####
res_folders<- list.files(res_root, pattern = "GAM_Naive", full.names = TRUE) 
res_folders<- res_folders[which(grepl(species_get, res_folders))]

for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelpred.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Sample data...
    t1_sampdat<- t1 %>%
      select(., DECDEG_BEGLAT, DECDEG_BEGLON, SUM_BIOMASS, EST_YEAR)
    names(t1_sampdat)<- c("lat", "lon", "wtcpue", "year")
    t1_sampdat$surveyfact<- ifelse(grepl("FALL", res_folder_use[i]), "NEFSC_NEUSFall", "NEFSC_NEUSSpring")
    t1_sampdat$ForecastScenario<- paste(fore_challenges[j], " to 2019", sep = "")
    t1_sampdat$commonname<- sub('.*_', '', res_folder_use)
    t1_sampdat$presenceabsence<- ifelse(t1_sampdat$wtcpue > 0, 1, 0)
    
    # Model predictions
    t1_out<- data.frame(t1_sampdat, "predicted.prob.presence" = t1$PRED_PRESENCE, "predicted.bio" = t1$PRED_BIO)
    
    # Arrange things
    t1_out_nest<- t1_out %>%
      select(., commonname, surveyfact, ForecastScenario, lat, lon, year, wtcpue, presenceabsence, predicted.prob.presence, predicted.bio) %>%
      group_by(commonname, ForecastScenario, surveyfact) %>%
      nest()
    
    if(i == 1 && j == 1){
      gam_fits_out<- t1_out_nest
      print(paste(res_use, " is done!", sep = ""))
    } else {
      gam_fits_out<- bind_rows(gam_fits_out, t1_out_nest)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}

# Name clean up
names(gam_fits_out)[4]<- "PredictionDF"
gam_fits_out$Model<- rep("GAM_Naive", nrow(gam_fits_out))
forecast_lead_time<- data.frame(expand.grid("surveyfact" = c("NEFSC_NEUSFall", "NEFSC_NEUSSpring"), "ForecastScenario" = unique(gam_fits_out$ForecastScenario), "LeadTime" = seq(from = 0, to = 9, by = 1)))
gam_fits_out_naive<- gam_fits_out %>%
  left_join(., forecast_lead_time)

#####
## GAM -- ForeYearM1
#####
res_folders<- list.files(res_root, pattern = "GAM_Fore_YearM_1", full.names = TRUE) 
res_folders<- res_folders[which(grepl("species_get", res_folders))]

for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelpred.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Sample data...
    t1_sampdat<- t1 %>%
      select(., DECDEG_BEGLAT, DECDEG_BEGLON, SUM_BIOMASS, EST_YEAR)
    names(t1_sampdat)<- c("lat", "lon", "wtcpue", "year")
    t1_sampdat$surveyfact<- ifelse(grepl("FALL", res_folder_use[i]), "NEFSC_NEUSFall", "NEFSC_NEUSSpring")
    t1_sampdat$ForecastScenario<- paste(fore_challenges[j], " to 2019", sep = "")
    t1_sampdat$commonname<- sub('.*_', '', res_folder_use)
    t1_sampdat$presenceabsence<- ifelse(t1_sampdat$wtcpue > 0, 1, 0)
    
    # Model predictions
    t1_out<- data.frame(t1_sampdat, "predicted.prob.presence" = t1$PRED_PRESENCE, "predicted.bio" = t1$PRED_BIO)
    
    # Arrange things
    t1_out_nest<- t1_out %>%
      select(., commonname, surveyfact, ForecastScenario, lat, lon, year, wtcpue, presenceabsence, predicted.prob.presence, predicted.bio) %>%
      group_by(commonname, ForecastScenario, surveyfact) %>%
      nest()
    
    if(i == 1 && j == 1){
      gam_fits_out<- t1_out_nest
      print(paste(res_use, " is done!", sep = ""))
    } else {
      gam_fits_out<- bind_rows(gam_fits_out, t1_out_nest)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}

# Name clean up
names(gam_fits_out)[4]<- "PredictionDF"
gam_fits_out$Model<- rep("GAM_YearM1", nrow(gam_fits_out))
forecast_lead_time<- data.frame(expand.grid("surveyfact" = c("NEFSC_NEUSFall", "NEFSC_NEUSSpring"), "ForecastScenario" = unique(gam_fits_out$ForecastScenario), "LeadTime" = seq(from = 0, to = 9, by = 1)))
gam_fits_out_yearm1<- gam_fits_out %>%
  left_join(., forecast_lead_time)

#####
## GAM -- ForeYearM2
#####
res_folders<- list.files(res_root, pattern = "GAM_Fore_YearM_2", full.names = TRUE) 
res_folders<- res_folders[which(grepl(species_get, res_folders))]

for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelpred.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Sample data...
    t1_sampdat<- t1 %>%
      select(., DECDEG_BEGLAT, DECDEG_BEGLON, SUM_BIOMASS, EST_YEAR)
    names(t1_sampdat)<- c("lat", "lon", "wtcpue", "year")
    t1_sampdat$surveyfact<- ifelse(grepl("FALL", res_folder_use[i]), "NEFSC_NEUSFall", "NEFSC_NEUSSpring")
    t1_sampdat$ForecastScenario<- paste(fore_challenges[j], " to 2019", sep = "")
    t1_sampdat$commonname<- sub('.*_', '', res_folder_use)
    t1_sampdat$presenceabsence<- ifelse(t1_sampdat$wtcpue > 0, 1, 0)
    
    # Model predictions
    t1_out<- data.frame(t1_sampdat, "predicted.prob.presence" = t1$PRED_PRESENCE, "predicted.bio" = t1$PRED_BIO)
    
    # Arrange things
    t1_out_nest<- t1_out %>%
      select(., commonname, surveyfact, ForecastScenario, lat, lon, year, wtcpue, presenceabsence, predicted.prob.presence, predicted.bio) %>%
      group_by(commonname, ForecastScenario, surveyfact) %>%
      nest()
    
    if(i == 1 && j == 1){
      gam_fits_out<- t1_out_nest
      print(paste(res_use, " is done!", sep = ""))
    } else {
      gam_fits_out<- bind_rows(gam_fits_out, t1_out_nest)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}

# Name clean up
names(gam_fits_out)[4]<- "PredictionDF"
gam_fits_out$Model<- rep("GAM_YearM2", nrow(gam_fits_out))
forecast_lead_time<- data.frame(expand.grid("surveyfact" = c("NEFSC_NEUSFall", "NEFSC_NEUSSpring"), "ForecastScenario" = unique(gam_fits_out$ForecastScenario), "LeadTime" = seq(from = 0, to = 9, by = 1)))
gam_fits_out_yearm2<- gam_fits_out %>%
  left_join(., forecast_lead_time)

#####
## Combine and calculate prediction statistics
#####
all_fits_out<- bind_rows(vast_fits_out, gam_fits_out_naive, gam_fits_out_yearm1, gam_fits_out_yearm2)

# Prediction skill statistics...
all_stats<- all_fits_out %>%
  mutate(., "AUC" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(auc_func, NA)),
         "MaxKappa" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(maxkappa_func, NA)),
         "Precision" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa, LeadTime = LeadTime), possibly(precision_func, NA)),
         "Specificity" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa, LeadTime = LeadTime), possibly(spec_func, NA)),
         "F1Measure" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa, LeadTime = LeadTime), possibly(fmeasure_func, NA)),
         "CorrCoeff.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(corr_coeff_func, NA)),
         "CoeffDet.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(coeff_det_func, NA)),
         "RMSE.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(rmse_func, NA)),
         "SDBias.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", LeadTime = LeadTime), possibly(sd_bias_func, NA)),
         "CorrCoeff.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(corr_coeff_func, NA)),
         "CoeffDet.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(coeff_det_func, NA)),
         "RMSE.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(rmse_func, NA)),
         "MAE.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(mae_func, NA)),
         "MASE.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(mase_func, NA)),
         "SDBias.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio", LeadTime = LeadTime), possibly(sd_bias_func, NA)))


# Formatting...
all_stats_long<- all_stats %>%
  select(., commonname, surveyfact, ForecastScenario, LeadTime, Model, CorrCoeff.Biomass:SDBias.Biomass) %>%
  gather(., key = "PredStatVariable", value = "PredStatValue", -commonname, -surveyfact, -ForecastScenario, - LeadTime, -Model)
all_stats_long$TrainingDataExtent<-  as.numeric(stri_extract_first(all_stats_long$ForecastScenario, regex="\\w+")) - 1980
all_stats_long$ForecastLeadTimePlot<- factor(all_stats_long$LeadTime + 1, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), labels = c("1y", "2y", "3y", "4y", "5y", "6y", "7y", "8y", "9y", "10y"))
all_stats_long$surveyfactplot<- factor(all_stats_long$surveyfact, levels = c("NEFSC_NEUSFall", "NEFSC_NEUSSpring"), labels = c("Fall", "Spring"))

# Summaries...
# all_stats_long<- all_stats_long %>%
#   filter(., Model %in% c("GAM_Naive", "GAM_YearM1", "VAST"))
all_stats_long<- all_stats_long %>%
  filter(., Model %in% c("GAM_Naive", "VAST"))
## Plots...
# Ribbon plot with multiple models 
all_scenarios<- TRUE
ribbon<- FALSE
focal_scenarios<- NULL
plot_name<- NULL
if(all_sencarios & ribbon){
  colors_all<- c("gray", "#1b9e77", "#d95f02")
  colors_use<- colors_all[1:length(unique(all_stats_long$Model))]
  corr_coeff_dat<- all_stats_long %>%
    filter(., PredStatVariable == "CorrCoeff.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  CorrCoeff_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "Correlation coefficient\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CorrCoeffBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  coeff_det_dat<- all_stats_long %>%
    filter(., PredStatVariable == "CoeffDet.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  CoeffDet_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "Coefficient of determination\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CoeffDetBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  rmse_dat<- all_stats_long %>%
    filter(., PredStatVariable == "RMSE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  RMSE_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = rmse_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "RMSE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/RMSEBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  mae_dat<- all_stats_long %>%
    filter(., PredStatVariable == "MAE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  MAE_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = mae_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "MAE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/MAEBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  mase_dat<- all_stats_long %>%
    filter(., PredStatVariable == "MASE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  MASE_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = mase_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = mase_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = mase_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "MASE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/MASEBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  
  sd_bias_dat<- all_stats_long %>%
    filter(., PredStatVariable == "SDBias.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  SDBias_Bio_Plot<- ggplot() + 
    #geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_ribbon(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "SD Ratio\n Predicted:Observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/SDBiasBiomass.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
}

# Error bar plot, with 2012
all_scenarios<- TRUE
ribbon<- FALSE
focal_scenarios<- NULL
plot_name<- NULL

scenarios_2012<- data.frame("ForecastScenario" = unique(all_stats_long$ForecastScenario)[1:9],
                            "LeadTime" = c(8, 7, 6, 5, 4, 3, 2, 1, 0))
stats_2012<- scenarios_2012 %>%
  left_join(., all_stats_long, by = c("ForecastScenario" = "ForecastScenario", "LeadTime" = "LeadTime"))
if(all_scenarios & !ribbon){
  colors_all<- c("#00608A", "#EA4F12")
  colors_use<- colors_all[1:length(unique(all_stats_long$Model))]
  stats_2012_only<- all_stats_long %>%
    filter(., )
  corr_coeff_dat<- all_stats_long %>%
    filter(., PredStatVariable == "CorrCoeff.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  corr_coeff_2012<- stats_2012 %>%
    filter(., PredStatVariable == "CorrCoeff.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  
  CorrCoeff_Bio_Plot<- ggplot() + 
    geom_errorbar(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    #geom_ribbon(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, fill = Model, group = Model), alpha = 0.4, width = .2) +
    geom_point(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = corr_coeff_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "Correlation coefficient\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CorrCoeffBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  coeff_det_dat<- all_stats_long %>%
    filter(., PredStatVariable == "CoeffDet.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  coeff_det_2012<- stats_2012 %>%
    filter(., PredStatVariable == "CoeffDet.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  CoeffDet_Bio_Plot<- ggplot() + 
    geom_errorbar(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_point(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = coeff_det_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "Coefficient of determination\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CoeffDetBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  rmse_dat<- all_stats_long %>%
    filter(., PredStatVariable == "RMSE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  rmse_2012<- stats_2012 %>%
    filter(., PredStatVariable == "RMSE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  RMSE_Bio_Plot<- ggplot() + 
    geom_errorbar(data = rmse_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_point(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = rmse_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "RMSE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/RMSEBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  mae_dat<- all_stats_long %>%
    filter(., PredStatVariable == "MAE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  mae_2012<- stats_2012 %>%
    filter(., PredStatVariable == "MAE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  MAE_Bio_Plot<- ggplot() + 
    geom_errorbar(data = mae_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_point(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = mae_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "MAE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/MAEBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  mase_dat<- all_stats_long %>%
    filter(., PredStatVariable == "MASE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  mase_2012<- stats_2012 %>%
    filter(., PredStatVariable == "MASE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  MASE_Bio_Plot<- ggplot() + 
    geom_errorbar(data = mase_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_point(data = mase_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = mase_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = mase_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "MASE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/MASEBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  
  sd_bias_dat<- all_stats_long %>%
    filter(., PredStatVariable == "SDBias.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  sd_bias_2012<- stats_2012 %>%
    filter(., PredStatVariable == "SDBias.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  
  SDBias_Bio_Plot<- ggplot() + 
    geom_errorbar(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, ymin = Mean-SD, ymax = Mean+SD, color = Model), width = .2, alpha = 0.4) +
    geom_point(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 21, size = 3) +
    geom_point(data = sd_bias_2012, aes(x = ForecastLeadTimePlot, y = Mean, fill = Model), shape = 23, size = 3) +
    geom_path(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Model", values = colors_use) +
    scale_fill_manual(name = "Model", values = colors_use) +
    scale_y_continuous(name = "SD Ratio\n Predicted:Observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/SDBiasBiomassWith2012.jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
}

if(!all_scenarios){
  corr_coeff_dat<- all_stats_long %>%
    filter(., ForecastScenario %in% focal_scenarios) %>%
    filter(., PredStatVariable == "CorrCoeff.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  CorrCoeff_Bio_Plot<- ggplot() + 
    geom_point(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = corr_coeff_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = c("#00608A", "#EA4F12")) +
    scale_y_continuous(name = "Correlation coefficient\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CorrCoeffBiomass", plot_name, ".jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  coeff_det_dat<- all_stats_long %>%
    filter(., ForecastScenario %in% focal_scenarios) %>%
    filter(., PredStatVariable == "CoeffDet.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  CoeffDet_Bio_Plot<- ggplot() + 
    geom_point(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = coeff_det_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = c("#00608A", "#EA4F12")) +
    scale_y_continuous(name = "Coefficient of determination\n predicted vs. observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/CoeffDetBiomass", plot_name, ".jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  rmse_dat<- all_stats_long %>%
    filter(., ForecastScenario %in% focal_scenarios) %>%
    filter(., PredStatVariable == "RMSE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  RMSE_Bio_Plot<- ggplot() + 
    geom_point(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = rmse_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = c("#00608A", "#EA4F12")) +
    scale_y_continuous(name = "RMSE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/RMSEBiomass", plot_name, ".jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  mae_dat<- all_stats_long %>%
    filter(., ForecastScenario %in% focal_scenarios) %>%
    filter(., PredStatVariable == "MAE.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  MAE_Bio_Plot<- ggplot() + 
    geom_point(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = mae_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = c("#00608A", "#EA4F12")) +
    scale_y_continuous(name = "MAE (kg/tow)") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/MAEBiomass", plot_name, ".jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
  
  
  sd_bias_dat<- all_stats_long %>%
    filter(., ForecastScenario %in% focal_scenarios) %>%
    filter(., PredStatVariable == "SDBias.Biomass") %>%
    group_by(., commonname, surveyfactplot, ForecastLeadTimePlot, Model) %>%
    summarize_at("PredStatValue", .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  SDBias_Bio_Plot<- ggplot() + 
    geom_point(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, color = Model), size = 3) +
    geom_path(data = sd_bias_dat, aes(x = ForecastLeadTimePlot, y = Mean, group = Model, color = Model, order = ForecastLeadTimePlot), lwd = 0.5) +
    scale_color_manual(name = "Model", values = c("#00608A", "#EA4F12")) +
    scale_y_continuous(name = "SD Ratio\n Predicted:Observed") +
    xlab("Forecast challenge scenario lead time") +
    facet_wrap(~surveyfactplot) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste(proj_box_path, "Temp Results/SDBiasBiomass", plot_name, ".jpg", sep = ""), dpi = 300, width = 12, height = 6, units = "in")
}

#####
## Effect plots
#####
# Getting fitted results ----------------------------------------------------------
fore_challenges<- seq(from = 2004, to = 2019, by =  1)

res_root<- paste(shared.path(os.use = os_use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results"))

#####
## VAST results
#####
res_folders<- list.files(res_root, pattern = "VAST", full.names = TRUE)
res_folders<- res_folders[which(grepl(species_get, res_folders))]
for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelfit.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Must add data-frames to global environment (hope to fix in future)
    if(is.null(t1$effects)){
      t1$effects<- list()
      if(!is.null(t1$catchability_data)) {
        catchability_data_full = data.frame(t1$catchability_data, 
                                            linear_predictor = 0)
        Q1_formula_full = update.formula(t1$Q1_formula, linear_predictor ~ 
                                           . + 0)
        call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
        Q2_formula_full = update.formula(t1$Q2_formula, linear_predictor ~ 
                                           . + 0)
        call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
        t1$effects = c(Return$effects, list(call_Q1 = call_Q1, 
                                                 call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
      }
      if (!is.null(t1$covariate_data)) {
        covariate_data_full = data.frame(t1$covariate_data, linear_predictor = 0)
        X1_formula_full = update.formula(t1$X1_formula, linear_predictor ~ 
                                           . + 0)
        call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
        X2_formula_full = update.formula(t1$X2_formula, linear_predictor ~ 
                                           . + 0)
        call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
        t1$effects = c(t1$effects, list(call_X1 = call_X1, 
                                                  call_X2 = call_X2, covariate_data_full = covariate_data_full))
      }
    }
   
    # Define formula.
    X1_formula = t1$X1_formula
    X2_formula = t1$X2_formula
    
    # Get effects...
    depth_plot<- data.frame(Effect.fit_model(focal.predictors = c("AVGDEPTH"), mod = t1, which_formula = "X2", xlevels = 100))
    depth_plot$Variable<- rep("Depth", nrow(depth_plot))
    names(depth_plot)[1]<- "VarValue"
    sst_plot<- data.frame(Effect.fit_model(focal.predictors = c("SODA_SST"), mod = t1, which_formula = "X2", xlevels = 100))
    sst_plot$Variable<- rep("SST", nrow(sst_plot))
    names(sst_plot)[1]<- "VarValue"
    bt_plot<-  data.frame(Effect.fit_model(focal.predictors = c("SODA_BT"), mod = t1, which_formula = "X2", xlevels = 100))
    bt_plot$Variable<- rep("BT", nrow(bt_plot))
    names(bt_plot)[1]<- "VarValue"
    
    plot_all<- bind_rows(depth_plot, sst_plot, bt_plot)
    plot_all$Variable<- factor(plot_all$Variable, levels = c("Depth", "SST", "BT"))
    plot_all$Model<- rep("VAST", nrow(plot_all))
    plot_all$ForeChallenge<- rep(fore_challenge_use, nrow(plot_all))
    plot_all$Season<- ifelse(grepl("SPRING", res_use), "Spring", "Fall")
    
    if(i == 1 & j == 1){
      vast_coeff_plot<- plot_all
      print(paste(res_use, " is done!", sep = ""))
    } else {
      vast_coeff_plot<- bind_rows(vast_coeff_plot, plot_all)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}

#####
## GAM results
#####
res_folders<- list.files(res_root, pattern = "GAM_Naive", full.names = TRUE)
res_folders<- res_folders[which(grepl("ATLANTIC COD", res_folders))]
for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelfit.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)$LogBiomass
    
    # Get effects...
    gam_smooth_data<- plot(t1, plot.me = FALSE)
    
    depth_plot<- data.frame(effects(focal.predictors = c("AVGDEPTH"), mod = t1, which_formula = "X2", xlevels = 100, sources = args))
    depth_plot$Variable<- rep("Depth", nrow(depth_plot))
    names(depth_plot)[1]<- "VarValue"
    sst_plot<- data.frame(Effect.fit_model(focal.predictors = c("SODA_SST"), mod = t1, which_formula = "X2", xlevels = 100))
    sst_plot$Variable<- rep("SST", nrow(sst_plot))
    names(sst_plot)[1]<- "VarValue"
    bt_plot<-  data.frame(Effect.fit_model(focal.predictors = c("SODA_BT"), mod = t1, which_formula = "X2", xlevels = 100))
    bt_plot$Variable<- rep("BT", nrow(bt_plot))
    names(bt_plot)[1]<- "VarValue"
    
    plot_all<- bind_rows(depth_plot, sst_plot, bt_plot)
    plot_all$Variable<- factor(plot_all$Variable, levels = c("Depth", "SST", "BT"))
    plot_all$Model<- rep("VAST", nrow(plot_all))
    plot_all$ForeChallenge<- rep(fore_challenge_use, nrow(plot_all))
    plot_all$Season<- ifelse(grepl("SPRING", res_use), "Spring", "Fall")
    
    if(i == 1 & j == 1){
      vast_coeff_plot<- plot_all
      print(paste(res_use, " is done!", sep = ""))
    } else {
      vast_coeff_plot<- bind_rows(vast_coeff_plot, plot_all)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}



fore_challenges<- seq(from = 2004, to = 2019, by =  1)

res_root<- paste(shared.path(os.use = os_use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results"))


# Prediction Maps ---------------------------------------------------------
#####
## VAST results
#####
res_folders<- list.files(res_root, pattern = "VAST", full.names = TRUE)
res_folders<- res_folders[which(grepl(species_get, res_folders))]
for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelfit.rds")
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    # Get biomass index stuff
    pred_temp<- as_tibble(data.frame("Scenario" = fore_challenge_use, "Season" = ifelse(grepl("SPRING", res_use), "Spring", "Fall")))
    pred_temp<- bind_cols(pred_temp, as_tibble_col(list(t1), column_name = "Mod")) %>%
      mutate(., "Biomass_Index" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), biomass_index_zz),
             "Density" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), Variable = list("Density"), out_folder = list(res_folder_use), var_maps_zz),
             "Epsilon_1" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), Variable = list("Epsilon_1"), out_folder = list(res_folder_use), var_maps_zz),
             "Epsilon_2" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), Variable = list("Epsilon_2"), out_folder = list(res_folder_use), var_maps_zz),
             "Omega_1" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), Variable = list("Omega_1"), out_folder = list(res_folder_use), var_maps_zz),
             "Omega_2" = pmap(list(Scenario = Scenario, Season = Season, Species = list(species_get), Fit = Mod), Variable = list("Omega_2"), out_folder = list(res_folder_use), var_maps_zz))
    
    if(i == 1 & j == 1){
      vast_preds_all<- pred_temp
      print(paste(res_use, " is done!", sep = ""))
    } else {
      vast_preds_all<- bind_rows(vast_preds_all, pred_temp)
      print(paste(res_use, " is done!", sep = ""))
    }
  }
}
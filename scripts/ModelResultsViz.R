#####
## NSF C-Accel Forecasting Challenge (11/12/2020): Visualizing results
#####
# Preliminaries ----------------------------------------------------------
library(tidyverse)
library(gmRi)
os.use<- "unix"

# Getting haul information ------------------------------------------------
dat_all<-  read.csv(paste(fish_dat_path, "ForecastingChallengeModelData.csv", sep = "")) 

# Going to want some other information from the dat.all as these just have haulid
hauls<- dat_all %>%
  dplyr::select(., ROW_ID, EST_YEAR, EST_MONTH, DECDEG_BEGLON, DECDEG_BEGLAT) %>%
  distinct() 
names(hauls)<- c("haulid", "year", "month", "lon", "lat")

# Getting fitted results ----------------------------------------------------------
fore_challenges<- seq(from = 2004, to = 2019, by =  1)
fore_challenges<- fore_challenges[c(1, 6, 11, 15)]
res_root<- paste(shared.path(os.use = os.use, group = "Mills Lab", folder = "Projects/ForecastingChallenge/Temp Results"))
res_folders<- list.files(res_root, pattern = "VAST", full.names = TRUE)

for(i in seq_along(res_folders)){
  # Get folder
  res_folder_use<- res_folders[i]
  
  # Get vector of potential result files
  res_files<- list.files(res_folder_use, full.names = TRUE, pattern = "modelfit.rds")
  
  for(j in seq_along(fore_challenges)){
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenges[j], res_files))]
    
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
    
    if(j == 1){
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

# Prediction skill statistics...
# Going to need things in a slightly different format...
# VAST models
vast_preds<- vast_fits_out %>%
  dplyr::select(commonname, ForecastScenario, PredictionDF) %>%
  ungroup() %>%
  unnest(cols = "PredictionDF")


# Subsetting...
vast_preds_2004<- vast_preds %>%
  dplyr::filter(., ForecastScenario == "2004 to 2019") %>%
  dplyr::select(., lat, lon, year, commonname, ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  mutate("Species_ForecastScenario" = paste(commonname, ForecastScenario, sep = "_")) %>%
  dplyr::select(lat, lon, year, Species_ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  gather(., "Response", "Prediction", -lat, -lon, -year, -Species_ForecastScenario) %>%
  pivot_wider(names_from = c(Species_ForecastScenario, Response), values_from = Prediction)

vast_preds_2009<- vast_preds %>%
  dplyr::filter(., ForecastScenario == "2009 to 2019") %>%
  dplyr::select(., lat, lon, year, commonname, ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  mutate("Species_ForecastScenario" = paste(commonname, ForecastScenario, sep = "_")) %>%
  dplyr::select(lat, lon, year, Species_ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  gather(., "Response", "Prediction", -lat, -lon, -year, -Species_ForecastScenario) %>%
  pivot_wider(names_from = c(Species_ForecastScenario, Response), values_from = Prediction)

vast_preds_2014<- vast_preds %>%
  dplyr::filter(., ForecastScenario == "2014 to 2019") %>%
  dplyr::select(., lat, lon, year, commonname, ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  mutate("Species_ForecastScenario" = paste(commonname, ForecastScenario, sep = "_")) %>%
  dplyr::select(lat, lon, year, Species_ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  gather(., "Response", "Prediction", -lat, -lon, -year, -Species_ForecastScenario) %>%
  pivot_wider(names_from = c(Species_ForecastScenario, Response), values_from = Prediction)

vast_preds_2018<- vast_preds %>%
  dplyr::filter(., ForecastScenario == "2018 to 2019") %>%
  dplyr::select(., lat, lon, year, commonname, ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  mutate("Species_ForecastScenario" = paste(commonname, ForecastScenario, sep = "_")) %>%
  dplyr::select(lat, lon, year, Species_ForecastScenario, predicted.prob.presence, predicted.bio) %>%
  gather(., "Response", "Prediction", -lat, -lon, -year, -Species_ForecastScenario) %>%
  pivot_wider(names_from = c(Species_ForecastScenario, Response), values_from = Prediction)

vast_preds_all<- vast_preds_2004 %>%
  left_join(., vast_preds_2009) %>%
  left_join(., vast_preds_2014) %>%
  left_join(., vast_preds_2018)

# Add haul info...
hauls<- vast_preds %>%
  dplyr::select(., lat, lon, year, wtcpue, presenceabsence)
vast_preds_all<- hauls %>%
  left_join(., vast_preds_all)
vast_preds_all$Model<- rep("VAST", nrow(vast_preds_all))
all_out<- vast_preds_all %>%
  dplyr::select(., lat, lon, year, Model, colnames(vast_preds_all)[5:ncol(vast_preds_all)-1])
all_out$ID<- seq(from = 1, to = nrow(all_out), by = 1)

## Now, likely going to work in the long format....
# First the observations....
obs_keep_cols<- c("ID", colnames(all_out)[which(grepl("wtcpue", colnames(all_out)))], colnames(all_out)[which(grepl("presenceabsence", colnames(all_out)))])
dat_obs<- all_out %>%
  dplyr::select(., one_of(obs_keep_cols)) %>%
  distinct() %>%
  pivot_longer(-ID, names_to = "Species_Response", values_to = "Observation") %>%
  #separate(., "Species_Response", into = c("Species", "Response"), sep = "_") %>%
  pivot_wider(names_from = Species_Response, values_from = Observation)

vast_keep_cols<- c("ID", colnames(all_out)[which(grepl("predicted.prob.presence", colnames(all_out)))], colnames(all_out)[which(grepl("predicted.bio", colnames(all_out)))])
vast_preds_w<- all_out %>%
  dplyr::select(., one_of(vast_keep_cols)) %>%
  distinct() %>%
  pivot_longer(., -ID, names_to = "Species_ForecastScenario", values_to = "Predictions") %>%
  separate(., "Species_ForecastScenario", into = c("Species", "ForecastScenario", "ModeledResponse"), sep = "_") %>%
  pivot_wider(names_from = c(ModeledResponse), values_from = Predictions)

vast_preds_w<- dat_obs %>%
  left_join(., vast_preds_w, by = c("ID")) %>%
  distinct() %>%
  drop_na(predicted.prob.presence, predicted.bio)

vast_stats<- vast_preds_w %>%
  group_by(., Species, ForecastScenario) %>%
  nest() %>%
  rename(PredictionDF = data) %>%
  mutate(., "PredRange.ProbPresence" = map2(PredictionDF, "predicted.prob.presence", possibly(pred_ranges_func, NA)),
         "AUC" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(auc_func, NA)),
         "MaxKappa" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(maxkappa_func, NA)),
         "Precision" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa), possibly(precision_func, NA)),
         "Specificity" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa), possibly(spec_func, NA)),
         "F1Measure" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence", maxkappa = MaxKappa), possibly(fmeasure_func, NA)),
         "CorrCoeff.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(corr_coeff_func, NA)),
         "CoeffDet.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(coeff_det_func, NA)),
         "RMSE.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(rmse_func, NA)),
         "SDBias.ProbPresence" = pmap_dbl(list(df = PredictionDF, obs = "presenceabsence", mod = "predicted.prob.presence"), possibly(sd_bias_func, NA)),
         "PredRange.Bio" = map2(PredictionDF, "predicted.bio", possibly(pred_ranges_func, NA)),
         "CorrCoeff.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio"), possibly(corr_coeff_func, NA)),
         "CoeffDet.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio"), possibly(coeff_det_func, NA)),
         "RMSE.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio"), possibly(rmse_func, NA)),
         "SDBias.Biomass" = pmap_dbl(list(df = PredictionDF, obs = "wtcpue", mod = "predicted.bio"), possibly(sd_bias_func, NA)))

#####
## Delta log normal wrapper functions
#####

deltalog_mod_fit_naive<- function(season, comname, data, fore_challenge, covs = c("EST_YEAR", "AVGDEPTH", "SODA_BT", "SODA_SST"), out_folder){
  
  # For debugging
  if(FALSE){
    season = dat_nest$SEASON[[1]]
    comname = dat_nest$COMNAME[[1]]
    data = dat_nest$data[[1]]
    fore_challenge = dat_nest$fore_challenge[[1]]
    covs = c("AVGDEPTH", "SODA_BT", "SODA_SST")
    out_folder = out_folder
  }
  
  # Some things for where to store the results
  fore_dates<- paste(fore_challenge, "to2019", sep = "")
  
  # Create output folder
  outfolder<- paste(out_folder, paste("GAM_Naive", season, comname, sep = "_"), sep = "")
  if(!file.exists(outfolder)){
    dir.create(outfolder)
  }
  
  # Create sub-folder
  outfile<- paste(outfolder, "/", fore_dates, sep = "")
  
  # Get specific data
  dat_run<- data
  
  # Rescale covariates to SD of 1
  dat_run<- dat_run %>%
    mutate_at(., {{covs}}, vast_scale_func, type = "AJA")
  
  # Model fitting data only
  dat_mod<- dat_run %>%
    filter(., EST_YEAR <  fore_challenge) %>%
    drop_na({{covs}})
  
  # Fit presence/absence and then the log(biomass) model
  pa_mod<- gam(PRESENCE ~ s(AVGDEPTH, bs = "cs") + s(SODA_BT, bs = "cs") + s(SODA_SST, bs = "cs"), family = "binomial", method = "REML", data = dat_mod)
  
  # Now, log positive biomass model
  dat_mod_logbio<- dat_mod %>%
    filter(., PRESENCE == 1)
  logbio_mod<- gam(LOGBIO ~ s(AVGDEPTH, bs = "cs") + s(SODA_BT, bs = "cs") + s(SODA_SST, bs = "cs"), family = "gaussian", method = "REML", data = dat_mod_logbio)
  
  # Save and return it
  mods_out<- list("PresenceAbsenceModel" = pa_mod, "LogBiomass" = logbio_mod)
  saveRDS(mods_out, file = paste(outfile, "modelfit.rds", sep = ""))
  return(mods_out)
  
}

deltalog_mod_fit_fore<- function(season, comname, data, fore_challenge, covs = c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST"), covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST"), out_folder, year_m = 1){
  
  # For debugging
  if(FALSE){
    season = dat_nest$SEASON[[1]]
    comname = dat_nest$COMNAME[[1]]
    data = dat_nest$data[[1]]
    fore_challenge = dat_nest$fore_challenge[[1]]
    covs = c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST")
    covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")
    out_folder = out_folder
  }
  
  # Some things for where to store the results
  fore_dates<- paste(fore_challenge, "to2019", sep = "")
  
  # Create output folder
  outfolder<- paste(out_folder, paste("GAM_ForeYearM", year_m, season, comname, sep = "_"), sep = "")
  if(!file.exists(outfolder)){
    dir.create(outfolder)
  }
  
  # Create sub-folder
  outfile<- paste(outfolder, "/", fore_dates, sep = "")
  
  # Get specific data
  dat_run<- data
  
  # Rescale covariates to SD of 1
  dat_run<- dat_run %>%
    mutate_at(., {{covs_scale}}, vast_scale_func, type = "AJA")
  
  # Model fitting data only
  dat_mod<- dat_run %>%
    filter(., EST_YEAR <  fore_challenge) %>%
    drop_na({{covs}})
  
  # Fit presence/absence and then the log(biomass) model
  pa_mod<- gam(PRESENCE ~ s(EST_YEAR, m = year_m) + s(DECDEG_BEGLON, DECDEG_BEGLAT) + s(AVGDEPTH, bs = "cs") + s(SODA_BT, bs = "cs") + s(SODA_SST, bs = "cs"), select = TRUE, family = "binomial", method = "REML", data = dat_mod)
  
  # Now, log positive biomass model
  dat_mod_logbio<- dat_mod %>%
    filter(., PRESENCE == 1)
  logbio_mod<- gam(LOGBIO ~ s(EST_YEAR, m = year_m) + s(DECDEG_BEGLON, DECDEG_BEGLAT) + s(AVGDEPTH, bs = "cs") + s(SODA_BT, bs = "cs") + s(SODA_SST, bs = "cs"), select = TRUE, family = "gaussian", method = "REML", data = dat_mod_logbio)
  
  # Save and return it
  mods_out<- list("PresenceAbsenceModel" = pa_mod, "LogBiomass" = logbio_mod)
  saveRDS(mods_out, file = paste(outfile, "modelfit.rds", sep = ""))
  return(mods_out)
  
}

deltalog_mod_predict<- function(season, comname, data, mod_list, fore_challenge, covs = c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST"), covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST"), mod_form, out_folder){
  
  # For debugging
  if(FALSE){
    season = dat_nest$SEASON[[1]]
    comname = dat_nest$COMNAME[[1]]
    data = dat_nest$data[[1]]
    mod_list = dat_nest$GAM_Fits_Fore_M1[[1]]
    fore_challenge = dat_nest$fore_challenge[[1]]
    covs = c("EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH", "SODA_BT", "SODA_SST") 
    covs_scale = c("AVGDEPTH", "SODA_BT", "SODA_SST")
    mod_form = "YearM1"
    out_folder = out_folder
  }
  
  # Some things for where to store the results
  fore_dates<- paste(fore_challenge, "to2019", sep = "")
  
  # Create output folder
  outfolder<- paste(out_folder, paste("GAM", mod_form, season, comname, sep = "_"), sep = "")
  if(!file.exists(outfolder)){
    dir.create(outfolder)
  }
  
  # Create sub-folder
  outfile<- paste(outfolder, "/", fore_dates, sep = "")
  
  # Get specific data
  dat_run<- data
  
  # Rescale covariates to SD of 1
  dat_run<- dat_run %>%
    mutate_at(., {{covs_scale}}, vast_scale_func, type = "AJA")
  
  # Model fitting data only
  dat_pred<- dat_run %>%
    filter(., EST_YEAR >=  fore_challenge) %>%
    drop_na({{covs}})
  
  # Predictions
  pa_mod<- mod_list[[1]]
  pa_pred<- as.numeric(predict.gam(pa_mod, newdata = dat_pred, type = "response"))
  
  # Log biomass prediction
  logbio_mod<- mod_list[[2]]
  logbio_pred<- as.numeric(predict.gam(logbio_mod, newdata = dat_pred, type = "response"))
  
  # Create prediction dataframe
  pred_df_out<- dat_pred 
  pred_df_out$PRED_PRESENCE<- pa_pred
  pred_df_out$PRED_LOGBIO<- logbio_pred
  pred_df_out$PRED_BIO<- pred_df_out$PRED_PRESENCE * exp(pred_df_out$PRED_LOGBIO)
  
  # Save and then return it
  saveRDS(pred_df_out, file = paste(outfile, "modelpred.rds", sep = ""))
  return(pred_df_out)
}
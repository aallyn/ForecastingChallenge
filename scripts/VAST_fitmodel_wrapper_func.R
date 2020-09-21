VAST_fitmodel_wrapper<- function(settings = settings_start, strata.limits = strat_limits, Lat_i = samp_dat[,'Lat'], Lon_i = samp_dat[,'Lon'], t_i = as.vector(samp_dat[,'Year']), b_i = samp_dat[,'Catch_KG'], c_iz = rep(0, nrow(samp_dat)), v_i = rep(0, nrow(samp_dat)), Q_ik = NULL, a_i = rep(area_swept, nrow(samp_dat)), PredTF_i = samp_dat[,'PRED_TF'], covariate_data = cov_dat, formula = formula_use, observations_LL = cbind("Lat" = samp_dat[,'Lat'], Lon = samp_dat[, 'Lon']), maximum_distance_from_sample = max_dist_from_sample, grid_dim_km = grid_dim_km, working_dir = outfile){
  
  ## Details
  # This function is a wrapper around the VAST fit_model function that helps users overcome potential issues when trying to fit the VAST model, particularly with spatial, spatio-temporal variability or temporal structure in the intercept or spatio-temporal variability. In a "one-off" situation, it is relatively simple to just add successive fits in the code while trying potential solutions. For example, turning spatio-temporal variability off when the model returns issues with the L_epsilonN_z component. However, this becomes increasingly problematic if we want to run things "in a loop" across multiple species, or even more so, if we want to run things in parallel. Thankfully, Alexa Fredston has already outlined a workflow for addressing these pieces. Here, I am adapting some of her code and putting it into one wrapper function that will make successive tries to fit nested models given potential issues. 
  
  # Args:
  # settings = 
  # strata.limits = 
  # Lat_i = 
  # Lon_i = 
  # t_i = 
  # b_i = 
  # c_iz = 
  # v_i = 
  # Q_ik = 
  # a_i = 
  # PredTF_i = 
  # covariate_data = 
  # formula = 
  # observations_LL = 
  # maximum_distance_from_sample = 
  # grid_dim_km = 
  # working_dir = 
  
  # Returns: 
  
  # Start function ----------------------------------------------------------
  ##########
  ## Preliminaries
  ##########
  # For parallel processing, need to make sure that each core has the necessary libraries
  library(VAST)
  library(TMB)
  
  # For debugging
  if(FALSE){
    settings = settings_start
    strata.limits = strat_limits
    Lat_i = samp_dat[,'Lat']
    Lon_i = samp_dat[,'Lon']
    t_i = as.vector(samp_dat[,'Year'])
    b_i = samp_dat[,'Catch_KG']
    c_iz = rep(0, nrow(samp_dat))
    v_i = rep(0, nrow(samp_dat))
    Q_ik = NULL
    a_i = rep(area_swept, nrow(samp_dat))
    PredTF_i = samp_dat[,'PRED_TF']
    covariate_data = cov_dat
    formula = formula_use
    observations_LL = cbind("Lat" = samp_dat[,'Lat'], Lon = samp_dat[, 'Lon'])
    maximum_distance_from_sample = max_dist_from_sample
    grid_dim_km = grid_dim_km
    working_dir = outfile
  }
  
  ##########
  ## Model fits
  ##########
  
  # Run the model with settings_start, but do NOT calculate SE and covariance pieces. This is just to get some idea for starting parameters that might help us later on. 
  fit_startparams<- try(fit_model("settings" = settings_start, strata.limits = strat.limits, "Lat_i" = Lat_i, "Lon_i" = Lon_i, "t_i" = t_i, "b_i" = b_i, "c_iz" = c_iz, "v_i" = v_i, "Q_ik" = Q_ik, "a_i" = a_i, "PredTF_i" = PredTF_i, covariate_data = covariate_data, formula = formula, "observations_LL" = observations_LL, "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "working_dir" = working_dir, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, "lower" = -Inf, "upper" = Inf, "newtonsteps" = 0, "getsd" = FALSE, "test_fit" = FALSE, "run_model" = FALSE, "Use_REML" = TRUE))
                       
  # Now, run model with same settings, but turning on SE and covariance pieces. If all behaves nicely, this would be the only model we need to run. 
  fit<- try(fit_model("settings" = settings_start, strata.limits = strat.limits, "Lat_i" = Lat_i, "Lon_i" = Lon_i, "t_i" = t_i, "b_i" = b_i, "c_iz" = c_iz, "v_i" = v_i, "Q_ik" = Q_ik, "a_i" = a_i, "PredTF_i" = PredTF_i, covariate_data = covariate_data, formula = formula, "observations_LL" = observations_LL, "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "working_dir" = working_dir, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, "lower" = -Inf, "upper" = Inf, "newtonsteps" = 1, "getsd" = TRUE, "test_fit" = FALSE, "run_model" = TRUE, "Use_REML" = TRUE))
  
  # Check to see if there were any issues. If no issues, add this note to the fitted object
  if(!class(fit)=="try-error"){
    fit$parameter_estimates$adjustments_for_convergence <- "none"
  } 
  
  # If there was an error OR if there were issues with the maximum gradient value, such that they are not below 0.01, we need to try something else. Our first attempt will be supplying the intial parameter estimates from fit_startparams
  # check that all maximum gradients have an absolute value below 0.01 and the model didn't throw an error 
  if(class(fit)=="try-error" || (!class(fit)=="try-error" & max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
    
    # Re-fit with supplied parameter starting values
    fit<- try(fit_model("settings" = settings_start, strata.limits = strat.limits, "Lat_i" = Lat_i, "Lon_i" = Lon_i, "t_i" = t_i, "b_i" = b_i, "c_iz" = c_iz, "v_i" = v_i, "Q_ik" = Q_ik, "a_i" = a_i, "PredTF_i" = PredTF_i, covariate_data = covariate_data, formula = formula, "observations_LL" = observations_LL, "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "working_dir" = working_dir, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, "lower" = -Inf, "upper" = Inf, "newtonsteps" = 1, "getsd" = TRUE, "test_fit" = FALSE, "run_model" = TRUE, "Use_REML" = TRUE, parameters = fit_startparams$ParHat))
    
    # Check fit again, if all okay, then make a note that we used fit_startparams for parameters
    if(!class(fit)=="try-error") {
      fit$parameter_estimates$adjustments_for_convergence <- "used fit0 parameters"
    }
    
    # If that still doesn't work, now we are going to need to start making some targetted adjustments and fitting subsequently simpler, nested models
    if(class(fit)=="try-error" || (!class(fit)=="try-error" & max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
      
      # If no error, but there were convergence issues
      if(!class(fit)=="try-error") {
        if("L_beta1_cf" %in% fit$parameter_estimates$diagnostics$Param & "L_beta2_ct" %in% fit$parameter_estimates$diagnostics$Param){ # AND if these parameters exist in the model... 
          if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001) {
            newRhoConfig = c("Beta1"=3, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
          }
          if(abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
            newRhoConfig = c("Beta1"=4, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
          }
          if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001 & 
             abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
            newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
          }
        }
      }
      
      # If an error and not convergence issues
      if(class(fit)=="try-error"){
        newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
      }
      
      # Refit with new field/rho configuration settings
      fit<- try(fit_model("settings" = settings_start, "FieldConfig" = new_fieldconfig, "RhoConfig" = new_rhoconfig, strata.limits = strat.limits, "Lat_i" = Lat_i, "Lon_i" = Lon_i, "t_i" = t_i, "b_i" = b_i, "c_iz" = c_iz, "v_i" = v_i, "Q_ik" = Q_ik, "a_i" = a_i, "PredTF_i" = PredTF_i, covariate_data = covariate_data, formula = formula, "observations_LL" = observations_LL, "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "working_dir" = working_dir, "getReportCovariance" = FALSE, "getJointPrecision" = TRUE, "lower" = -Inf, "upper" = Inf, "newtonsteps" = 1, "getsd" = TRUE, "test_fit" = FALSE, "run_model" = TRUE, "Use_REML" = TRUE, parameters = fit_startparams$ParHat))
                  
      if(!class(fit)=="try-error"){
        fit$parameter_estimates$adjustments_for_convergence <- "changed field or rho config"
      }
    }
  }
  
  
  # check that all maximum gradients have an absolute value below 0.01 and the model didn't throw an error 
  
  if(class(fit)=="try-error" || (!class(fit)=="try-error" &  max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
    
    # attempt 3: if fit still failed / didn't converge, it could be because L_beta1_cf or L_beta2_ct is approaching zero; use a different RhoConfig 
    
    # if it's a model, but didn't converge... 
    if(!class(fit)=="try-error") {
      if("L_beta1_cf" %in% fit$parameter_estimates$diagnostics$Param & "L_beta2_ct" %in% fit$parameter_estimates$diagnostics$Param){ # AND if these parameters exist in the model... 
        if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001) {
          newRhoConfig = c("Beta1"=3, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
        }
        if(abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
          newRhoConfig = c("Beta1"=4, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
        }
        if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001 & 
           abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
          newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
        }}
    }
    
    # if it's an error...
    if(class(fit)=="try-error"){
      newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
    }
    
    fit = try(
      fit_model("settings"=settings, 
                "RhoConfig"= newRhoConfig,
                "Lat_i"=Data_Geostat[,'Lat'],
                "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                lower=-Inf, upper=Inf,
                test_fit = FALSE,
                fine_scale=fine_scale,
                anisotropy=FALSE,
                Use_REML=TRUE,
                Z_gm = Z_gm,
                getsd=TRUE,
                newtonsteps=1
                
      ))
    if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "changed RhoConfig"}
    
  }
}

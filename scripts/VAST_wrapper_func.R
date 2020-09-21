#####
## VAST wrapper functions to apply to each species
#####
vast_scale_func<- function(x, type){
  if(type == "JT"){
    x.out<- x/100
    return(x.out)
  } else {
    x.out<- as.numeric(scale(abs(x)))
    return(x.out)
  }
}

VAST_data_func<- function(run.name, all.data, train.test.year, species, survey.season, covariates, output.dir){
  # Details -----------------------------------------------------------------
  # This function processes a fisheries data set into two components, a sample dataset and a covariate dataset, which are used by the VAST_fit_func wrapper

  # Arguments:
  # run.name = Character string. A unique name to help you identify the output data. For example: "NEFSCFall_WinterFlounder."
  # all.data = Fisheries survey data set, either path to an RDS file or an existing data frame. At a minimum, this data set needs to include columns for species, Year, lat, lon, Catch_KG. Additional columns should be for covariates included in the VAST model.
  # train.test.year = Numeric. The cutoff year to split the data into training and testing subsets. Training data will be created as all samples taken in years less than or equal to the train.test.year.
  # species = Character vector. The species to model -- species names should be specified in sentence case (e.g., "Atlantic Cod"). If "All" all groundfish are included.
  # survey.season = Character vector. The survey and seasons to include in the model.
  # covariates = Character vector. The covariates to include in the output covariate dataset. In the output covariate dataset, these are scaled and centered.
  # output.dir = Character string. The path to store the sample data and covariate data.

  # Debugging ---------------------------------------------------------------
  if(FALSE){
    run.name = run.name = paste("NEFSC", season, gsub(" ", "", paste("_", to_any_case(spp, case = "sentence"), sep = "")), sep = "")
    all.data = dat.sub
    train.test.year = 1999
    species = spp
    survey.season = survey.season
    covariates = coeff.get.covs
    output.dir = OutFile
  }

  # Read in the data
  dat<- if(class(all.data) == "character"){
    readRDS(all.data) } else {
      all.data
    }

  # Filtering based on train.test.year, species, and survey.season
  species.mod<- if(species == "All"){
    unique(dat$commonname)
  } else {
    species
  }

  dat.mod<- dat %>%
    dplyr::filter(., commonname %in% species.mod & surveyfact %in% survey.season)

  # We won't be able to have NA values in any of the locations OR in any of the covariates that we want to model.
  dat.mod<- dat.mod %>%
    drop_na(lat, lon, paste0(covariates))

  # Some data processing work
  dat.mod$Season<- ifelse(dat.mod$surveyfact == "NEFSC_NEUSSpring", "Spring", "Fall")

  #dat.mod$Season.Num<- ifelse(dat.mod$Season == "Spring", 0, 1)
  # dat.mod<- dat.mod %>%
  #   arrange(., year, Season.Num)
  # dat.mod$Season<- factor(dat.mod$Season, levels = c("Spring", "Fall"))

  dat.mod$Season.Num<- rep(1, nrow(dat.mod))
  dat.mod<- dat.mod %>%
    arrange(., year, Season.Num)
  dat.mod$Season<- factor(dat.mod$Season)
  dat.mod$Year.Fac<- factor(dat.mod$year, levels = unique(dat.mod$year))

  # Finally, an ordered vector, with year plus decimal for each of the seasons (or months) within a year.
  dat.mod<- dat.mod %>%
    group_by(year) %>%
    mutate(., "Yearly.Season.Levels" = length(unique(Season)),
           "Ordered.Time" = round(year + (Season.Num/Yearly.Season.Levels), digits = 1)) %>%
    arrange(Ordered.Time)
  date.vec<- unique(paste(dat.mod$year, dat.mod$Season, sep = "-"))
  dat.mod$Date<- factor(paste(dat.mod$year, dat.mod$Season, sep = "-"), levels = date.vec[order(dat.mod$Ordered.Time)])

  # Organizing response/sample data
  samp.dat<- data.frame("spp" = dat.mod$commonname, "Year" = dat.mod$year, "Year.Fac" = dat.mod$Year.Fac, "Season" = dat.mod$Season, "Date" = dat.mod$Date, "Ordered.Time" = dat.mod$Ordered.Time, "Survey" = dat.mod$surveyfact, "Strata" = dat.mod$STRATA, "Lat" = dat.mod$lat, "Lon" = dat.mod$lon, "Catch_KG" = dat.mod$wtcpue, "AreaSwept_km2" = rep(0.01, nrow(dat.mod)))

  # Add predictionTF column?
  if(!is.null(train.test.year)){
    samp.dat$Pred_TF<- ifelse(samp.dat$Year <= train.test.year, 0, 1)
  }

  # Organizing the covariate data
  cov.dat.base<- data.frame("Year" = dat.mod$year, "Date" = dat.mod$Date, "Lat" = dat.mod$lat, "Lon" = dat.mod$lon)
  cov.dat<- data.frame(cov.dat.base, dat.mod[,which(colnames(dat.mod) %in% covariates)])

  # VAST prefers scaled covariates
  cov.dat.scaled<- cov.dat %>%
    mutate_at(covariates, vast_scale_func, "Scale")

  # Need to keep mean and sd from rescale to use when we predict or project to other time periods
  rescale.df<- cov.dat %>%
    mutate_at(covariates, .funs = c(mean, sd))

  # Save these and return them
  out.list<- list("SampleData" = samp.dat, "CovariateData" = cov.dat.scaled, "RescaleData" = rescale.df)
  saveRDS(out.list, paste(output.dir, "/data.rds", sep = ""))
  return(out.list)
}

VAST_fit_func<- function(sample.data = out.list$SampleData,
                         covariate.data = out.list$CovariateData,
                         region.use = "northwest_atlantic", utm.zone.use = NA, strata.limits.use = data.frame('STRATA' = "All_areas"), knots.use = 100, spde.method.use = "Mesh",
                         obs.model.use = cbind("PostDist" = 1, "Link" = 1),
                         field.config.use = c("Omega1" = length(unique(samp.dat$spp)), "Epsilon1" = length(unique(samp.dat$spp)), "Omega2" = length(unique(samp.dat$spp)), "Epsilon2" = length(unique(samp.dat$spp))), rho.config.use = c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 4, "Epsilon2" = 4), aniso.use = TRUE,
                         mod.formula = ~ Depth.scale + I(Depth.scale^2) + SBT.seasonal.scale + I(SBT.seasonal.scale^2),
                         vessels.factors.use, overdispersion.config.use = c("Eta1" = 0, "Eta2" = 0),
                         derived.quants = c('SD_site_logdensity' = FALSE, 'Calculate_Range' = FALSE, 'Calculate_effective_area' = FALSE, 'Calculate_Cov_SE' = FALSE, 'Calculate_Synchrony' = FALSE, 'Calculate_proportion'=FALSE), bias.correct.use = TRUE,
                         output.dir.use){

  # Details -----------------------------------------------------------------
  # This is a wrapper function around Jim Thorson's VAST model to hopefully facilitate fitting running the model across multiple species/time ranges.

  # Arguments:
    # sample.data = Data frame of the response data, with minimum columns for spp, Year, lat, lon, Catch_KG and AreaSwept_km2
    # covariate.data =  Data frame of the covariates data associated with each sample with mimumum columns for Year, lat, lon. All covariates should be scaled/centered.
    # region.use = "northwest_atlantic" (default). Character vector for help with detemrining extrapolation grid. See "FishStatsUtils::make_extrapolation_grid" for more details.
    # utm.zone.use = NA (default). UTM zone for projecting lat/long to km distances. Zone = NA will automatically detect utm zone based on location of extrapolation grid samples.
    # strata.limits.use = data.frame('STRATA' = "All_areas") (default). Input to determine startifications of indices. Not clear how this works when region = "Other."
    # knots.use = Numeric. Specify number of stations (a.k.a. "knots"). This defines the grain of the SPDE mesh when method = "Mesh"
    # spde.method.use = "Mesh" (default). The approximation method to use when estimating spatial and spatio-temporal variation, which uses stochastic partial differential equation approximation.
    # obs.model.use = cbind("PostDist" = 1, "Link" = 1). Specifies response distribution and link function. See "?VAST::make_data" for more options.
    # field.config.use = c("Omega1" = length(unique(samp.dat$spp)), "Epsilon1" = length(unique(samp.dat$spp)), "Omega2" = length(unique(samp.dat$spp)), "Epsilon2" = length(unique(samp.dat$spp))) (default). Spatial and spatio-temporal factors, where omega = spatial variation and epsilon = spatio-temporal variation. Numeric value signals how many factors to include in factor analysis covariance.
    # rho.config.use = c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 4, "Epsilon2" = 4) (default). Temporal correlation of intercepts (Beta) or spatio-temporal variation (Epsilon). See "?VAST::make_data" for other options.
    # aniso.use = TRUE (default). Whether to assume isotropy (aniso.use = FALSE), or geometric anisotropy (aniso.use = TRUE)
    # mod.formula = Object of class formula describing relationship between response and covariates.
    # vessels.factors.use =
    # overdispersion.config.use = overdispersion.config.use = c("Eta1" = 0, "Eta2" = 0) (default). Determining if there is correlated overdispersion among vessels.
    # derived.quants = Built in options for calculating values from fitted models. See "?VAST::make_data" for options
    # bias.correct.use = TRUE (default). Whether or not to bias correct derived quantity values for systematic differences that arise when non-linearly transforming response and predicted values.
    # output.dir.use = Output directory to store model settings and fitted model objects.


  # Debugging ---------------------------------------------------------------
  if(FALSE){
    sample.data = samp.dat
    covariate.data = cov.dat
    region.use = "northwest_atlantic"
    utm.zone.use = NA
    strata.limits.use = data.frame('STRATA' = "All_areas")
    knots.use = 100
    spde.method.use = "Mesh"
    obs.model.use = cbind("PostDist" = 1, "Link" = 1)
    field.config.use = c("Omega1" = length(unique(samp.dat$spp)), "Epsilon1" = length(unique(samp.dat$spp)), "Omega2" = length(unique(samp.dat$spp)), "Epsilon2" = length(unique(samp.dat$spp)))
    aniso.use = TRUE
    rho.config.use = c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 4, "Epsilon2" = 4)
    mod.formula = ~ Depth.scale + I(Depth.scale^2) + SBT.seasonal.scale + I(SBT.seasonal.scale^2)
    vessels.factors.use = 0
    overdispersion.config.use = c("Eta1" = 0, "Eta2" = 0)
    derived.quants = c('SD_site_logdensity' = FALSE, 'Calculate_Range' = FALSE, 'Calculate_effective_area' = FALSE, 'Calculate_Cov_SE' = FALSE, 'Calculate_Synchrony' = FALSE, 'Calculate_proportion'=FALSE)
    bias.correct.use = TRUE
    output.dir.use = "~/GitHub/SDM-convergence/VAST_WrapperTest/"
  }


  # Compiling data and VAST components --------------------------------------
  # Spatial extrapolation grid
  extrap.dat<- make_extrapolation_info(Region = region.use, zone = utm.zone.use, strata.limits = strata.limits.use, observations_LL = cbind("Lat" = sample.data$Lat, "Lon" = sample.data$Lon))
  extrap.dat$a_el$STRATA<- extrap.dat$Area_km2_x

  # Spatial and spatio-temporal factor information
  spat.dat<- make_spatial_info(n_x = knots.use, Lon_i = sample.data$Lon, Lat_i = sample.data$Lat, knot_method = "samples", Method = spde.method.use, fine_scale = TRUE, Extrapolation_List = extrap.dat)

  # We need to add the knot info to our sampling_data dataset...maybe?
  sample.data<- cbind(sample.data, knot_i = spat.dat$knot_i)

  # Now, making the data
  b.i.length<- length(sample.data$Catch_KG)
  vast.dat<- make_data(b_i = sample.data$Catch_KG,
                       a_i = sample.data$AreaSwept_km2,
                       t_iz = cbind(sample.data$Year),
                       c_iz = as.numeric(sample.data$spp)-1,
                       v_i = vessels.factors.use,
                       spatial_list = spat.dat,
                       FieldConfig = field.config.use,
                       RhoConfig = rho.config.use,
                       ObsModel_ez = obs.model.use,
                       OverdispersionConfig = overdispersion.config.use,
                       VamConfig = c(Method = 0, Rank = 0, Timing = 0),
                       Aniso = aniso.use,
                       PredTF_i = rep(0, b.i.length),
                       Xconfig_zcp = NULL,
                       covariate_data = covariate.data,
                       formula = mod.formula,
                       Q_ik = NULL,
                       F_ct = NULL,
                       t_yz = NULL,
                       CheckForErrors = TRUE,
                       yearbounds_zz = NULL,
                       Options = derived.quants,
                       Expansion_cz = NULL)

  # Saving
  OutFile<- output.dir.use
  if(!file.exists(OutFile)){
    dir.create(OutFile)
  }

  Version<- get_latest_version(package="VAST")
  Record<- list("Version" = Version,
                "Method" = spde.method.use,
                "grid_size_km" = NA,
                "n_x" = knots.use,
                "FieldConfig" = field.config.use,
                "RhoConfig" = rho.config.use,
                "OverdispersionConfig" = overdispersion.config.use,
                "ObsModel" = obs.model.use,
                "Region" = region.use,
                "Species_set" = spp.set.use,
                "strata.limits" = strata.limits.use)
  save(Record, file = file.path(OutFile, "Record.RData"))
  capture.output(Record, file = paste0(OutFile, "Record.txt"))

  # Build and fit VAST model ------------------------------------------------
  TmbList<- make_model("TmbData" = vast.dat,
                       "RunDir" = OutFile,
                       "Version" = Version,
                       "RhoConfig" = rho.config.use,
                       "loc_x" = spat.dat$loc_x,
                       "Method" = spat.dat$Method)
  Obj<- TmbList[["Obj"]]

  TMBfit<- TMBhelper::fit_tmb(obj = Obj,
                              lower = TmbList[["Lower"]],
                              upper = TmbList[["Upper"]],
                              getsd = TRUE,
                              savedir = OutFile,
                              bias.correct = bias.correct.use,
                              newtonsteps = 1)
  Report<- Obj$report()

  # Diagnostics -------------------------------------------------------------
  # Knot locations and sample locations
  plot_data(Extrapolation_List = extrap.dat, Spatial_List = spat.dat, Data_Geostat = sample.data, PlotDir = OutFile)

  # Convergence
  cov.check.table<- (TMBfit$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')])
  write.csv(cov.check.table, file = paste(OutFile, "ConvergenceParamaterCheck.csv", sep = ""))

  # Encounter probability diagnostics
  enc.prob<- plot_encounter_diagnostic(Report = Report, Data_Geostat = sample.data, DirName = OutFile)

  # Positive catch rate diagnostics
  Q<- plot_quantile_diagnostic(TmbData = vast.dat, Report = Report, FileName_PP = "Posterior_Predictive", FileName_Phist = "Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile = OutFile)

  # Residuals
  # Get region-specific settings for plots
  MapDetails_List<- make_map_info("Region" = region.use, "spatial_list" = spat.dat, "Extrapolation_List" = extrap.dat)

  # Decide which years to plot
  Year_Set<- seq(min(sample.data[,'Year']), max(sample.data[,'Year']))
  Years2Include<- which(Year_Set %in% sort(unique(sample.data[,'Year'])))

  plot_residuals(Lat_i = sample.data[,'Lat'], Lon_i = sample.data[,'Lon'], TmbData = vast.dat, Report = Report, Q = Q, working_dir = OutFile, spatial_list = spat.dat, extrapolation_list = extrap.dat, Year_Set = Year_Set, Years2Include = Years2Include, mar = c(0,0,2,0), oma = c(3.5,3.5,0,0), cex = 1.8)

  # Inferences --------------------------------------------------------------
  plot_anisotropy(FileName = paste0(OutFile, "Aniso.png"), Report = Report, TmbData = vast.dat)

  #Cov_List<- summarize_covariance(Report = Report, ParHat = Obj$env$parList(), Data = vast.dat, SD = TMBfit$SD, plot_cor = TRUE, category_names = levels(sample.data[,'spp']), plotdir = OutFile, plotTF = field.config.use, mgp = c(2,0.5,0), tck = -0.02, oma = c(0,5,2,2))

  plot_maps(plot_set = c(3), Report = Report, Sdreport = TMBfit$SD, category_names = levels(sample.data[,'spp']), working_dir = OutFile, PlotDF = MapDetails_List[["PlotDF"]])

  Index<- plot_biomass_index(DirName = OutFile, TmbData = vast.dat, Sdreport = TMBfit[["SD"]], Year_Set = Year_Set, Years2Include = Years2Include, strata_names = strata.limits.use[,1], use_biascorr = bias.correct.use, category_names = levels(sample.data[,'spp']))
  if(length(levels(sample.data[,'spp'])) == 1){
    Index.table.out<- Index$Table[,c("Year", "Estimate_metric_tons", "SD_mt")]
    Index.table.out$Category<- rep(levels(sample.data[,'spp']), nrow(Index.table.out))
  } else {
    Index.table.out<- Index$Table[,c("Category", "Year", "Estimate_metric_tons", "SD_mt")]
  }

  write.csv(Index.table.out, file = paste(OutFile, "YearlyAbundanceEstimates.csv", sep = ""))

  if(derived.quants['Calculate_effective_area']){
    plot_range_index(Report = Report, TmbData = vast.dat, Sdreport = TMBfit[["SD"]], Znames = colnames(vast.dat$Z_xm), PlotDir = OutFile, category_names = levels(sample.data[,'spp']), Year_Set = Year_Set)
  }

  if(any(overdispersion.config.use != 0)){
    plot_overdispersion(filename1 = paste0(OutFile,"Overdispersion"), filename2 = paste0(OutFile, "Overdispersion--panel"), Data = vast.dat, ParHat = ParHat, Report = Report, ControlList1 = list("Width" = 5, "Height" = 10, "Res" = 200, "Units" = 'in'), ControlList2 = list("Width" = vast.dat$n_c, "Height" = vast.dat$n_c, "Res" = 200, "Units" = 'in'))
  }

  #plot_factors(Report = Report, ParHat = TMBfit$env$parList(), Data = vast.dat, SD = TMBfit$SD, mapdetails_list = MapDetails_List, Year_Set = Year_Set, category_names = levels(sample.data[,'spp']), plotdir = OutFile)

}

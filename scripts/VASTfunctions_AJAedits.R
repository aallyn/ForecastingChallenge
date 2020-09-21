#####
## Andrew's VAST function edits
#####


# Overview ----------------------------------------------------------------
# I ran into a bit of trouble fitting the VAST model and visualizing the results using the packaged functions. To overcome these issues, I had to make a few edits to a variety of Jim's functions. Those edits are here. I also include a simple function to calculate the number of covariates and then a "recaling" function for rescaling covariate values to meet the SD < 10 for stability requirement.
# Function for calculating number of coefficients from model formula
n_coeffs_func<- function(mod.formula, cov.data){
  
  # First get the coefficient names from the model formula
  coef.names<- attributes(terms(mod.formula))$term.labels
  
  # Now get class of each of these
  cov.classes<- sapply(cov.data, "class")
  
  n.coeffs<- 0
  
  for(i in seq_along(coef.names)){
    # Need to account for the quadratic terms as these are not going to be in the data frame...
    if(str_detect(coef.names[i], "\\^")){
      coef.adj<- str_replace_all(coef.names[i], c("I" = "", "\\(" = "", "\\^2\\)" = ""))
      coef.use<- cov.classes[names(cov.classes) == coef.adj]
    } else {
      coef.use<- cov.classes[names(cov.classes) == coef.names[i]]
    }
    if(coef.use == "factor"){
      n.coeffs.add<- length(levels(cov.data[,names(coef.use)]))
      n.coeffs<- n.coeffs + n.coeffs.add
    } else {
      n.coeffs.add<- 1
      n.coeffs<- n.coeffs + n.coeffs.add
    }
  }
  return(n.coeffs)
}

# Rescaling covariates
vast_scale_func<- function(x, type){
  if(type == "JT"){
    x.out<- x/100
    return(x.out)
  } else {
    x.out<- as.numeric(scale(abs(x)))
    return(x.out)
  }
}

# VAST data building functions (eventually called by “fit_model_zz") -------
make_covariates_zz<- function(formula = ~0, covariate_data, Time_i = t_iz, spatial_list, extrapolation_list, design = design, ...) {
  # For debugging, walk through "make_data" first..
  if(FALSE){
    formula = formula
    covariate_data = covariate_data
    Time_i = t_iz
    spatial_list = spatial_list
    extrapolation_list = extrapolation_list
    design = design
  }
  if (!is.data.frame(covariate_data))
    stop("Please ensure that `covariate_data` is a data frame")
  if (!all(c("Lat", "Lon", "Year") %in% names(covariate_data))) {
    stop("`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`")
  }

  # Formatting sample data and getting a sequence of dates to match to knots/extrapolation grid cells? Need to have a similar column in covariate data, too...
  if(ncol(Time_i) == 2){
    sample_data = data.frame(Year = Time_i$Year, Date = Time_i$Date, Lat = spatial_list$latlon_i[, "Lat"], Lon = spatial_list$latlon_i[, "Lon"])
    Date_Set = unique(sample_data$Date) # Andrew edit
  } else if(ncol(Time_i) == 3 & design == "Season"){
    sample_data = data.frame(Year = Time_i$Year, Season = Time_i$Season, Date = Time_i$Date, Lat = spatial_list$latlon_i[, "Lat"], Lon = spatial_list$latlon_i[, "Lon"])
    Date_Set = unique(sample_data$Date) # Andrew edit
  }

  covariate_names_all = attr(terms.formula(formula), "term.labels")
  covariate_names_adj = attr(terms.formula(formula), "term.labels")

  # Some modification for quadratic terms...
  adjust<- str_which(covariate_names_all, "\\^")

  if(length(adjust)>0){
    coef.adj<- str_replace_all(covariate_names_adj, c("I" = "", "\\(" = "", "\\^2\\)" = ""))
    covariate_names_adj<- covariate_names_all[covariate_names_adj == coef.adj]
  }

  # Which do we need to add?
  covariate_names_to_add = setdiff(covariate_names_adj, names(sample_data))
  latlon_g = spatial_list$latlon_g

  # Creating a dataframe with same rows as observed points, but now all covariates are set to NA?
  DF_zp = NULL

  # Do this based on covariate names
  if(is_empty(covariate_names_to_add) == TRUE){
    DF_ip = data.frame(sample_data)
  } else {
    DF_ip = data.frame(sample_data, covariate_data[rep(1, nrow(sample_data)), covariate_names_to_add])
    names.new<- ncol(DF_ip) - ncol(sample_data)
    names(DF_ip)[(ncol(sample_data)+1):(ncol(sample_data)+names.new)]<- covariate_names_to_add
  }

  # Set covariates to NA
  DF_ip[, covariate_names_to_add] = NA

  # Adjust covariate column types...
  cov.classes<- sapply(covariate_data, "class")
  cov.factors<- intersect(attributes(cov.classes[which(cov.classes == "factor")])$names, covariate_names_adj)

  for(covI in seq_along(cov.factors)){
    DF.col.ind<- which(colnames(DF_ip) == cov.factors[covI])
    samp.dat.ind<- which(colnames(covariate_data) == cov.factors[covI])
    DF_ip[,DF.col.ind]<- factor(DF_ip[,DF.col.ind], levels = levels(covariate_data[,samp.dat.ind]))
  }

  # Loop over each date...
  for (tI in seq_along(Date_Set)) {
    tmp_covariate_data = covariate_data[which(Date_Set[tI] == covariate_data[, "Date"] | is.na(covariate_data[, "Date"])), , drop = FALSE]
    if (nrow(tmp_covariate_data) == 0) {
      stop("Date ", Date_Set[tI], " not found in `covariate_data` please specify covariate values for all dates")
    }
    Which = which(Date_Set[tI] == sample_data[, "Date"])
    if (length(Which) > 0) {
      NN = RANN::nn2(data = tmp_covariate_data[, c("Lat", "Lon")], query = sample_data[Which, c("Lat", "Lon")], k = 1)
      nearest_covariates = tmp_covariate_data[NN$nn.idx[, 1], covariate_names_adj, drop = FALSE]
      DF_ip[Which, covariate_names_adj] = nearest_covariates
    }
    NN = RANN::nn2(data = tmp_covariate_data[, c("Lat", "Lon")], query = latlon_g[, c("Lat", "Lon")], k = 1)
    nearest_covariates = tmp_covariate_data[NN$nn.idx[, 1], covariate_names_adj, drop = FALSE]
    newrows = data.frame(Date = rep(as.character(unique(tmp_covariate_data$Date)), nrow(latlon_g)), latlon_g, nearest_covariates)
    DF_zp = rbind(DF_zp, newrows)
  }
  if (any(is.na(DF_ip)))
    stop("Problem with `DF_ip` in `make_covariates(.)")
  DF_ip<- DF_ip[, intersect(names(DF_zp), names(DF_ip))]
  DF = rbind(DF_ip, DF_zp)

  # This was dropping one of the year levels....Do I want that??? I don't think so, seems like this still needs to have the baseline year....
  X = model.matrix(update.formula(formula, ~. + 0), data = DF, contrasts.arg = lapply(DF[,cov.factors], contrasts, contrasts=FALSE))[, , drop = FALSE]
  X_ip = X[1:nrow(DF_ip), , drop = FALSE]
  X_itp = aperm(X_ip %o% rep(1, length(Date_Set)), perm = c(1, 3, 2))
  X_gpt = NULL
  indices = nrow(X_ip)
  for (tI in seq_along(Date_Set)) {
    indices = max(indices) + 1:nrow(latlon_g)
    if (max(indices) > nrow(X))
      stop("Check problem in `make_covariates`")
    X_gpt = abind::abind(X_gpt, X[indices, , drop = FALSE], along = 3)
  }
  X_gtp = aperm(X_gpt, perm = c(1, 3, 2))
  if (any(apply(X_gtp, MARGIN = 2:3, FUN = sd) > 10 | apply(X_itp, MARGIN = 2:3, FUN = sd) > 10)) {
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }
  Return = list(X_gtp = X_gtp, X_itp = X_itp, covariate_names = covariate_names_all)
  return(Return)
}

make_data_zz<- function(b_i, a_i, t_iz, c_iz = rep(0, length(b_i)), e_i = c_iz[,1], v_i = rep(0, length(b_i)), FieldConfig, spatial_list, ObsModel_ez = c(PosDist = 1, Link = 0), OverdispersionConfig = c(eta1 = 0, eta2 = 0), RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0), VamConfig = c(Method = 0, Rank = 0, Timing = 0), Aniso = TRUE, PredTF_i = rep(0, length(b_i)), Xconfig_zcp = NULL, covariate_data = NULL, formula = ~ 0, design = design, Q_ik = NULL, Network_sz = NULL, F_ct = NULL, F_init = 1, t_yz = NULL, CheckForErrors = TRUE, yearbounds_zz = NULL, Options = c(), Expansion_cz = NULL, Z_gm = NULL, Version = FishStatsUtils::get_latest_version(package = "VAST"), ...) {

  if(FALSE){
    Version = data_args_input$Version
    FieldConfig = data_args_input$FieldConfig
    OverdispersionConfig = data_args_input$OverdispersionConfig
    RhoConfig = data_args_input$RhoConfig
    VamConfig = data_args_input$VamConfig
    ObsModel = data_args_input$ObsModel
    ObsModel_ez = ObsModel
    design = data_args_input$design
    c_iz = data_args_input$c_iz
    e_i = c_iz
    b_i = data_args_input$b_i
    a_i = data_args_input$a_i
    v_i = data_args_input$v_i
    t_iz = data_args_input$t_iz
    Options = data_args_input$Options
    Aniso = data_args_input$Aniso
    Xconfig_zcp = data_args_input$Xconfig_zcp
    covariate_data = data_args_input$covariate_data
    formula = data_args_input$formula
    Q_ik = data_args_input$Q_ik
    spatial_list = data_args_input$spatial_list
  }

  deprecated_inputs = list(...)
  X_xtp = deprecated_inputs[["X_xtp"]]
  X_gtp = deprecated_inputs[["X_gtp"]]
  X_itp = deprecated_inputs[["X_itp"]]
  X_xj = deprecated_inputs[["X_xj"]]

  if (missing(spatial_list)) {
    warning("Consider changing use of `make_data` to include `spatial_list` as input")
    a_xl = a_gl = deprecated_inputs[["a_xl"]]
    MeshList = deprecated_inputs[["MeshList"]]
    GridList = deprecated_inputs[["GridList"]]
    Method = deprecated_inputs[["Method"]]
    s_i = deprecated_inputs[["s_i"]]
  } else {
    MeshList = spatial_list[["MeshList"]]
    GridList = spatial_list[["GridList"]]
    Method = spatial_list[["Method"]]
    a_xl = a_gl = spatial_list[["a_gl"]]
    s_i = spatial_list[["knot_i"]] - 1
  }
  Options2use = c(SD_site_density = FALSE, SD_site_logdensity = FALSE,
                  Calculate_Range = FALSE, SD_observation_density = FALSE,
                  Calculate_effective_area = FALSE, Calculate_Cov_SE = FALSE,
                  Calculate_Synchrony = FALSE, Calculate_Coherence = FALSE,
                  Calculate_proportion = FALSE, normalize_GMRF_in_CPP = TRUE,
                  Calculate_Fratio = FALSE, Estimate_B0 = FALSE, Project_factors = FALSE, treat_nonencounter_as_zero = FALSE, simulate_random_effects = TRUE)
  for (i in seq_along(Options)) {
    if (tolower(names(Options)[i]) %in% tolower(names(Options2use))) {
      Options2use[[match(tolower(names(Options)[i]), tolower(names(Options2use)))]] = Options[[i]]
    }
  }
  if (FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0")) {
    if (!is.null(X_xtp)) {
      stop("`X_xtp` is not used in version >= 8.0.0. If you'd like to specify covariates using input `X_xtp` please use `Version='VAST_v7_0_0'`")
    }
  }
  if (FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v7_0_0")) {
    if (!is.null(X_gtp) | !is.null(X_itp)) {
      stop("`X_gtp` and `X_itp` are not used in version <= 7.0.0. If you'd like to specify covariates using input `X_gtp` and `X_itp` please use `Version='VAST_v8_0_0'` or higher")
    }
  }
  if(is.vector(FieldConfig) && length(FieldConfig) == 4) {
    # ORDER MATTERS!
    order.use<- c("Omega1", "Omega2", "Epsilon1", "Epsilon2")
    FieldConfig<- FieldConfig[order.use]
    FieldConfig = rbind(matrix(FieldConfig, nrow = 2, ncol = 2, dimnames = list(c("Omega", "Epsilon"), c("Component_1", "Component_2")), byrow = TRUE), Beta = c(Beta1 = "IID", Beta2 = "IID"))
  } else {
    if (!is.matrix(FieldConfig) || !all(dim(FieldConfig) == c(3, 2))) {
      stop("`FieldConfig` has the wrong dimensions in `make_data`")
    } else {
      dimnames(FieldConfig) = list(c("Omega", "Epsilon", "Beta"), c("Component_1", "Component_2"))
    }
  }

  # Time vector
  tprime_iz<- as.numeric(t_iz$Date)-1
  tprime_iz = matrix(tprime_iz, ncol = 1)

  if (Options2use[12] == 1) {
    tprime_iz = tprime_iz + 1
    F_ct = cbind(0, F_ct)
  }
  if (!is.matrix(c_iz))
    c_iz = matrix(c_iz, ncol = 1)
  n_t = max(tprime_iz, na.rm = TRUE) + 1
  n_c = max(c_iz, na.rm = TRUE) + 1
  n_e = max(e_i) + 1
  n_v = length(unique(v_i))
  n_i = length(b_i)
  n_x = nrow(a_gl)
  n_l = ncol(a_gl)
  n_g = ifelse(is.null(spatial_list), 1, spatial_list$n_g)
  if (!is.matrix(ObsModel_ez))
    ObsModel_ez = matrix(ObsModel_ez, ncol = 2, nrow = n_e, byrow = TRUE)
  if (Options2use["treat_nonencounter_as_zero"] == TRUE) {
    Index = list(factor(c_iz[, 1], levels = 0:max(c_iz[, 1])), factor(tprime_iz[, 1], levels = 0:max(tprime_iz[, 1])))
    Num_ct = tapply(b_i, INDEX = Index, FUN = function(vec) {
      sum(vec > 0)
    })
    Num_ct = ifelse(is.na(Num_ct), 0, Num_ct)
    b_i = ifelse(Num_ct[cbind(as.numeric(Index[[1]]), as.numeric(Index[[2]]))] == 0, NA, b_i)
  }
  if (is.null(X_xj)) {
    X_xj = matrix(0, nrow = n_x, ncol = 1)
  } else {
    if (!is.array(X_xj) || !(dim(X_xj)[1] == n_x)) {
      stop("`X_xj` has wrong dimensions")
    }
  }
  if (is.null(X_xtp)) {
    X_xtp = array(0, dim = c(n_x, n_t, 1))
  } else {
    if (!is.array(X_xtp) || !(all(dim(X_xtp)[1:2] == c(n_x, n_t)))) {
      stop("`X_xtp` has wrong dimensions")
    }
  }
  Works = FALSE
  if (is.null(X_gtp) & is.null(X_itp)) {
    if (is.null(covariate_data)) {
      X_gtp = array(0, dim = c(n_g, n_t, 1))
      X_itp = array(0, dim = c(n_i, n_t, 1))
      Works = TRUE
    }
    if (!is.null(covariate_data)) {
      covariate_list = make_covariates_zz(formula = formula, covariate_data = covariate_data, Time_i = t_iz, spatial_list = spatial_list, extrapolation_list = extrapolation_list, design = design)
      X_gtp = covariate_list$X_gtp
      X_itp = covariate_list$X_itp
      if (dim(X_gtp)[3] == 0)
        X_gtp = array(0, dim = c(n_g, n_t, 1))
      if (dim(X_itp)[3] == 0)
        X_itp = array(0, dim = c(n_i, n_t, 1))
      Works = TRUE
    }
  }
  if (!is.null(X_gtp) & !is.null(X_itp)) {
    if (!is.array(X_gtp) || !(all(dim(X_gtp)[1:2] == c(n_g, n_t)))) {
      stop("`X_gtp` has wrong dimensions")
    }
    if (!is.array(X_itp) || !(all(dim(X_itp)[1:2] == c(n_i, n_t)))) {
      stop("`X_itp` has wrong dimensions")
    }
    Works = TRUE
  }
  if (Works == FALSE)
    stop("Report problem with covariate asembly to package developer")
  if (is.null(Q_ik)) {
    Q_ik = matrix(0, nrow = n_i, ncol = 1)
  } else {
    if (!is.array(Q_ik) || !(all(dim(Q_ik)[1] == c(n_i)))) {
      stop("`Q_ik` has wrong dimensions")
    }
  }
  if (is.null(yearbounds_zz)) {
    yearbounds_zz = matrix(c(0, n_t - 1), nrow = 1)
  } else {
    if (!is.array(yearbounds_zz) || !(all(dim(yearbounds_zz)[2] == 2))) {
      stop("`yearbounds_zz` has wrong dimensions")
    }
  }
  if (is.null(t_yz)) {
    t_yz = matrix(0:max(tprime_iz[, 1], na.rm = TRUE), ncol = 1)
    for (cI in seq(2, ncol(tprime_iz), length = ncol(tprime_iz) - 1)) t_yz = cbind(t_yz, min(tprime_iz[, cI], na.rm = TRUE))
  }
  n_j = ncol(X_xj)
  if (FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0")) {
    n_p = dim(X_gtp)[3]
  } else {
    n_p = dim(X_xtp)[3]
  }
  n_k = ncol(Q_ik)
  n_y = nrow(t_yz)
  if (is.null(F_ct)) {
    F_ct = matrix(0, nrow = n_c, ncol = n_t)
  } else {
    if (!is.array(F_ct) || !(all(dim(F_ct) == c(n_c, n_t)))) {
      stop("`F_ct` has wrong dimensions")
    }
  }
  if (is.null(Expansion_cz)) {
    Expansion_cz = matrix(0, nrow = n_c, ncol = 2)
  } else {
    if (!is.array(Expansion_cz) || !(all(dim(Expansion_cz) == c(n_c, 2)))) {
      stop("`Expansion_cz` has wrong dimensions")
    }
  }
  if (is.null(Xconfig_zcp)) {
    Xconfig_zcp = array(1, dim = c(2, n_c, n_p))
  } else {
    if (!is.array(Xconfig_zcp) || !(all(dim(Xconfig_zcp) == c(2, n_c, n_p)))) {
      stop("`Xconfig_zcp` has wrong dimensions")
    }
    if (!all(Xconfig_zcp %in% c(0, 1, 2, 3))) {
      stop("`Xconfig_zcp` has some wrong element(s)")
    }
  }
  FieldConfig_input = array(NA, dim = dim(FieldConfig), dimnames = dimnames(FieldConfig))
  g = function(mat) suppressWarnings(array(as.numeric(mat), dim = dim(mat)))
  FieldConfig_input[] = ifelse(FieldConfig == "AR1", 0, FieldConfig_input)
  FieldConfig_input[] = ifelse(FieldConfig == "IID", -2, FieldConfig_input)
  FieldConfig_input[] = ifelse(!is.na(g(FieldConfig)) & g(FieldConfig) > 0 & g(FieldConfig) <= n_c, g(FieldConfig), FieldConfig_input)
  FieldConfig_input[] = ifelse(!is.na(g(FieldConfig)) & g(FieldConfig) == 0, -1, FieldConfig_input)
  if (any(is.na(FieldConfig_input)))
    stop("'FieldConfig' must be: 0 (turn off overdispersion); 'IID' (independent for each factor); 'AR1' (use AR1 structure); or 0<n_f<=n_c (factor structure)")
  message("FieldConfig_input is:")
  print(FieldConfig_input)
  OverdispersionConfig_input = rep(NA, length(OverdispersionConfig))
  names(OverdispersionConfig_input) = names(OverdispersionConfig)
  g = function(vec) suppressWarnings(as.numeric(vec))
  OverdispersionConfig_input[] = ifelse(OverdispersionConfig == "AR1", 0, OverdispersionConfig_input)
  OverdispersionConfig_input[] = ifelse(!is.na(g(OverdispersionConfig)) & g(OverdispersionConfig) > 0 & g(OverdispersionConfig) <= n_c, g(OverdispersionConfig), OverdispersionConfig_input)
  OverdispersionConfig_input[] = ifelse(!is.na(g(OverdispersionConfig)) & g(OverdispersionConfig) == 0, -1, OverdispersionConfig_input)
  if (all(OverdispersionConfig_input < 0)) {
    v_i = rep(0, length(b_i))
    n_v = 1
  }
  if (any(is.na(OverdispersionConfig_input)))
    stop("'OverdispersionConfig' must be: 0 (turn off overdispersion); 'AR1' (use AR1 structure); or 0<n_f<=n_c (factor structure)")
  message("OverdispersionConfig_input is:")
  print(OverdispersionConfig_input)
  if (Options2use["Calculate_Range"] == FALSE) {
    Z_gm = Z_xm = matrix(0, nrow = nrow(a_gl), ncol = ncol(a_gl))
  } else {
    if (is.null(spatial_list)) {
      Z_xm = Z_gm = MeshList$loc_x
    } else {
      Z_xm = spatial_list$loc_x
      if (is.null(Z_gm))
        Z_gm = spatial_list$loc_g
    }
    message("Calculating range shift for stratum #1:", colnames(a_gl[1]))
  }
  if (CheckForErrors == TRUE) {
    if (ncol(ObsModel_ez) != 2 | nrow(ObsModel_ez) != n_e)
      stop("Wrong dimensions for ObsModel_ez")
    if (!is.matrix(a_gl) | !is.matrix(X_xj) | !is.matrix(Q_ik))
      stop("a_gl, X_xj, and Q_ik should be matrices")
    if (!is.array(X_xtp) | !is.array(X_gtp) | !is.array(X_itp))
      stop("X_xtp, X_gtp, and X_itp should be arrays")
    if ((max(s_i) - 1) > n_x | min(s_i) < 0)
      stop("s_i exceeds bounds in MeshList")
    if (any(a_i <= 0))
      stop("a_i must be greater than zero for all observations, and at least one value of a_i is not")
    Prop_nonzero = tapply(b_i, INDEX = list(tprime_iz[, 1], c_iz[, 1]),
                          FUN = function(vec) {
                            mean(vec > 0)})
    if (Options2use[12] == 1) {
      Prop_nonzero = Prop_nonzero[-1, ]
    }
    if (any(ObsModel_ez[, 2] %in% c(0, 1))) {
      if (RhoConfig[1] == 0) {
        if (any(!is.na(Prop_nonzero) & (Prop_nonzero == 1))) {
          print(Prop_nonzero)
          stop("Some years and/or categories have 100% encounters, and this requires either temporal structure of a different link-function")
        }
        if (any(!is.na(Prop_nonzero) & (Prop_nonzero == 0)) & Options2use["treat_nonencounter_as_zero"] == FALSE) {
          print(Prop_nonzero)
          stop("Some years and/or categories have 0% encounters, and this requires either temporal structure of a different link-function")
        }
      }
    }
    if (length(OverdispersionConfig) != 2)
      stop("length(OverdispersionConfig)!=2")
    if (ncol(yearbounds_zz) != 2)
      stop("yearbounds_zz must have two columns")
    if (Options2use["Calculate_Coherence"] == 1 & any(ObsModel_ez[, 2] == 0))
      stop("Calculating coherence only makes sense when 'ObsModel_ez[,2]=1'")
    if (any(yearbounds_zz < 0) | any(yearbounds_zz >= max(n_t)))
      stop("yearbounds_zz exceeds bounds for modeled years")
    if (ncol(t_yz) != ncol(tprime_iz))
      stop("t_yz and tprime_iz must have same number of columns")
    if (n_c != length(unique(na.omit(as.vector(c_iz)))))
      stop("n_c doesn't equal the number of levels in c_i")
    if (any(X_xj != 0))
      stop("X_xj is deprecated, please use X_xtp to specify static or dynamic density covariates (which by default have constant effect among years but differ among categories)")
    if (any(ObsModel_ez[, 1] == 9) & !all(b_i %in% 0:3))
      stop("If using 'ObsModel_ez[e,1]=9', all 'b_i' must be in {0,1,2,3}")
    if (length(unique(ObsModel_ez[, 2])) > 1)
      stop("All `ObsModel_ez[,2]` must have the same value")
    if (any(OverdispersionConfig > 0) & length(unique(v_i)) == 1)
      stop("It doesn't make sense to use use `OverdispersionConfig` when using only one level of `v_i`")
    if (any(ObsModel_ez[, 1] %in% c(12, 13, 14))) {
      if (any(ObsModel_ez[, 2] != 1))
        stop("If using `ObsModel_ez[e,1]` in {12,13,14} then must use `ObsModel_ez[e,2]=1`")
      if (!any(ObsModel_ez[, 1] %in% c(0, 1, 2, 3)))
        stop("Using `ObsModel_ez[e,1]` in {12,13,14} is only intended when combining data with biomass-sampling data")
    }
    if (all(b_i > 0) & all(ObsModel_ez[, 1] == 0) & !all(FieldConfig_input[1:2, 1] == -1))
      stop("All data are positive and using a conventional delta-model, so please turn off `Omega1` and `Epsilon1` terms")
    if (!(all(ObsModel_ez[, 1] %in% c(0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))))
      stop("Please check `ObsModel_ez[,1]` input")
    if (!(all(ObsModel_ez[, 2] %in% c(0, 1, 2, 3, 4))))
      stop("Please check `ObsModel_ez[,2]` input")
    if (!all(RhoConfig[1] %in% c(0, 1, 2, 3, 4)) | !all(RhoConfig[2] %in% c(0, 1, 2, 3, 4, 6)) | !all(RhoConfig[3] %in% c(0, 1, 2, 4, 5)) | !all(RhoConfig[4] %in% c(0, 1, 2, 4, 5, 6)))
      stop("Check `RhoConfig` inputs")
    if (any(is.na(X_xtp)))
      stop("Some `X_xtp` is NA, and this is not allowed")
    if (any(is.na(X_gtp)))
      stop("Some `X_gtp` is NA, and this is not allowed")
    if (any(is.na(X_itp)))
      stop("Some `X_itp` is NA, and this is not allowed")
    if (n_c == 1 && !all(FieldConfig_input %in% c(-2, -1, 1)))
      stop("If using a univariate model, `FieldConfig` must be 0, 1, or `IID` for all entries")
  }
  if (CheckForErrors == TRUE) {
    if (any(c(length(b_i), length(a_i), nrow(c_iz), length(s_i), nrow(tprime_iz), length(v_i), length(PredTF_i)) != n_i))
      stop("b_i, a_i, c_i, s_i, v_i, or tprime_i doesn't have length n_i")
    if (nrow(a_gl) != n_x | ncol(a_gl) != n_l)
      stop("a_xl has wrong dimensions")
    if (nrow(X_xj) != n_x | ncol(X_xj) != n_j)
      stop("X_xj has wrong dimensions")
    if (nrow(Q_ik) != n_i | ncol(Q_ik) != n_k)
      stop("Q_ik has wrong dimensions")
    if (FishStatsUtils::convert_version_name(Version) >=
        FishStatsUtils::convert_version_name("VAST_v8_0_0")) {
      if (dim(X_gtp)[1] != n_g | dim(X_gtp)[2] != n_t |
          dim(X_gtp)[3] != n_p)
        stop("X_gtp has wrong dimensions")
      if (dim(X_itp)[1] != n_i | dim(X_itp)[2] != n_t |
          dim(X_itp)[3] != n_p)
        stop("X_itp has wrong dimensions")
    } else {
      if (dim(X_xtp)[1] != n_x | dim(X_xtp)[2] != n_t |
          dim(X_xtp)[3] != n_p)
        stop("X_xtp has wrong dimensions")
    }
    if (ncol(c_iz) > 1 & any(ObsModel_ez[, 2] != 1))
      stop("Using multiple columnns in `c_iz` only makes sense using a Poisson-link delta model via `ObsModel[2]=1`")
    if (nrow(F_ct) != n_c | ncol(F_ct) != n_t)
      stop("F_ct has wrong dimensions")
  }
  if (FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0")) {
    if (is.null(spatial_list))
      stop("Must provide `spatial_list` for Version >= 8.0.0")
  }
  if (any(ObsModel_ez[, 2] == 1)) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v2_2_0")) {
      stop("Problem with Poisson-link model for VAST versions 1.6.0 through 2.2.0")
    }
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v1_5_0")) {
      stop("Poisson-link model for VAST versions 1.0.0 through 1.5.0")
    }
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v2_8_0")) {
      if (length(unique(ObsModel_ez[, 1])) > 1)
        stop("Can't use multiple error distributions prior to version 3.0.0")
    }
  }
  if (any(FieldConfig_input == -2)) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v2_5_0")) {
      stop("Problem with 'IID' option for factors for VAST versions 1.0.0 through 2.6.0")
    }
  }
  if (ncol(c_iz) > 1) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v3_0_0")) {
      stop("Can't have observations assigned to more than one category prior to version 4.0.0")
    }
  }
  if (any(RhoConfig[1:2] == 3)) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v4_1_0")) {
      stop("There was bug in specifying fixed intercepts among years for versions prior to V4.2.0")
    }
  }
  if (Options2use["SD_observation_density"] == 1) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v4_1_0")) {
      stop("Calculating 'SD_observation_density' is not possible prior to V4.2.0")
    }
  }
  if (any(ObsModel_ez[, 1] == 3)) {
    if (FishStatsUtils::convert_version_name(Version) <=
        FishStatsUtils::convert_version_name("VAST_v8_2_0")) {
      stop("Inverse-gaussian distribution only available for CPP version >= 8_3_0")
    }
  }
  if (ncol(tprime_iz) >= 2 & any(RhoConfig[1:2] != 0)) {
    stop("Temporal structure on intercepts is not implemented for seasonal models")
  }
  if (ncol(tprime_iz) >= 2 & any(VamConfig[1] != 0)) {
    stop("Species interactions are not implemented for seasonal models")
  }
  if ((FieldConfig_input[2, 1] == (-1) & RhoConfig[3] != 0) |
      (FieldConfig_input[2, 2] == (-1) & RhoConfig[4] != 0)) {
    stop("Spatio-temporal variation is turned off for a component with temporal structure, and this combination doesn't make sense")
  }
  if (!is.null(spatial_list$fine_scale) && spatial_list$fine_scale ==
      TRUE) {
    if (Options2use["SD_site_density"] == TRUE | Options2use["SD_site_logdensity"] == TRUE) {
      warning("'SD_site_density' and 'SD_site_logdensity' are very slow when using `fine_scale=TRUE`")
    }
  }
  if (VamConfig[1] != 0) {
    if (any(ObsModel_ez[, 2] != 1)) {
      stop("Must use Poisson-link delta model when estimating interactions")
    }
    if (any(RhoConfig[1:2] != 3)) {
    }
    if (!(FieldConfig_input[2, 1] %in% c(-2, -1, 0, n_c)) & VamConfig[3] == 1) {
      stop("Spatio-temporal variation must either have full rank covariance or be turned off or IID for the 1st linear predictor for when using interactions and when VamConfig[`Timing`]==1")
    }
    if (VamConfig[2] > n_c | VamConfig[2] < 0) {
      stop("Rank for interactions must be an integer between 0 and `n_c`, inclusive")
    }
    if (VamConfig[2] == FieldConfig_input[2, 1] & RhoConfig[3] != 0) {
      stop("Can't simultaneously identify full-rank interactions and temporal correlation on spatio-temporal component for 1st linear predictor")
    }
  }
  if (RhoConfig[4] == 6) {
    if (FieldConfig_input[2, 1] != FieldConfig_input[2, 2]) {
      stop("To fix 'Epsilon_rho2_f` equal to 'Epsilon_rho2_f`, you must specify the same rank using `FieldConfig_input[2,1]` equal to `FieldConfig_input[2,2]`")
    }
  }
  if (Options2use[12] != 0) {
    if (FieldConfig_input[2, 1] != n_c) {
      stop("Must have full-rank spatio-temporal component matrix to estimate B0 using `Options2use[12]=1`")
    }
    if (!(Options2use[12] %in% c(0, 1))) {
      stop("`Options2use[12]` must be either 0 or 1")
    }
    if (any(ObsModel_ez[, 2] != 1)) {
      stop("Must use Poisson-link delta model when estimating interactions")
    }
  }
  if (Options2use[11] != 0) {
    if (FieldConfig_input[2, 2] != 0 & !(RhoConfig[4] %in% c(6))) {
      stop("To estimate Fratio, either Epsilon2 must be turned off (i.e., `FieldConfig_input[2,2]=0`) or B2 must equal B1_cc (i.e., `RhoConfig[4]=6`)")
    }
  }
  if (any(F_ct != 0)) {
    if (any(ObsModel_ez[, 2] != 1)) {
      stop("Must use Poisson-link delta model when estimating the impact of fishing mortality")
    }
    if (any(RhoConfig[1:2] != 3)) {
      stop("Must use constant intercepts when estimating the impact of fishing mortality")
    }
    if (!(F_init %in% c(1, 2))) {
      stop("`F_init` must be either 1 or 2")
    }
    if (F_init == 2) {
      if (FieldConfig_input[2, 2] != 0 & !(RhoConfig[4] %in% c(6))) {
        stop("To estimate stationary distribution for initial F, either Epsilon2 must be turned off (i.e., `FieldConfig_input[2,2]=0`) or B2 must equal B1_cc (i.e., `RhoConfig[4]=6`)")
      }
    }
  }
  if (Method == "Stream_network") {
    if (is.null(Network_sz)) {
      stop("Must specify 'Network_sz' when using Method=='Stream_network'")
    }
    if (ncol(Network_sz) != 3 | !all(c("parent_s", "child_s", "dist_s") %in% colnames(Network_sz))) {
      stop("'Network_sz' must have three columns, 'parent_s', 'child_s', and 'dist_s'")
    }
  } else {
    Network_sz = matrix(c(1, 1, 1), nrow = 1, dimnames = list(NULL, c("parent_s", "child_s", "dist_s")))
  }
  if (any(RhoConfig[1:2] != 0) & any(ObsModel_ez[, 2] == 3)) {
    stop("RhoConfig[1:2] must be 0 when using ObsModel[2]=3:  Other options are not coded to work together")
  }
  if (any(RhoConfig[1:2] != 0) & any(ObsModel_ez[, 2] == 4)) {
    stop("RhoConfig[1:2] must be 0 when using ObsModel[2]=4: Other options are not coded to work together")
  }
  if (any(FieldConfig_input[3, 1:2] != -2) & any(ObsModel_ez[, 2] %in% c(3, 4))) {
    stop("Factor model for intercepts is incompatible  with ObsModel_ez[,2] being 3 or 4")
  }
  if (RhoConfig[1] == 0) {
    if (FieldConfig_input[3, 1] != -2)
      stop("Using a factor model doesn't make sense using fixed-effect intercepts.  If you want to use a factor model without temporal structure, please change `RhoConfig[1]=1` for covariance that is independent in each year, or use some other temporal structure on intercepts.")
  }
  if (RhoConfig[2] == 0) {
    if (FieldConfig_input[3, 2] != -2)
      stop("Using a factor model doesn't make sense using fixed-effect intercepts.  If you want to use a factor model without temporal structure, please change `RhoConfig[2]=1` for covariance that is independent in each year, or use some other temporal structure on intercepts.")
  }
  if (Method == "Grid") {
    Aniso = 0
    message("Using isotropic 2D AR1 hyperdistribution, so switching to Aniso=0")
  }
  if (Method == "Spherical_mesh") {
    Aniso = 0
    message("Using spherical projection for SPDE approximation, so switching to Aniso=0")
  }
  if (Method == "Stream_network") {
    Aniso = 0
    message("Using stream network correlations, so switching to Aniso=0")
  }
  if (VamConfig[2] == 0 & VamConfig[1] != 0) {
    VamConfig[1] = 0
    message("Using interactions with zero rank (`VamConfig[2]==0`), so turning off interactions (`VamConfig[1]=0`)")
  }
  if (n_c > 1 & any(FieldConfig_input == 1)) {
    warning("Using 1 factor for more than one category:  Please note that this is non-standard, and it is more common to use multiple factors (often as many as the number of categories)")
  }
  SD_p = apply(X_xtp, MARGIN = 3, FUN = sd)
  if (any(SD_p > 3)) {
    warning("I highly recommend that you standardize each density covariate `X_xtp` to have a low standard deviation, to avoid numerical under/over-flow")
  }
  Options_vec = c(Aniso = Aniso, R2_interpretation = 0, Rho_beta1_TF = ifelse(RhoConfig[1] %in% c(1, 2, 4), 1, 0), Rho_beta2_TF = ifelse(RhoConfig[2] %in% c(1, 2, 4), 1, 0), AreaAbundanceCurveTF = 0, CMP_xmax = 200, CMP_breakpoint = 1, Method = switch(Method, Mesh = 0, Grid = 1, Spherical_mesh = 0, Stream_network = 2), Include_F = ifelse(all(F_ct == 0), 0, F_init))
  Return = NULL
  if (Version %in% c("VAST_v1_1_0", "VAST_v1_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2,]), ObsModel = ObsModel_ez[1, ], Options = Options2use, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list())
  }
  if (Version %in% c("VAST_v1_4_0", "VAST_v1_3_0", "VAST_v1_2_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_f_input = OverdispersionConfig_input[1], n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), ObsModel = ObsModel_ez[1, ], Options = Options2use, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], v_i = match(v_i, sort(unique(v_i))) - 1, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list())
  }
  if (Version %in% c("VAST_v1_6_0", "VAST_v1_5_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_f_input = OverdispersionConfig_input[1], n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), ObsModel = ObsModel_ez[1, ], Options = Options2use, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list())
  }
  if (Version %in% c("VAST_v1_7_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel = ObsModel_ez[1, ], Options = Options2use, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list())
  }
  if (Version %in% c("VAST_v1_8_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel = ObsModel_ez[1, ], Options = Options2use, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v1_9_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel = ObsModel_ez[1, ], Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_i = tprime_iz[, 1], v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v2_8_0", "VAST_v2_7_0", "VAST_v2_6_0",
                     "VAST_v2_5_0", "VAST_v2_4_0", "VAST_v2_3_0", "VAST_v2_2_0",
                     "VAST_v2_1_0", "VAST_v2_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel = ObsModel_ez[1, ], Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v3_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_i = c_iz[, 1], e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v4_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_j = n_j, n_p = n_p, n_k = n_k, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v4_4_0", "VAST_v4_3_0", "VAST_v4_2_0",
                     "VAST_v4_1_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, include_data = TRUE, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v5_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, include_data = TRUE, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v5_1_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, include_data = TRUE, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v5_2_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, include_data = TRUE, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v5_4_0", "VAST_v5_3_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_vec = Options_vec, FieldConfig = as.vector(FieldConfig_input[1:2, ]), RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, include_data = TRUE, Options = Options2use, yearbounds_zz = yearbounds_zz, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v5_5_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x,n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_list = list(Options_vec = Options_vec, Options = Options2use, yearbounds_zz = yearbounds_zz, Expansion_cz = Expansion_cz), FieldConfig = as.vector(FieldConfig_input[1:2,]), RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, include_data = TRUE, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xj = X_xj, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v6_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_list = list(Options_vec = Options_vec, Options = Options2use, yearbounds_zz = yearbounds_zz, Expansion_cz = Expansion_cz), FieldConfig = as.vector(FieldConfig_input[1:2, ]), RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, Xconfig_zcp = Xconfig_zcp, include_data = TRUE, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v7_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_x = n_x, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_xm), Options_list = list(Options_vec = Options_vec, Options = Options2use, yearbounds_zz = yearbounds_zz, Expansion_cz = Expansion_cz), FieldConfig = FieldConfig_input, RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, Xconfig_zcp = Xconfig_zcp, include_data = TRUE, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, s_i = s_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_xl = a_gl, X_xtp = X_xtp, Q_ik = Q_ik, t_yz = t_yz, Z_xm = Z_xm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2)
  }
  if (Version %in% c("VAST_v8_5_0", "VAST_v8_3_0", "VAST_v8_2_0", "VAST_v8_1_0",
                     "VAST_v8_0_0")) {
    Return = list(n_i = n_i, n_s = c(MeshList$anisotropic_spde$n.spde, n_x, n_x)[Options_vec["Method"] + 1], n_g = n_g, n_t = n_t, n_c = n_c, n_e = n_e, n_p = n_p, n_v = n_v, n_l = n_l, n_m = ncol(Z_gm), Options_list = list(Options_vec = Options_vec, Options = Options2use, yearbounds_zz = yearbounds_zz, Expansion_cz = Expansion_cz), FieldConfig = FieldConfig_input, RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig_input, ObsModel_ez = ObsModel_ez, VamConfig = VamConfig, Xconfig_zcp = Xconfig_zcp, include_data = TRUE, b_i = b_i, a_i = a_i, c_iz = c_iz, e_i = e_i, t_iz = tprime_iz, v_i = match(v_i, sort(unique(v_i))) - 1, PredTF_i = PredTF_i, a_gl = a_gl, X_itp = X_itp, X_gtp = X_gtp, Q_ik = Q_ik, t_yz = t_yz, Z_gm = Z_gm, F_ct = F_ct, parent_s = Network_sz[, "parent_s"] - 1, child_s = Network_sz[, "child_s"] - 1, dist_s = Network_sz[, "dist_s"], spde = list(), spde_aniso = list(), M0 = GridList$M0, M1 = GridList$M1, M2 = GridList$M2, Ais_ij = cbind(spatial_list$A_is@i, spatial_list$A_is@j), Ais_x = spatial_list$A_is@x, Ags_ij = cbind(spatial_list$A_gs@i, spatial_list$A_gs@j), Ags_x = spatial_list$A_gs@x)
  }
  if (is.null(Return))
    stop("`Version` provided does not match the list of possible values")
  if ("spde" %in% names(Return))
    Return[["spde"]] = MeshList$isotropic_spde$param.inla[c("M0", "M1", "M2")]
  if ("spde_aniso" %in% names(Return))
    Return[["spde_aniso"]] = list(n_s = MeshList$anisotropic_spde$n.spde, n_tri = nrow(MeshList$anisotropic_mesh$graph$tv), Tri_Area = MeshList$Tri_Area, E0 = MeshList$E0, E1 = MeshList$E1, E2 = MeshList$E2, TV = MeshList$TV - 1, G0 = MeshList$anisotropic_spde$param.inla$M0, G0_inv = INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)))
  if (CheckForErrors == TRUE) {
    NoNAs = setdiff(names(Return), c("t_iz", "t_yz", "c_iz", "Network_sz", "b_i"))
    if (any(sapply(Return[NoNAs], FUN = function(Array) {
      any(is.na(Array))
    }) == TRUE))
      stop("Please find and eliminate the NA from your inputs")
  }
  class(Return) = "make_data"
  return(Return)
}

make_map_info_zz<- function (Region, Extrapolation_List, spatial_list = NULL, NN_Extrap = spatial_list$PolygonList$NN_Extrap, fine_scale = spatial_list$fine_scale, Include)
{
  if (is.null(fine_scale))
    fine_scale = FALSE
  if (is.null(spatial_list)) {
    if (fine_scale == FALSE) {
      warning("Consider updating inputs to `make_map_info` to enable future use of feature `fine_scale=TRUE`")
    } else {
      stop("Must update inputs to `make_map_info` to enable feature `fine_scale=TRUE`")
    }
  }
  if (missing(Include)) {
    Include = Extrapolation_List[["Area_km2_x"]] > 0
  }
  PlotDF = NULL
  if (tolower(Region) == "california_current") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("state", c("alabama", "arizona", "arkansas", "california", "colorado", "connecticut", "delaware", "district of columbia", "florida", "georgia", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts:martha's vineyard", "massachusetts:main", "massachusetts:nantucket", "michigan:north", "michigan:south", "minnesota", "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", "new jersey", "new mexico", "new york:manhattan", "new york:main", "new york:statenisland", "new york:longisland", "north carolina:knotts", "north carolina:main", "north carolina:spit", "north dakota", "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", "south dakota", "tennessee", "texas", "utah", "vermont", "virginia:chesapeake", "virginia:chincoteague", "virginia:main", "washington:san juan island", "washington:lopez island", "washington:orcas island", "washington:whidbey island", "washington:main", "west virginia", "wisconsin", "wyoming"))
    Xlim = c(-126, -117)
    Ylim = c(32, 49)
    Rotate = 20
    Cex = 0.01
    Legend = list(use = TRUE, x = c(65, 75), y = c(35, 65))
  }
  if (tolower(Region) == "british_columbia") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = c(-133, -126)
    Ylim = c(50, 55)
    Rotate = 0
    Cex = 0.1
    Legend = list(use = FALSE, x = c(5, 10), y = c(5, 45))
  }
  if (tolower(Region) == "eastern_bering_sea") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = c(-180, -158)
    Ylim = c(54, 63)
    Rotate = 0
    Cex = 0.01
    Legend = list(use = TRUE, x = c(76, 86), y = c(48, 83))
  }
  if (tolower(Region) == "aleutian_islands") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    PlotDF[, "Lon"] = PlotDF[, "Lon"]%%360
    MappingDetails = list("world2Hires", NULL)
    Xlim = c(170, 195)
    Ylim = c(51, 55)
    Rotate = 0
    Cex = 0.2
    Legend = list(use = FALSE, x = c(5, 10), y = c(5, 45))
  }
  if (tolower(Region) == "gulf_of_alaska") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("world", NULL)
    Xlim = c(-171, -132)
    Ylim = c(52, 61)
    Rotate = 0
    Cex = 0.01
    Legend = list(use = TRUE, x = c(5, 10), y = c(30, 65))
  }
  if (tolower(Region) == "northwest_atlantic") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("state", c("alabama", "arizona", "arkansas", "california", "colorado", "connecticut", "delaware", "district of columbia", "florida", "georgia", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts:martha's vineyard", "massachusetts:main", "massachusetts:nantucket", "michigan:north", "michigan:south", "minnesota", "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", "new jersey", "new mexico", "new york:manhattan", "new york:main", "new york:statenisland", "new york:longisland", "north carolina:knotts", "north carolina:main", "north carolina:spit", "north dakota", "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", "south dakota", "tennessee", "texas", "utah", "vermont", "virginia:chesapeake", "virginia:chincoteague", "virginia:main", "washington:san juan island", "washington:lopez island", "washington:orcas island", "washington:whidbey island", "washington:main", "west virginia", "wisconsin", "wyoming"))
    Xlim = c(-80, -65)
    Ylim = c(32, 45)
    Rotate = 0
    Cex = 0.01
    Legend = list(use = TRUE, x = c(70, 80), y = c(5, 35))
  }
  if (tolower(Region) == "south_africa") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = c(14, 26)
    Ylim = c(-37, -28)
    Rotate = 0
    Cex = 0.1
    Legend = list(use = FALSE, x = c(5, 10), y = c(4, 45))
  }
  if (tolower(Region) == "gulf_of_st_lawrence") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", "Canada")
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lon"])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lat"])
    Rotate = 0
    Cex = 1
    Legend = list(use = FALSE, x = c(5, 10), y = c(4, 45))
  }
  if (tolower(Region) == "new_zealand") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = c(172, 187)
    Ylim = c(-46, -41)
    Rotate = 0
    Cex = 0.01
    Legend = list(use = FALSE, x = c(5, 10), y = c(5, 45))
  }
  if (tolower(Region) == "habcam") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lon"])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lat"])
    Rotate = 20
    Cex = 0.01
    Legend = list(use = TRUE, x = c(70, 90), y = c(5, 35))
  }
  if (tolower(Region) == "gulf_of_mexico") {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lon"])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lat"])
    Rotate = 0
    Cex = 0.6
    Legend = list(use = TRUE, x = c(76, 86), y = c(48, 83))
  }
  if (is.null(PlotDF)) {
    PlotDF = cbind(Extrapolation_List[["Data_Extrap"]][, c("Lat", "Lon")], x2i = NA, Include = Include)
    MappingDetails = list("worldHires", NULL)
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lon"])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(Extrapolation_List[["Area_km2_x"]] > 0), "Lat"])
    Rotate = 0
    Cex = 0.1
    Legend = list(use = FALSE, x = c(5, 10), y = c(5, 45))
  }
  if (fine_scale == TRUE | spatial_list$Method == "Stream_network") {
    PlotDF[, "x2i"] = NA
    PlotDF[which(Extrapolation_List[["Area_km2_x"]] > 0), "x2i"] = 1:length(which(Extrapolation_List[["Area_km2_x"]] > 0))
  } else {
    PlotDF[, "x2i"] = NN_Extrap$nn.idx[, 1]
  }
  if (tolower(Region) == "eastern_bering_sea") {
    PlotDF = PlotDF[which(PlotDF[, "Lon"] < 0), ]
  }
  if (is.numeric(Extrapolation_List$zone)) {
    Zone = Extrapolation_List[["zone"]] - ifelse(Extrapolation_List$flip_around_dateline == TRUE, 30, 0)
  } else {
    Zone = Extrapolation_List$zone
  }
  if (all(!is.na(Extrapolation_List$Data_Extrap[, c("N_km", "E_km")]))) {
    MapSizeRatio = c(`Height(in)` = diff(range(Extrapolation_List$Data_Extrap[, "N_km"])), `Width(in)` = diff(range(Extrapolation_List$Data_Extrap[, "E_km"])))
  } else {
    MapSizeRatio = c(`Height(in)` = diff(range(Extrapolation_List$Data_Extrap[, "Lat"])), `Width(in)` = diff(range(Extrapolation_List$Data_Extrap[, "Lon"])))
  }
  MapSizeRatio = MapSizeRatio/sqrt(prod(MapSizeRatio)) * 4
  mapdetails_list = list(PlotDF = PlotDF, MappingDetails = MappingDetails, Xlim = Xlim, Ylim = Ylim, MapSizeRatio = MapSizeRatio, Rotate = Rotate, Cex = Cex, Legend = Legend, Zone = Zone, fine_scale = fine_scale)
  return(mapdetails_list)
}

# “Plot” function edits ---------------------------------------------------
plot_results_zz<- function (fit, settings = fit$settings, plot_set = 3, working_dir = paste0(OutFile, "/"), year_labels = fit$year_labels, years_to_plot = fit$years_to_plot, use_biascorr = TRUE, map_list, category_names, check_residuals = TRUE, projargs = "+proj=longlat", projargs2 = "+init=epsg:32619", zrange, n_samples = 100, ...) {

  if(FALSE){
    fit<- fit.seasind.spattemp.aja
    settings = fit$settings
    plot_set = c(3, 6, 7, 13, 14, 15, 16)
    working_dir = paste(OutFile.aja, "/", sep = "")
    year_labels = fit$year_labels
    years_to_plot = fit$years_to_plot
    use_biascorr = TRUE
    check_residuals = FALSE
    projargs = "+proj=longlat"
    projargs2 = "+init=epsg:32619"
    n_samples = 100

    ... = NULL
  }

  if (is.null(fit$Report))
    stop("`fit$Report` is missing, please check inputs")
  if (missing(category_names))
    category_names = 1:fit$data_list$n_c
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  message("\n### Making plots of data availability and knots")
  plot_data_zz(Extrapolation_List = fit$extrapolation_list, Spatial_List = fit$spatial_list, Lat_i = fit$data_frame[, "Lat_i"], Lon_i = fit$data_frame[, "Lon_i"], Year_i = fit$data_frame[, "t_i.Date"], projargs = "+proj=longlat",  projargs2 = "+init=epsg:32619", PlotDir = working_dir, Year_Set = year_labels)
  if (missing(map_list)) {
    message("\n### Obtaining default settings for plotting maps")
    map_list = make_map_info_zz(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
  }
  message("\n### Making plot of encounter probability")
  Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = cbind(Catch_KG = fit$data_frame[, "b_i"]), DirName = working_dir)
  message("\n### Making plot of anisotropy")
  plot_anisotropy(FileName = paste0(working_dir, "Aniso.png"), Report = fit$Report, TmbData = fit$data_list)
  if (!is.null(fit$parameter_estimates$SD)) {
    message("\n### Making plot of abundance index")
    Index = plot_biomass_index(DirName = working_dir, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Year_Set = year_labels, Years2Include = years_to_plot, use_biascorr = use_biascorr, category_names = category_names)
  } else {
    Index = "Not run"
    message("\n### Skipping plot of abundance index; must re-run with standard errors to plot")
  }
  if (!is.null(fit$parameter_estimates$SD)) {
    message("\n### Making plot of spatial indices")
    Range = plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = working_dir, Year_Set = year_labels, Years2Include = years_to_plot, use_biascorr = use_biascorr, category_names = category_names)
  } else {
    Range = "Not run"
    message("\n### Skipping plot of spatial indices; must re-run with standard errors to plot")
  }
  if ("jointPrecision" %in% names(fit$parameter_estimates$SD) & n_samples > 0) {
    message("\n### Making plot of spatial indices")
    Edge = plot_range_edge(Obj = fit$tmb_list$Obj, Sdreport = fit$parameter_estimates$SD, working_dir = working_dir, Year_Set = year_labels, Years2Include = years_to_plot, category_names = category_names, n_samples = n_samples, quantiles = c(0.05, 0.5, 0.95))
  } else {
    Edge = "Not run"
    message("\n### Skipping plot of range edge; only possible if `getJointPrecision=TRUE` and `n_samples`>0")
  }
  message("\n### Making plots of spatial predictions")
  plot_maps_args = list(...)
  plot_maps_args = combine_lists(input = plot_maps_args, default = list(plot_set = plot_set, category_names = category_names, TmbData = fit$data_list, Report = fit$Report, Sdreport = fit$parameter_estimates$SD, PlotDF = map_list[["PlotDF"]], MapSizeRatio = map_list[["MapSizeRatio"]], working_dir = working_dir, Year_Set = year_labels, Years2Include = years_to_plot, legend_x = map_list[["Legend"]]$x/100, legend_y = map_list[["Legend"]]$y/100))
  Dens_xt = do.call(what = plot_maps_zz, args = plot_maps_args)
  if (check_residuals == TRUE) {
    message("\n### Making Q-Q plot")
    Q = plot_quantile_diagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive", FileName_Phist = "Posterior_Predictive-Histogram", FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist", save_dir = working_dir)
    message("\n### Making plot of Pearson residuals")
    plot_residuals_zz(Lat_i = fit$data_frame[, "Lat_i"], Lon_i = fit$data_frame[, "Lon_i"], TmbData = fit$data_list, Report = fit$Report, Q = Q, working_dir = working_dir, spatial_list = fit$spatial_list, extrapolation_list = fit$extrapolation_list, Year_Set = year_labels, Years2Include = years_to_plot, zrange = zrange, legend_x = map_list[["Legend"]]$x/100, legend_y = map_list[["Legend"]]$y/100)
  } else {
    Q = "Not run"
    message("\n### Skipping Q-Q plot")
    message("\n### Skipping plot of Pearson residuals")
  }
  Return = list(Q = Q, Index = Index, Range = Range, Dens_xt = Dens_xt, Edge = Edge, map_list = map_list, plot_maps_args = plot_maps_args)
  return(invisible(Return))
}

plot_data_zz<- function (Extrapolation_List, Spatial_List, Data_Geostat = NULL, Lat_i = Data_Geostat[, "Lat"], Lon_i = Data_Geostat[, "Lon"], Year_i = Data_Geostat[, "Year"], PlotDir = working_dir, Plot1_name = "Data_and_knots.png", Plot2_name = "Data_by_year.png", Year_Set, projargs = "+proj=longlat", projargs2 = "+init=epsg:32619", map_resolution = "medium", land_color = "grey", ...) {
  if(FALSE){
    Extrapolation_List = fit$extrapolation_list
    Spatial_List = fit$spatial_list
    Data_Geostat = NULL
    Lat_i = fit$data_frame[, "Lat_i"]
    Lon_i = fit$data_frame[, "Lon_i"]
    Year_i = fit$data_frame[, "t_i.Date"]
    PlotDir = working_dir
    Year_Set = year_labels
    Plot1_name = "Data_and_knots.png"
    Plot2_name = "Data_by_year.png"
    Year_Set = NULL
    projargs = "+proj=longlat"
    projargs2 = "+init=epsg:26919"
    map_resolution = "medium"
    land_color = "#d9d9d9"
  }

  if (is.null(Lat_i) | is.null(Lon_i) | is.null(Year_i)) {
    stop("Problem with inputs")
  }

  CRS_proj = sp::CRS(projargs)
  CRS_proj2 = st_crs(projargs2)
  map_data = ne_countries(scale='medium', returnclass = 'sf') %>%
    filter(., admin %in% c("United States of America", "Canada"))
  map_data.proj = st_transform(map_data, crs = projargs2)

  if(all(Extrapolation_List[["a_el"]] == 0)){
    which_rows = seq(from = 1, to = nrow(Extrapolation_List[["a_el"]]))
  } else {
    which_rows = which(Extrapolation_List[["Area_km2_x"]] > 0 & rowSums(Extrapolation_List[["a_el"]]) > 0)
  }
  extrap.lonlat<- ggplot() +
    geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
    geom_point(data = Extrapolation_List$Data_Extrap, aes(x = Lon, y = Lat), pch = 21, fill = NA, color = "black", size = 0.15, alpha = 0.5) +
    coord_sf(xlim = range(Extrapolation_List$Data_Extrap[, c("Lon")]), ylim = range(Extrapolation_List$Data_Extrap[, c("Lat")]), expand = FALSE) +
    ggtitle("Extrapolation (Lat-Lon)") +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA))

  extrap.sf<- st_as_sf(Extrapolation_List$Data_Extrap, coords = c("Lon", "Lat"), crs = projargs) %>%
    st_transform(., crs = projargs2)
  extrap.sf.coords<- data.frame(st_coordinates(extrap.sf))

  knots.sf<- st_as_sf(data.frame(Spatial_List$latlon_x), coords = c("Lon", "Lat"), crs = projargs) %>%
    st_transform(., crs = projargs2)

  proj.base<- ggplot() +
    geom_sf(data = map_data.proj, fill = land_color, lwd = 0.2) +
    coord_sf(xlim = range(extrap.sf.coords$X), ylim = range(extrap.sf.coords$Y), expand = FALSE, datum = sf::st_crs(projargs2)) +
    ggtitle("Extrapolation (North-East)") +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA))

  extrap.proj<- proj.base +
    geom_sf(data = extrap.sf, pch = 21, fill = NA, color = "black", size = 0.15, alpha = 0.5) +
    coord_sf(xlim = range(extrap.sf.coords$X), ylim = range(extrap.sf.coords$Y), expand = FALSE, datum = sf::st_crs(projargs2)) +
    ggtitle("Extrapolation (North-East)") +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA))

  knots.proj<- extrap.proj +
    geom_sf(data = knots.sf, pch = 21, fill = "#238b45", size = 1) +
    coord_sf(xlim = range(extrap.sf.coords$X), ylim = range(extrap.sf.coords$Y), expand = FALSE, datum = sf::st_crs(projargs2)) +
    ggtitle("Extrapolation (North-East)") +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA))

  extrap.plot.out<- extrap.lonlat + extrap.proj + knots.proj
  ggsave(filename = paste(working_dir, Plot1_name, sep = ""), extrap.plot.out, height = 8, width = 11, units = "in")

  samp.df<- data.frame("Lon" = Lon_i, "Lat" = Lat_i, "Date" = Year_i)
  samp.plot<- ggplot() +
    geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
    geom_point(data = samp.df, aes(x = Lon, y = Lat), pch = 21, fill = "#2171b5", alpha = 0.5, size = 0.75) +
    coord_sf(xlim = range(Extrapolation_List$Data_Extrap[, c("Lon")]), ylim = range(Extrapolation_List$Data_Extrap[, c("Lat")]), expand = FALSE) +
    scale_x_continuous(breaks = c(-74, -68)) +
    scale_y_continuous(breaks = c(36, 40, 44)) +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA)) +
    facet_wrap(~Date, nrow = 4, ncol = 10)
  ggsave(filename = paste(working_dir, Plot2_name, sep = ""), samp.plot, height = 8, width = 11, units = "in")
}

plot_maps_zz<- function (plot_set = 3, Report, PlotDF, Sdreport = NULL, TmbData = NULL, projargs = "+proj=longlat", projargs2 = "+init=epsg:32619", Panel = "Category", Year_Set = NULL, Years2Include = NULL, category_names = NULL, quiet = FALSE, working_dir = paste0(OutFile, "/"), ...) {
  if(FALSE){
    plot_set = plot_maps_args$plot_set
    category_names = plot_maps_args$category_names
    TmbData = plot_maps_args$data_list
    Report = plot_maps_args$Report
    Sdreport = plot_maps_args$Sdreport
    PlotDF = plot_maps_args$PlotDF
    MapSizeRatio = plot_maps_args$MapSizeRatio
    working_dir = plot_maps_args$working_dir
    Year_Set = plot_maps_args$Year_Set
    Years2Include = plot_maps_args$years_to_plot
    quiet = FALSE
    Panel = "Category"
  }
  if ("D_xt" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:ncol(Report$D_xt)
    if (is.null(Years2Include))
      Years2Include = 1:ncol(Report$D_xt)
    category_names = "singlespecies"
    Ncategories = length(category_names)
    Nyears = dim(Report$D_xt)[2]
  }
  if ("D_xct" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$D_xct)[3]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$D_xct)[3]
    if (is.null(category_names))
      category_names = 1:dim(Report$D_xct)[2]
    Ncategories = dim(Report$D_xct)[2]
    Nyears = dim(Report$D_xct)[3]
  }
  if ("D_xcy" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$D_xcy)[3]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$D_xcy)[3]
    if (is.null(category_names))
      category_names = 1:dim(Report$D_xcy)[2]
    Ncategories = dim(Report$D_xcy)[2]
    Nyears = dim(Report$D_xcy)[3]
  }
  if ("D_gcy" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$D_gcy)[3]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$D_gcy)[3]
    if (is.null(category_names))
      category_names = 1:dim(Report$D_gcy)[2]
    Ncategories = dim(Report$D_gcy)[2]
    Nyears = dim(Report$D_gcy)[3]
  }
  if ("dhat_ktp" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$dhat_ktp)[2]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$dhat_ktp)[2]
    if (is.null(category_names))
      category_names = 1:dim(Report$dhat_ktp)[3]
    Ncategories = dim(Report$dhat_ktp)[3]
    Nyears = dim(Report$dhat_ktp)[2]
  }
  if ("dpred_ktp" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$dpred_ktp)[2]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$dpred_ktp)[2]
    if (is.null(category_names))
      category_names = 1:dim(Report$dpred_ktp)[3]
    Ncategories = dim(Report$dpred_ktp)[3]
    Nyears = dim(Report$dpred_ktp)[2]
  }
  if (Nyears != length(Year_Set)) {
    stop("Problem with `Year_Set`")
  }
  if (Ncategories != length(category_names)) {
    stop("Problem with `category_names`")
  }
  Return = NULL
  for (plot_num in plot_set) {
    Array_xct = NULL
    plot_code <- c("encounter_prob", "pos_catch", "density", "", "", "epsilon_1", "epsilon_2", "linear_predictor_1", "linear_predictor_2", "density_CV", "covariates", "total_density", "covariate_effects_1", "covariate_effects_2", "omega_1", "omega_2")[plot_num]
    if (plot_num == 1) {
      if (quiet == FALSE)
        message(" # Plotting presence/absense maps")
      if ("D_xt" %in% names(Report))
        Array_xct = Report$R1_xt
      if ("D_xct" %in% names(Report))
        Array_xct = Report$R1_xct
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$R1_xcy
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$R1_gcy
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
      message("`plot_num=1` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells")
    }
    if (plot_num == 2) {
      if (quiet == FALSE)
        message(" # Plotting positive catch rate maps")
      if ("D_xt" %in% names(Report))
        Array_xct = log(Report$R2_xt)
      if ("D_xct" %in% names(Report))
        Array_xct = log(Report$R2_xct)
      if ("D_xcy" %in% names(Report))
        Array_xct = log(Report$R2_xcy)
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$R2_gcy
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
      message("`plot_num=2` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells")
    }
    if (plot_num == 3) {
      if (quiet == FALSE)
        message(" # Plotting density maps")
      if ("D_xt" %in% names(Report))
        Array_xct = log(Report$D_xt)
      if ("D_xct" %in% names(Report))
        Array_xct = log(Report$D_xct)
      if ("D_xcy" %in% names(Report))
        Array_xct = log(Report$D_xcy)
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$D_gcy
      if ("dhat_ktp" %in% names(Report))
        Array_xct = aperm(Report$dhat_ktp[, , cI], c(1, 3, 2))
      if ("dpred_ktp" %in% names(Report))
        Array_xct = aperm(Report$dpred_ktp[, , cI], c(1, 3, 2))
    }
    if (plot_num == 6) {
      if (quiet == FALSE)
        message(" # Plotting spatio-temporal effects (Epsilon) in 1st linear predictor")
      if ("D_xt" %in% names(Report))
        Array_xct = Report$Epsilon1_st
      if ("D_xct" %in% names(Report))
        Array_xct = Report$Epsilon1_sct
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$Epsilon1_sct
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$Epsilon1_gct
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
    }
    if (plot_num == 7) {
      if (quiet == FALSE)
        message(" # Plotting spatio-temporal effects (Epsilon) in 2nd linear predictor")
      if ("D_xt" %in% names(Report))
        Array_xct = Report$Epsilon2_st
      if ("D_xct" %in% names(Report))
        Array_xct = Report$Epsilon2_sct
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$Epsilon2_sct
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$Epsilon2_gct
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
    }
    if (plot_num == 8) {
      if (quiet == FALSE)
        message(" # Plotting 1st predictor after action of link function")
      if ("D_xt" %in% names(Report))
        Array_xct = Report$P1_xt
      if ("D_xct" %in% names(Report))
        Array_xct = Report$P1_xct
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$P1_xcy
      if ("D_gcy" %in% names(Report))
        stop("`plot_maps` not implemented for requested plot_num")
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
    }
    if (plot_num == 9) {
      if (quiet == FALSE)
        message(" # Plotting 2nd predictor after action of link function")
      if ("D_xt" %in% names(Report))
        Array_xct = Report$P2_xt
      if ("D_xct" %in% names(Report))
        Array_xct = Report$P2_xct
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$P2_xcy
      if ("D_gcy" %in% names(Report))
        stop("`plot_maps` not implemented for requested plot_num")
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("Not implemented for SpatialVAM")
    }
    if (plot_num == 10) {
      if (quiet == FALSE)
        message(" # Plotting density maps")
      if (is.null(Sdreport))
        stop("Must supply 'Sdreport' if 'plot_num=10'")
      if ("D_xt" %in% names(Report)) {
        if (!("log(Index_xtl)" %in% rownames(TMB::summary.sdreport(Sdreport))))
          stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'SpatialDeltaGLMM'")
        Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == "log(Index_xtl)"), ], dim = c(dim(Report$D_xt), ncol(Report$Index_tl), 2), dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error")))[, , 1, "Std. Error"]
      }
      if ("D_xct" %in% names(Report)) {
        if (!("log(Index_xctl)" %in% rownames(TMB::summary.sdreport(Sdreport))))
          stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == "log(Index_xctl)"), ], dim = c(dim(Report$D_xct), dim(Report$Index_ctl)[3], 2), dimnames = list(NULL, NULL, NULL, NULL, c("Estimate", "Std. Error")))[, , , 1, "Std. Error"]
      }
      if ("D_xcy" %in% names(Report)) {
        if (!("log(Index_xcyl)" %in% rownames(TMB::summary.sdreport(Sdreport))))
          stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == "log(Index_xcyl)"), ], dim = c(dim(Report$D_xcy), dim(Report$Index_cyl)[3], 2), dimnames = list(NULL, NULL, NULL, NULL, c("Estimate", "Std. Error")))[, , , 1, "Std. Error"]
      }
      if (any(c("dhat_ktp", "dpred_ktp") %in% names(Report)))
        stop("'plot_num=10' not implemented for 'SpatialVAM'")
      Array_xct = sqrt(exp(Array_xct^2) - 1)
      if ("D_gcy" %in% names(Report))
        stop("`plot_maps` not implemented for requested plot_num")
    }
    if (plot_num == 11) {
      if (quiet == FALSE)
        message(" # Plotting covariates")
      if (is.null(TmbData))
        stop("Must provide `TmbData` to plot covariates")
      if ("X_xtp" %in% names(TmbData))
        Array_xct = aperm(TmbData$X_xtp, perm = c(1, 3, 2))
      if ("X_gtp" %in% names(TmbData))
        Array_xct = aperm(TmbData$X_gtp, perm = c(1, 3, 2))
      category_names = 1:dim(Array_xct)[2]
    }
    if (plot_num == 12) {
      if (quiet == FALSE)
        message(" # Plotting total density")
      if ("D_xt" %in% names(Report))
        Array_xct = log(Report$D_xt)
      if ("D_xct" %in% names(Report))
        Array_xct = log(apply(Report$D_xct, FUN = sum, MARGIN = c(1, 3)))
      if ("D_xcy" %in% names(Report))
        Array_xct = log(apply(Report$D_xcy, FUN = sum, MARGIN = c(1, 3)))
      if ("D_gcy" %in% names(Report))
        Array_xct = log(apply(Report$D_gcy, FUN = sum, MARGIN = c(1, 3)))
      logsum = function(vec) {
        max(vec) + log(sum(exp(vec - max(vec))))
      }
      if ("dhat_ktp" %in% names(Report))
        Array_xct = apply(aperm(Report$dhat_ktp, c(1, 3, 2)), FUN = logsum, MARGIN = c(1, 3))
      if ("dpred_ktp" %in% names(Report))
        Array_xct = apply(aperm(Report$dpred_ktp, c(1, 3, 2)), FUN = logsum, MARGIN = c(1, 3))
    }
    if (plot_num == 13) {
      if (quiet == FALSE)
        message(" # Plotting covariate effects for 1st linear predictor")
      if ("D_xt" %in% names(Report))
        stop()
      if ("D_xct" %in% names(Report))
        stop()
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$eta1_xct
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$eta1_gct
      if ("dhat_ktp" %in% names(Report))
        stop()
      if ("dpred_ktp" %in% names(Report))
        stop()
    }
    if (plot_num == 14) {
      if (quiet == FALSE)
        message(" # Plotting covariate effects for 2nd linear predictor")
      if ("D_xt" %in% names(Report))
        stop()
      if ("D_xct" %in% names(Report))
        stop()
      if ("D_xcy" %in% names(Report))
        Array_xct = Report$eta2_xct
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$eta2_gct
      if ("dhat_ktp" %in% names(Report))
        stop()
      if ("dpred_ktp" %in% names(Report))
        stop()
    }
    if (plot_num == 15) {
      if (quiet == FALSE)
        message(" # Plotting spatial effects (Omega) in 1st linear predictor")
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$Omega1_gc
    }
    if (plot_num == 16) {
      if (quiet == FALSE)
        message(" # Plotting spatial effects (Omega) in 2nd linear predictor")
      if ("D_gcy" %in% names(Report))
        Array_xct = Report$Omega2_gc
    }
    if (is.null(Array_xct))
      stop("Problem with `plot_num` in `plot_maps(.)")
    if (tolower(Panel) == "category") {
      if (length(dim(Array_xct)) == 2)
        Nplot = 1
      if (length(dim(Array_xct)) == 3)
        Nplot = dim(Array_xct)[2]
      for (cI in 1:Nplot) {
        if (length(dim(Array_xct)) == 2) {
          Return = Mat_xt = Array_xct
          file_name = paste0(plot_code, ifelse(Nplot > 1, paste0("--", category_names[cI]), ""))
          plot_args = plot_variable_zz(Y_gt = Mat_xt[,], map_list = list(PlotDF = PlotDF, MapSizeRatio = MapSizeRatio), Year_Set = Year_Set, projargs = projargs, projargs2 = projargs2, working_dir = working_dir, file_name = file_name)
        } else {
          Return = Mat_xt = array(as.vector(Array_xct[, cI, ]), dim = dim(Array_xct)[c(1, 3)])
          file_name = paste0(plot_code, ifelse(Nplot > 1, paste0("--", category_names[cI]), ""))
          plot_args = plot_variable_zz(Y_gt = Mat_xt[, Years2Include, drop = FALSE], map_list = list(PlotDF = PlotDF, MapSizeRatio = MapSizeRatio), Year_Set = Year_Set, projargs = projargs, projargs2 = projargs2, working_dir = working_dir, panel_labels = Year_Set[Years2Include], file_name = file_name)
        }
      }
    }
    if (tolower(Panel) == "year") {
      Nplot = length(Years2Include)
      for (tI in 1:Nplot) {
        if (length(dim(Array_xct)) == 2)
          Mat_xc = Array_xct[, Years2Include[tI], drop = TRUE]
        if (length(dim(Array_xct)) == 3)
          Mat_xc = Array_xct[, , Years2Include[tI], drop = TRUE]
        Return = Mat_xc = array(as.vector(Mat_xc), dim = c(dim(Array_xct)[1], Ncategories))
        file_name = paste0(plot_code, ifelse(Nplot > 1, paste0("--", Year_Set[Years2Include][tI]), ""))
        plot_args = plot_variable(Y_gt = Mat_xc, map_list = list(PlotDF = PlotDF, MapSizeRatio = MapSizeRatio), projargs = projargs, working_dir = working_dir, panel_labels = category_names, file_name = file_name, n_cells = n_cells, ...)
      }
    }
  }
  if (is.null(Return) & quiet == FALSE)
    message(" # No plots selected in `plot_set`")
  return(invisible(Return))
}

plot_variable_zz<- function (Y_gt, map_list, Year_Set = Year_Set, panel_labels, projargs = "+proj=longlat", projargs2 = "+init=epsg:32619", map_resolution = "medium", file_name = "density", working_dir = paste0(OutFile, "/"), land_color = "grey") {
  if(FALSE){
    Y_gt = Mat_xt[, Years2Include, drop = FALSE]
    map_list = list(PlotDF = PlotDF, MapSizeRatio = MapSizeRatio)
    projargs = projargs
    working_dir = working_dir
    Year_Set = Year_Set
    file_name = file_name
    land_color = "grey"
  }

  if (is.vector(Y_gt)) {
    Y_gt = matrix(Y_gt, ncol = 1)
  }
  Y_gt = Y_gt[map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 0), "x2i"], , drop = FALSE]
  if (missing(panel_labels)) {
    panel_labels = rep("", ncol(Y_gt))
  }
  if (length(panel_labels) != ncol(Y_gt)) {
    warning("panel_labels and `ncol(Y_gt)` don't match: Changing panel_labels'")
    panel_labels = 1:ncol(Y_gt)
  }

  # Generating "nice" file names
  file_name.nice<- switch(file_name,
                          "epsilon_1" = "Spatio-temporal variation\nin encounter rate",
                          "epsilon_2" = "Spatio-temporal variation\nin positive catch rates",
                          "density" = "Catch (kg)",
                          "covariate_effects_1" = "Covariate effects\non encounter rates",
                          "covariate_effects_2" = "Covariate effects\non positive catch rates",
                          "omega_1" = "Spatial variation\nin encounter rate",
                          "omega_2" = "Spatial variation\nin positive catch rates")
  loc_g = map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 0), c("Lon", "Lat")]
  CRS_orig = sp::CRS("+proj=longlat")
  CRS_proj = sp::CRS(projargs)
  map_data = ne_countries(scale='medium', returnclass = 'sf') %>%
    filter(., admin %in% c("United States of America", "Canada"))
  map_data.proj = st_transform(map_data, crs = projargs2)

  rasts.out<- vector("list", length(ncol(Y_gt)))
  rasts.range<- Y_gt
  rast.lims<- c(round(min(rasts.range)-0.000001, 2), round(max(rasts.range) + 0.0000001, 2))

  if(ncol(Y_gt) == 1){
    data.df<- data.frame(loc_g, z = Y_gt[, tI])
    Points_LongLat = st_as_sf(data.df, coords = c("Lon", "Lat"), crs = CRS_orig)
    Points_proj = Points_LongLat %>%
      st_transform(., crs = projargs2)
    Points_bbox<- st_bbox(Points_proj)
    Raster_proj<- st_rasterize(Points_proj)

    plot.out<- ggplot() +
      geom_stars(data = Raster_proj, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis_c(name = file_name.nice, option = "viridis", na.value = "transparent", limits = rast.lims) +
      geom_sf(data = map_data.proj, fill = land_color, lwd = 0.2) +
      coord_sf(xlim = Points_bbox[c(1,3)], ylim = Points_bbox[c(2,4)], expand = FALSE, datum = sf::st_crs(projargs2)) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05))

    ggsave(filename = paste(working_dir, file_name, ".png", sep = ""), plot.out, width = 11, height = 8, units = "in")
  } else {
    for (tI in 1:ncol(Y_gt)) {
      data.df<- data.frame(loc_g, z = Y_gt[, tI])
      Points_LongLat = st_as_sf(data.df, coords = c("Lon", "Lat"), crs = CRS_orig)
      Points_proj = Points_LongLat %>%
        st_transform(., crs = projargs2)
      Points_bbox<- st_bbox(Points_proj)
      Raster_proj<- st_rasterize(Points_proj)

      rasts.out[[tI]]<- ggplot() +
        geom_stars(data = Raster_proj, aes(x = x, y = y, fill = z)) +
        scale_fill_viridis_c(name = file_name.nice, option = "viridis", na.value = "transparent", limits = rast.lims) +
        geom_sf(data = map_data.proj, fill = land_color, lwd = 0.2) +
        coord_sf(xlim = Points_bbox[c(1,3)], ylim = Points_bbox[c(2,4)], expand = FALSE, datum = sf::st_crs(projargs2)) +
        annotate("text", x = Points_bbox[3]-175000, y = Points_bbox[2] + 100000, label = parse_number(unique(as.character(year_labels))[tI])) +
        theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
    }

    all.plot<- wrap_plots(rasts.out, ncol = 10, nrow = 4, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
    ggsave(filename = paste(working_dir, file_name, ".png", sep = ""), all.plot, width = 11, height = 8, units = "in")
  }
  return(invisible(list(Par = Par, xlim = Points_bbox[c(1,3)], ylim = Points_bbox[c(2,4)])))
}

plot_residuals_zz<- function (Lat_i, Lon_i, TmbData, Report, Q, projargs = "+proj=longlat", working_dir = paste0(OutFile, "/"), spatial_list, extrapolation_list, Year_Set = NULL, Years2Include = NULL, zrange, ...) {
  if(FALSE){
    fit<- fit.basicfore
    settings = fit$settings
    plot_set = 3
    working_dir = OutFile
    year_labels = fit$year_labels
    years_to_plot = fit$years_to_plot
    use_biascorr = TRUE
    check_residuals = TRUE
    projargs = "+proj=longlat"
    n_samples = 100
    Lat_i = fit$data_frame[, "Lat_i"]
    Lon_i = fit$data_frame[, "Lon_i"]
    TmbData = fit$data_list
    Report = fit$Report
    Q = plot_quantile_diagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive", FileName_Phist = "Posterior_Predictive-Histogram", FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist", save_dir = working_dir)
    spatial_list = fit$spatial_list
    extrapolation_list = fit$extrapolation_list
    Year_Set = year_labels
    Years2Include = years_to_plot
    zrange = zrange
    legend_x = map_list[["Legend"]]$x/100
    legend_y = map_list[["Legend"]]$y/100
  }
  if (!("t_iz" %in% names(TmbData))) {
    TmbData$t_iz = matrix(TmbData$t_i, ncol = 1)
  }
  if (!("t_yz" %in% names(TmbData))) {
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol = 1)
  }
  exp_rate_xy = obs_rate_xy = total_num_xy = exp_num_xy = obs_num_xy = matrix(NA, nrow = spatial_list$n_x, ncol = nrow(TmbData$t_yz))
  for (yI in 1:nrow(TmbData$t_yz)) {
    which_i_in_y = (TmbData$t_iz == outer(rep(1, TmbData$n_i), TmbData$t_yz[yI, ]))
    which_i_in_y = which(apply(which_i_in_y, MARGIN = 1, FUN = all))
    if (length(which_i_in_y) > 0) {
      exp_rate_xy[, yI] = tapply(Report$R1_i[which_i_in_y], INDEX = factor(spatial_list$knot_i[which_i_in_y], levels = 1:spatial_list$n_x), FUN = mean)
      obs_rate_xy[, yI] = tapply(TmbData$b_i[which_i_in_y] > 0, INDEX = factor(spatial_list$knot_i[which_i_in_y], levels = 1:spatial_list$n_x), FUN = mean)
      total_num_xy[, yI] = tapply(TmbData$b_i[which_i_in_y], INDEX = factor(spatial_list$knot_i[which_i_in_y], levels = 1:spatial_list$n_x), FUN = length)
    }  else {
      total_num_xy[, yI] = 0
    }
    exp_num_xy = exp_rate_xy * total_num_xy
    obs_num_xy = obs_rate_xy * total_num_xy
  }
  Q1_xy = (obs_num_xy - exp_num_xy)/sqrt(exp_num_xy * (total_num_xy - exp_num_xy)/total_num_xy)
  which_pos = which(TmbData$b_i > 0)
  bvar_ipos = bpred_ipos = NULL
  if (all(c("var_y", "pred_y") %in% names(Q))) {
    bvar_ipos = Q[["var_y"]]
    bpred_ipos = Q[["pred_y"]]
  }
  if (all(c("var_y", "pred_y") %in% names(Q[[1]]))) {
    bvar_ipos = bpred_ipos = rep(NA, length = length(which_pos))
    for (i_e in 1:length(Q)) {
      which_pos_and_e = which(TmbData$e_i[which_pos] == (i_e - 1))
      bvar_ipos[which_pos_and_e] = Q[[i_e]][["var_y"]]
      bpred_ipos[which_pos_and_e] = Q[[i_e]][["pred_y"]]
    }
  }
  if (is.null(bvar_ipos) & is.null(bpred_ipos)) {
    stop("Something is wrong with `Q` input")
  }
  sum_obs_xy = sum_exp_xy = var_exp_xy = matrix(NA, nrow = spatial_list$n_x, ncol = nrow(TmbData$t_yz))
  for (yI in 1:nrow(TmbData$t_yz)) {
    which_i_in_y = (TmbData$t_iz == outer(rep(1, TmbData$n_i), TmbData$t_yz[yI, ]))
    which_i_in_y = which(apply(which_i_in_y, MARGIN = 1, FUN = all))
    which_i_in_y_and_pos = intersect(which_i_in_y, which_pos)
    which_ipos_in_y = (TmbData$t_iz[which_pos, ] == outer(rep(1, length(which_pos)), TmbData$t_yz[yI, ]))
    which_ipos_in_y = which(apply(which_ipos_in_y, MARGIN = 1, FUN = all))
    if (length(which_i_in_y_and_pos) > 0) {
      sum_obs_xy[, yI] = tapply(TmbData$b_i[which_i_in_y_and_pos], INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], levels = 1:spatial_list$n_x), FUN = sum)
      sum_exp_xy[, yI] = tapply(bpred_ipos[which_ipos_in_y], INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], levels = 1:spatial_list$n_x), FUN = sum)
      var_exp_xy[, yI] = tapply(bvar_ipos[which_ipos_in_y], INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], levels = 1:spatial_list$n_x), FUN = sum)
    }
  }
  Q2_xy = (sum_obs_xy - sum_exp_xy)/sqrt(var_exp_xy)
  if (!is.null(working_dir)) {
    for (zI in 1:2) {
      Q_xy = list(Q1_xy, Q2_xy)[[zI]]
      if (!missing(zrange)) {
        Q_xy = ifelse(Q_xy < zrange[1], zrange[1], Q_xy)
        Q_xy = ifelse(Q_xy > zrange[2], zrange[2], Q_xy)
        zlim = zrange
      } else {
        zlim = c(-1, 1) * ceiling(max(abs(Q_xy), na.rm = TRUE))
      }
      Col = colorRampPalette(colors = c("blue", "white", "red"))
      textmargin = "Pearson residual"
      plot_code = c("pearson_residuals_1", "pearson_residuals_2")[zI]
      x2i = spatial_list$NN_Extrap$nn.idx[, 1]

      if(all(extrapolation_list[["a_el"]] == 0)){
        Include = seq(from = 1, to = nrow(extrapolation_list[["a_el"]]))
      } else {
        Include = extrapolation_list[["Area_km2_x"]] > 0 & extrapolation_list[["a_el"]][,1] > 0
      }

      DF = cbind(extrapolation_list$Data_Extrap[, c("Lon", "Lat")], x2i = x2i, Include = Include)
      if (is.null(Year_Set))
        Year_Set = 1:ncol(Q_xy)
      if (is.null(Years2Include))
        Years2Include = 1:ncol(Q_xy)
      plot_args = plot_variable(Y_gt = ifelse(is.na(Q_xy), mean(zlim), Q_xy), map_list = list(PlotDF = DF), projargs = projargs, working_dir = working_dir, panel_labels = Year_Set[Years2Include], file_name = plot_code, zlim = zlim, col = Col, ...)
    }
  }
  Return = list(Q1_xy = Q1_xy, Q2_xy = Q2_xy)
  return(invisible(Return))
}

plot_biomass_index_zz<- function(TmbData, Sdreport, Year_Set = NULL, Years2Include = NULL, DirName = paste0(getwd(), "/"), PlotName = "Index", interval_width = 1, strata_names = NULL, category_names = NULL, use_biascorr = TRUE, plot_legend = TRUE, total_area_km2 = NULL, plot_log = FALSE, width = NULL, height = NULL, create_covariance_table = FALSE, ...)
{
  if (is.null(Sdreport))
    stop("Sdreport is NULL; please provide Sdreport")
  if (!is.null(category_names) && length(category_names) != TmbData$n_c)
    stop("`category_names` must have same length as `TmbData$n_c`")
  if (!is.null(Year_Set) && length(Year_Set) != TmbData$n_t)
    stop("`Year_Set` must have same length as `TmbData$n_t`")
  if (!is.null(strata_names) && length(strata_names) != TmbData$n_l)
    stop("`strata_names` must have same length as `TmbData$n_l`")
  if ("ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    ParName = "Index_tl"
    TmbData[["n_c"]] = 1
  }
  if ("ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    ParName = "Index_ctl"
  }
  if ("ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    ParName = "Index_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }
  if ("Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    ParName = "Index_tp"
    TmbData[["n_l"]] = 1
    TmbData[["n_c"]] = TmbData[["n_p"]]
  }
  if (!("t_iz" %in% names(TmbData))) {
    TmbData$t_iz = matrix(TmbData$t_i, ncol = 1)
  }
  if (!("t_yz" %in% names(TmbData))) {
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol = 1)
  }
  mfrow = c(ceiling(sqrt(TmbData$n_c)), ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c))))
  if (is.null(width))
    width = mfrow[2] * 3
  if (is.null(height))
    height = mfrow[1] * 3
  if (is.null(Year_Set))
    Year_Set = 1:TmbData$n_t
  if (is.null(Years2Include))
    Years2Include = 1:TmbData$n_t
  if (is.null(strata_names))
    strata_names = 1:TmbData$n_l
  if (is.null(category_names))
    category_names = 1:TmbData$n_c
  if ("unbiased" %in% names(Sdreport)) {
    if (all(is.na(Sdreport$unbiased$value))) {
      stop("You appear to be using bias-correction, but all values are NA. Please report problem to package author.")
    }
  }
  if ("treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options)) {
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  } else {
    treat_missing_as_zero = FALSE
  }
  SD = TMB::summary.sdreport(Sdreport)
  if (!"report" %in% names(as.list(args(TMB:::as.list.sdreport)))) {
    warning("package `TMB` should be updated to easily access standard errors")
  }
  SD_stderr = TMB:::as.list.sdreport(Sdreport, what = "Std. Error", report = TRUE)
  SD_estimate = TMB:::as.list.sdreport(Sdreport, what = "Estimate", report = TRUE)
  if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
    SD_estimate_biascorrect = TMB:::as.list.sdreport(Sdreport, what = "Est. (bias.correct)", report = TRUE)
  }
  if (any(is.na(SD_estimate)) | any(is.na(SD_stderr))) {
    stop("Problem: Standard errors contain NAs")
  }
  if (ParName %in% c("Index_tl", "Index_ctl", "Index_cyl")) {
    Index_ctl = log_Index_ctl = array(NA, dim = c(unlist(TmbData[c("n_c", "n_t", "n_l")]), 2), dimnames = list(category_names, Year_Set, strata_names, c("Estimate", "Std. Error")))
    if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
      Index_ctl[] = SD[which(rownames(SD) == ParName), c("Est. (bias.correct)", "Std. Error")]
    }
    if (!any(is.na(Index_ctl))) {
      message("Using bias-corrected estimates for abundance index (natural-scale)...")
    } else {
      message("Not using bias-corrected estimates for abundance index (natural-scale)...")
      Index_ctl[] = SD[which(rownames(SD) == ParName), c("Estimate", "Std. Error")]
    }
    if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
      log_Index_ctl[] = SD[which(rownames(SD) == paste0("ln_", ParName)), c("Est. (bias.correct)", "Std. Error")]
    }
    if (!any(is.na(log_Index_ctl))) {
      message("Using bias-corrected estimates for abundance index (log-scale)...")
    } else {
      message("Not using bias-corrected estimates for abundance index (log-scale)...")
      log_Index_ctl[] = SD[which(rownames(SD) == paste0("ln_", ParName)), c("Estimate", "Std. Error")]
    }
  }
  if (ParName %in% c("Index_tp")) {
    if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
      Index_ctl = aperm(array(c(Sdreport$unbiased$value[which(names(Sdreport$value) ==  ParName)], TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) ==  ParName), "Std. Error"]), dim = c(unlist(TmbData[c("n_t", "n_c", "n_l")]), 2), dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error"))), perm = c(2, 1, 3))
      if ("ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))) {
        log_Index_ctl = aperm(array(c(Sdreport$unbiased$value[which(names(Sdreport$value) == paste0("ln_", ParName))], TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) ==  paste0("ln_", ParName)), "Std. Error"]), dim = c(unlist(TmbData[c("n_t", "n_c", "n_l")]), 2), dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error"))), perm = c(2, 1, 3))
      } else {
        log_Index_ctl = log(Index_ctl)
        log_Index_ctl[, , , "Std. Error"] = log_Index_ctl[, , , "Std. Error"]/log_Index_ctl[, , , "Estimate"]
        warning("Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'")
      }
    } else {
      Index_ctl = aperm(array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == ParName), ], dim = c(unlist(TmbData[c("n_t", "n_c", "n_l")]), 2), dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error"))), perm = c(2, 1, 3, 4))
      if ("ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))) {
        log_Index_ctl = aperm(array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == paste0("ln_", ParName)), ], dim = c(unlist(TmbData[c("n_t", "n_c", "n_l")]), 2), dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error"))), perm = c(2, 1, 3, 4))
      } else {
        log_Index_ctl = log(Index_ctl)
        log_Index_ctl[, , , "Std. Error"] = log_Index_ctl[, , , "Std. Error"]/log_Index_ctl[, , , "Estimate"]
        warning("Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'")
      }
    }
  }
  if ("Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    Bratio_ctl = array(NA, dim = c(unlist(TmbData[c("n_c", "n_t", "n_l")]), 2), dimnames = list(category_names, Year_Set, strata_names, c("Estimate", "Std. Error")))
    if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
      Bratio_ctl[] = SD[which(rownames(SD) == "Bratio_cyl"), c("Est. (bias.correct)", "Std. Error")]
    }
    if (!any(is.na(Bratio_ctl))) {
      message("Using bias-corrected estimates for biomass ratio (natural-scale)...")
    } else {
      message("Not using bias-corrected estimates for biomass ratio (natural-scale)...")
      Bratio_ctl[] = SD[which(rownames(SD) == "Bratio_cyl"), c("Estimate", "Std. Error")]
    }
  } else {
    Bratio_ctl = NULL
  }
  if ("ln_Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    log_Bratio_ctl = array(NA, dim = c(unlist(TmbData[c("n_c", "n_t", "n_l")]), 2), dimnames = list(category_names, Year_Set, strata_names, c("Estimate", "Std. Error")))
    if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
      log_Bratio_ctl[] = SD[which(rownames(SD) == "ln_Bratio_cyl"), c("Est. (bias.correct)", "Std. Error")]
    }
    if (!any(is.na(log_Bratio_ctl))) {
      message("Using bias-corrected estimates for biomass ratio (log-scale)...")
    } else {
      message("Not using bias-corrected estimates for biomass ratio (log-scale)...")
      log_Bratio_ctl[] = SD[which(rownames(SD) == "ln_Bratio_cyl"), c("Estimate", "Std. Error")]
    }
  } else {
    log_Bratio_ctl = NULL
  }
  if ("Fratio_ct" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    Fratio_ct = array(NA, dim = c(unlist(TmbData[c("n_c", "n_t")]), 2), dimnames = list(category_names, Year_Set, c("Estimate", "Std. Error")))
    Fratio_ct[] = SD[which(rownames(SD) == "Fratio_ct"), c("Estimate", "Std. Error")]
  } else {
    Fratio_ct = NULL
  }
  if (!is.null(total_area_km2) & TmbData$n_c == 1) {
    message("Calculating naive design-based index -- do not use this, its intended only for comparison purposes")
    Calc_design = TRUE
    Design_t = tapply(TmbData$b_i/TmbData$a_i, INDEX = TmbData$t_i, FUN = mean) * total_area_km2/1000
    Design_t = cbind(Estimate = Design_t, `Std. Error` = sqrt(tapply(TmbData$b_i/TmbData$a_i, INDEX = TmbData$t_i, FUN = var)/tapply(TmbData$b_i/TmbData$a_i, INDEX = TmbData$t_i, FUN = length)) * total_area_km2/1000)
    Design_t = cbind(Design_t, CV = Design_t[, "Std. Error"]/Design_t[, "Estimate"])
  } else {
    Calc_design = FALSE
  }
  if (treat_missing_as_zero == TRUE) {
    Num_ct = tapply(TmbData$b_i, INDEX = list(factor(TmbData$c_i, levels = 1:TmbData$n_c - 1), factor(TmbData$t_i[, 1], levels = 1:TmbData$n_t - 1)), FUN = function(vec) {
      sum(!is.na(vec))
      })
    Num_ct = ifelse(is.na(Num_ct), 0, Num_ct)
    Index_ctl[, , , "Estimate"] = ifelse(Num_ct %o% rep(1, TmbData$n_l) == 0, 0, Index_ctl[, , , "Estimate"])
    Index_ctl[, , , "Std. Error"] = ifelse(Num_ct %o% rep(1, TmbData$n_l) == 0, NA, Index_ctl[, , , "Std. Error"])
    log_Index_ctl[, , , "Estimate"] = ifelse(Num_ct %o% rep(1, TmbData$n_l) == 0, -Inf, log_Index_ctl[, , , "Estimate"])
    log_Index_ctl[, , , "Std. Error"] = ifelse(Num_ct %o% rep(1, TmbData$n_l) == 0, NA, log_Index_ctl[, , ,  "Std. Error"])
  }
  Plot_suffix = "Biomass"
  if (!is.null(Bratio_ctl))
    Plot_suffix = c(Plot_suffix, "Bratio")
  for (plotI in 1:length(Plot_suffix)) {
    if (Plot_suffix[plotI] == "Biomass") {
      Array_ctl = Index_ctl
      log_Array_ctl = log_Index_ctl
    }
    if (Plot_suffix[plotI] == "Bratio") {
      Array_ctl = Bratio_ctl
      log_Array_ctl = log_Bratio_ctl
    }
    plot_index(Index_ctl = array(Index_ctl[, , , "Estimate"], dim(Index_ctl)[1:3]), sd_Index_ctl = array(log_Index_ctl[, , , "Std. Error"], dim(log_Index_ctl)[1:3]), Year_Set = Year_Set, Years2Include = Years2Include, strata_names = strata_names, category_names = category_names, DirName = DirName, PlotName = paste0(PlotName, "-", Plot_suffix[plotI],  ".png"), interval_width = interval_width, width = width, height = height, xlab = "Year", ylab = "Index", scale = "log", plot_args = list(log = ifelse(plot_log == TRUE, "y",  "")), Yrange = c(ifelse(plot_log == TRUE, NA, 0), NA))
  }
  if (!is.null(Fratio_ct)) {
    Array_ct = Fratio_ct
    Array_ct = ifelse(Array_ct == 0, NA, Array_ct)
    plot_index(Index_ctl = array(Array_ct[, , "Estimate"], dim(Array_ct)[1:2]), sd_Index_ctl = array(Array_ct[, , "Std. Error"], dim(Array_ct)[1:2]), Year_Set = Year_Set, Years2Include = Years2Include, strata_names = strata_names, category_names = category_names, DirName = DirName, PlotName = paste0(PlotName, "-Fratio.png"), scale = "uniform", interval_width = interval_width, width = width, height = height, xlab = "Year", ylab = "Fishing ratio")
  }
  if (!is.null(Bratio_ctl) & !is.null(Fratio_ct)) {
    Par = list(mar = c(2, 2, 1, 0), mgp = c(2, 0.5, 0), tck = -0.02, yaxs = "i", oma = c(1, 2, 0, 0), mfrow = mfrow, ...)
    Col = colorRampPalette(colors = c("blue", "purple", "red"))
    png(file = paste0(DirName, "/", PlotName, "-Status.png"), width = width, height = height, res = 200, units = "in")
    par(Par)
    Array1_ct = Bratio_ctl[, , 1, ]
    Array1_ct = ifelse(Array1_ct == 0, NA, Array1_ct)
    Array2_ct = Fratio_ct
    Array2_ct = ifelse(Array2_ct == 0, NA, Array2_ct)
    for (cI in 1:TmbData$n_c) {
      Xlim = c(0, max(1, Array1_ct[cI, Years2Include, "Estimate"] %o% c(1, 1) + Array1_ct[cI, Years2Include, "Std. Error"] %o%  c(-interval_width, interval_width), na.rm = TRUE))
      Ylim = c(0, max(2, Array2_ct[cI, Years2Include, "Estimate"] %o%  c(1, 1) + Array2_ct[cI, Years2Include, "Std. Error"] %o% c(-interval_width, interval_width), na.rm = TRUE))
      plot(1, type = "n", xlim = Xlim, ylim = Ylim, xlab = "", ylab = "", main = ifelse(TmbData$n_c > 1, category_names[cI], ""))
      points(x = Array1_ct[cI, Years2Include, "Estimate"], y = Array2_ct[cI, Years2Include, "Estimate"], col = Col(length(Year_Set))[Years2Include])
      for (tI in Years2Include) {
        lines(x = rep(Array1_ct[cI, tI, "Estimate"], 2), y = Array2_ct[cI, tI, "Estimate"] + Array2_ct[cI, tI, "Std. Error"] * c(-interval_width, interval_width), col = Col(length(Year_Set))[tI])
        lines(x = Array1_ct[cI, tI, "Estimate"] + Array1_ct[cI, tI, "Std. Error"] * c(-interval_width, interval_width), y = rep(Array2_ct[cI, tI, "Estimate"], 2), col = Col(length(Year_Set))[tI])
      }
      abline(v = 0.4, lty = "dotted")
      abline(h = 1, lty = "dotted")
    }
    legend("topright", bty = "n", fill = c(Col(length(Year_Set))[Years2Include[1]], Col(length(Year_Set))[rev(Years2Include)[1]]), legend = c(Year_Set[Years2Include[1]], Year_Set[rev(Years2Include)[1]]))
    mtext(side = 1:2, text = c("Biomass relative to unfished", "Fishing relative to F_40%"), outer = TRUE, line = c(0, 0))
    dev.off()
  }
  Table = NULL
  for (cI in 1:TmbData$n_c) {
    Tmp = data.frame(Year = Year_Set, Unit = 1, Fleet = rep(strata_names, each = TmbData$n_t), Estimate_metric_tons = as.vector(Index_ctl[cI, "Estimate"]), SD_log = as.vector(log_Index_ctl[cI, , , "Std. Error"]), SD_mt = as.vector(Index_ctl[cI, , , "Std. Error"]))
    if (TmbData$n_c > 1)
      Tmp = cbind(Category = category_names[cI], Tmp)
    Table = rbind(Table, Tmp)
  }
  if (!is.null(total_area_km2))
    Table = cbind(Table, `Naive_design-based_index` = Design_t)
  write.csv(Table, file = paste0(DirName, "/Table_for_SS3.csv"), row.names = FALSE)
  Return = list(Table = Table, log_Index_ctl = log_Index_ctl, Index_ctl = Index_ctl)
  if ("cov" %in% names(Sdreport) & create_covariance_table == TRUE) {
    DF = expand.grid(Category = 1:TmbData$n_c, Year = 1:TmbData$n_t, Stratum = 1:TmbData$n_l)
    Which = which(names(Sdreport$value) == ParName)
    Cov = Sdreport$cov[Which, Which]
    Corr = cov2cor(Cov) - diag(nrow(Cov))
    rowcolDF = cbind(RowNum = row(Corr)[lower.tri(Corr, diag = TRUE)],  ColNum = col(Corr)[lower.tri(Corr, diag = TRUE)])
    Table = cbind(DF[rowcolDF[, "ColNum"], ], DF[rowcolDF[, "RowNum"], ])
    colnames(Table) = paste0(colnames(Table), rep(c(1, 2), each = 3))
    Table = cbind(Table, Correlation = cov2cor(Cov)[lower.tri(Corr, diag = TRUE)], Covariance = Cov[lower.tri(Corr, diag = TRUE)])
    Table = cbind(Table, Index1 = Index_ctl[as.matrix(cbind(DF[rowcolDF[, "ColNum"], ], 1))], Index2 = Index_ctl[as.matrix(cbind(DF[rowcolDF[, "RowNum"], ], 1))])
    WhichZero = which((Table[, "Index1"] * Table[, "Index2"]) ==  0)
    Table[WhichZero, c("Correlation", "Covariance")] = 0
    Return = c(Return, Table_of_estimted_covariance = Table)
  }
  if (!is.null(Bratio_ctl))
    Return = c(Return, list(Bratio_ctl = Bratio_ctl))
  if (!is.null(log_Bratio_ctl))
    Return = c(Return, list(log_Bratio_ctl = log_Bratio_ctl))
  if (!is.null(Fratio_ct))
    Return = c(Return, list(Fratio_ct = Fratio_ct))
  return(invisible(Return))
}

# VAST fit_model edits ----------------------------------------------------
fit_model_zz<- function (settings,
                         Lat_i,
                         Lon_i,
                         t_iz,
                         b_i,
                         a_i,
                         c_iz = rep(0, length(b_i)),
                         v_i = rep(0, length(b_i)),
                         working_dir = paste0(getwd(), "/"),
                         Xconfig_zcp = NULL,
                         covariate_data,
                         formula = ~0,
                         design = "Season",
                         Q_ik = NULL,
                         newtonsteps = 1,
                         silent = TRUE, run_model = TRUE, test_fit = TRUE,
                         ...) {

  if(FALSE){
    "settings" = settings.use
    "design" = design
    "Lat_i" =  samp.dat.use[, "Lat"]
    "Lon_i" = samp.dat.use[,'Lon']
    "t_iz" = t_iz.aja.use
    "b_i" = samp.dat.use[,'Catch_KG']
    "a_i" = samp.dat.use[,'AreaSwept_km2']
    "v_i" = rep(0, nrow(samp.dat.use))
    "c_iz" = rep(0, nrow(samp.dat.use))
    "formula" = formula.use
    "covariate_data" = cov.dat.use
    "PredTF_i" = rep(0, nrow(samp.dat.use))
    working_dir = OutFile
    "Xconfig_zcp" = Xconfig_zcp.aja.use
    Q_ik = NULL
    "observations_LL" = cbind("Lat" = samp.dat.use[,'Lat'], "Lon" = samp.dat.use[, 'Lon'])
    "maximum_distance_from_sample" = maximum_distance_from_sample
    "grid_dim_km" = grid_dim_km
    newtonsteps = 1
    silent = TRUE
    run_model = FALSE
    test_fit = TRUE
    strata.limits = strata.limits


    "settings" = settings
    "design" = design
    "Lat_i" = vast.dat$SampleData[,'Lat']
    "Lon_i" = vast.dat$SampleData[,'Lon']
    "t_iz" = t_iz.use
    "b_i" = vast.dat$SampleData[,'Catch_KG']
    "c_iz" = rep(0, nrow(vast.dat$SampleData))
    "v_i" = rep(0, nrow(vast.dat$SampleData))
    "a_i" = vast.dat$SampleData[,'AreaSwept_km2']
    "formula" = formula.use
    "covariate_data" = vast.dat$CovariateData
    "PredTF_i" = vast.dat$SampleData[,'Pred_TF']
    Q_ik = NULL
    FieldConfig = FieldConfig.noomega2
    RhoConfig = RhoConfig.nobeta2
    Xconfig_zcp = Xconfig_zcp.use
    working_dir = OutFile
    run_model = TRUE
    newtonsteps = 1
    silent = TRUE
    run_model = FALSE
    test_fit = TRUE
    CompileDir = here::here()
  }

  extra_args = list(...)
  extra_args = c(extra_args, extra_args$extrapolation_args, extra_args$spatial_args, extra_args$optimize_args, extra_args$model_args)
  data_frame = data.frame(Lat_i = Lat_i, Lon_i = Lon_i, a_i = a_i, v_i = v_i, b_i = b_i, t_i = t_iz, c_iz = c_iz)

  # Time labels
  year_labels = unique(data_frame$t_i.Date) # Andrew edit
  years_to_plot = which(year_labels %in% data_frame$t_i.Date)

  # Set up and model fitting
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  capture.output(settings, file = file.path(working_dir, "settings.txt"))
  message("\n### Making extrapolation-grid")
  extrapolation_args_default = list(Region = settings$Region, strata.limits = NULL, zone = settings$zone)
  extrapolation_args_input = extra_args[intersect(names(extra_args), formalArgs(make_extrapolation_info))]
  extrapolation_args_input = combine_lists(input = extrapolation_args_input, default = extrapolation_args_default)
  extrapolation_list = do.call(what = make_extrapolation_info, args = extrapolation_args_input)
  message("\n### Making spatial information")
  spatial_args_default = list(grid_size_km = settings$grid_size_km, n_x = settings$n_x, Method = settings$Method, Lon_i = Lon_i, Lat_i = Lat_i, Extrapolation_List = extrapolation_list, DirPath = working_dir, Save_Results = TRUE, fine_scale = settings$fine_scale)
  spatial_args_input = extra_args[intersect(names(extra_args), formalArgs(make_spatial_info))]
  spatial_args_input = combine_lists(input = spatial_args_input, default = spatial_args_default)
  spatial_list = do.call(what = make_spatial_info, args = spatial_args_input)
  message("\n### Making data object")
  if (missing(covariate_data))
    covariate_data = NULL
  data_args_default = list(Version = settings$Version, FieldConfig = settings$FieldConfig, OverdispersionConfig = settings$OverdispersionConfig, RhoConfig = settings$RhoConfig, VamConfig = settings$VamConfig, ObsModel = settings$ObsModel, c_iz = c_iz, b_i = b_i, a_i = a_i, v_i = v_i, s_i = spatial_list$knot_i - 1, t_iz = t_iz, spatial_list = spatial_list, Options = settings$Options, Aniso = settings$use_anisotropy, Xconfig_zcp = Xconfig_zcp, covariate_data = covariate_data, formula = formula, Q_ik = Q_ik, design = design)
  data_args_input = combine_lists(input = extra_args, default = data_args_default)
  #data_list = do.call(what = make_data, args = data_args_input) # Error here... going into "make_covariates"
  data_list = do.call(what = make_data_zz, args = data_args_input)
  message("\n### Making TMB object")
  model_args_default = list(TmbData = data_list, RunDir = working_dir, Version = settings$Version, RhoConfig = settings$RhoConfig, loc_x = spatial_list$loc_x, Method = spatial_list$Method)
  model_args_input = extra_args[intersect(names(extra_args), formalArgs(make_model))]
  model_args_input = combine_lists(input = model_args_input, default = model_args_default)
  tmb_list = do.call(what = make_model, args = model_args_input)
  if (silent == TRUE)
    tmb_list$Obj$env$beSilent()
  if (run_model == FALSE) {
    Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings)
    class(Return) = "fit_model"
    return(Return)
  }
  message("\n### Estimating parameters")
  optimize_args_default1 = combine_lists(default = list(lower = tmb_list$Lower, upper = tmb_list$Upper, loopnum = 2), input = extra_args[intersect(names(extra_args), formalArgs(TMBhelper::fit_tmb))])
  optimize_args_input1 = list(obj = tmb_list$Obj, savedir = NULL, newtonsteps = 0, bias.correct = FALSE, control = list(eval.max = 10000, iter.max = 10000, trace = 1), quiet = TRUE, getsd = FALSE)
  optimize_args_input1 = combine_lists(default = optimize_args_default1, input = optimize_args_input1)
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, args = optimize_args_input1)
  if (exists("check_fit") & test_fit == TRUE) {
    problem_found = VAST::check_fit(parameter_estimates)
    if (problem_found == TRUE) {
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try\n",
           call. = FALSE)
    }
  }
  optimize_args_default2 = list(obj = tmb_list$Obj, lower = tmb_list$Lower, upper = tmb_list$Upper, savedir = working_dir, bias.correct = settings$bias.correct, newtonsteps = newtonsteps, bias.correct.control = list(sd = FALSE, split = NULL, nsplit = 1, vars_to_correct = settings$vars_to_correct), control = list(eval.max = 10000, iter.max = 10000, trace = 1), loopnum = 1)
  optimize_args_input2 = extra_args[intersect(names(extra_args), formalArgs(TMBhelper::fit_tmb))]
  optimize_args_input2 = combine_lists(input = optimize_args_input2, default = optimize_args_default2)
  optimize_args_input2 = combine_lists(input = list(startpar = parameter_estimates$par), default = optimize_args_input2)
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, args = optimize_args_input2)
  Report = tmb_list$Obj$report()
  ParHat = tmb_list$Obj$env$parList(parameter_estimates$par)
  input_args = list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, model_args_input = model_args_input, spatial_args_input = spatial_args_input, optimize_args_input1 = optimize_args_input1, optimize_args_input2 = optimize_args_input2)
  Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, parameter_estimates = parameter_estimates, Report = Report, ParHat = ParHat, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings, input_args = input_args)
  class(Return) = "fit_model"
  return(Return)
}

fit_tmb_zz<- function (obj, fn = obj$fn, gr = obj$gr, startpar = NULL, lower = -Inf, upper = Inf, getsd = TRUE, control = list(eval.max = 10000, iter.max = 10000, trace = 0), bias.correct = FALSE, bias.correct.control = list(sd = FALSE, split = NULL, nsplit = NULL, vars_to_correct = NULL), savedir = NULL, loopnum = 3, newtonsteps = 0, n = Inf, getReportCovariance = FALSE, getJointPrecision = FALSE, getHessian = FALSE, quiet = FALSE, ...)
{
  if (is.null(startpar))
    startpar = obj$par
  List = list(...)
  combine_lists = function(default, input) {
    output = default
    for (i in seq_along(input)) {
      if (names(input)[i] %in% names(default)) {
        output[[names(input)[i]]] = input[[i]]
      } else {
        output = c(output, input[i])
      }
    }
    return(output)
  }
  BS.control = list(sd = FALSE, split = NULL, nsplit = NULL, vars_to_correct = NULL)
  BS.control = combine_lists(default = BS.control, input = bias.correct.control)
  nlminb.control = list(eval.max = 10000, iter.max = 10000, trace = 0)
  nlminb.control = combine_lists(default = nlminb.control, input = control)
  start_time = Sys.time()
  parameter_estimates = nlminb(start = startpar, objective = fn, gradient = gr, control = nlminb.control, lower = lower, upper = upper)
  for (i in seq(2, loopnum, length = max(0, loopnum - 1))) {
    Temp = parameter_estimates[c("iterations", "evaluations")]
    parameter_estimates = nlminb(start = parameter_estimates$par, objective = fn, gradient = gr, control = nlminb.control, lower = lower, upper = upper)
    parameter_estimates[["iterations"]] = parameter_estimates[["iterations"]] +  Temp[["iterations"]]
    parameter_estimates[["evaluations"]] = parameter_estimates[["evaluations"]] + Temp[["evaluations"]]
  }
  for (i in seq_len(newtonsteps)) {
    g = as.numeric(gr(parameter_estimates$par))
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    parameter_estimates$par = parameter_estimates$par - solve(h, g)
    parameter_estimates$objective = fn(parameter_estimates$par)
  }
  parameter_estimates = parameter_estimates[c("par", "objective", "iterations", "evaluations")]
  parameter_estimates[["time_for_MLE"]] = Sys.time() - start_time
  parameter_estimates[["max_gradient"]] = max(abs(gr(parameter_estimates$par)))
  parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates[["max_gradient"]] <  1e-04, "There is no evidence that the model is not converged", "The model is likely not converged")
  parameter_estimates[["number_of_coefficients"]] = c(Total = length(unlist(obj$env$parameters)), Fixed = length(startpar), Random = length(unlist(obj$env$parameters)) - length(startpar))
  parameter_estimates[["AIC"]] = TMBhelper::TMBAIC(opt = parameter_estimates)
  if (n != Inf) {
    parameter_estimates[["AICc"]] = TMBhelper::TMBAIC(opt = parameter_estimates, n = n)
    parameter_estimates[["BIC"]] = TMBhelper::TMBAIC(opt = parameter_estimates, p = log(n))
  }
  parameter_estimates[["diagnostics"]] = data.frame(Param = names(startpar), starting_value = startpar, Lower = lower, MLE = parameter_estimates$par, Upper = upper, final_gradient = as.vector(gr(parameter_estimates$par)))
  if (getsd == TRUE) {
    sd_time = Sys.time()
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    if (is.character(try(chol(h), silent = TRUE))) {
      warning("Hessian is not positive definite, so standard errors are not available")
      if (!is.null(savedir)) {
        capture.output(parameter_estimates, file = file.path(savedir, "parameter_estimates.txt"))
      }
      return(list(opt = parameter_estimates, h = h))
    }
    if (bias.correct == FALSE | is.null(BS.control[["vars_to_correct"]])) {
      if (!is.null(BS.control[["nsplit"]])) {
        if (BS.control[["nsplit"]] == 1)
          BS.control[["nsplit"]] = NULL
      }
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h, bias.correct = bias.correct, bias.correct.control = BS.control[c("sd", "split", "nsplit")], getReportCovariance = getReportCovariance, getJointPrecision = getJointPrecision, ...)
    } else {
      if ("ADreportIndex" %in% names(obj$env)) {
        Which = as.vector(unlist(obj$env$ADreportIndex()[BS.control[["vars_to_correct"]]]))
      } else {
        parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h, bias.correct = FALSE, getReportCovariance = FALSE, getJointPrecision = FALSE, ...)
        Which = which(rownames(summary(parameter_estimates[["SD"]], "report")) %in% BS.control[["vars_to_correct"]])
      }
      if (!is.null(BS.control[["nsplit"]]) && BS.control[["nsplit"]] >  1) {
        Which = split(Which, cut(seq_along(Which), BS.control[["nsplit"]]))
      }
      Which = Which[sapply(Which, FUN = length) > 0]
      if (length(Which) == 0)
        Which = NULL
      message(paste0("Bias correcting ", length(Which), " derived quantities"))
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, par.fixed = parameter_estimates$par, hessian.fixed = h, bias.correct = TRUE, bias.correct.control = list(sd = BS.control[["sd"]], split = Which, nsplit = NULL), getReportCovariance = getReportCovariance, getJointPrecision = getJointPrecision, ...)
    }
    parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates$SD$pdHess == TRUE, parameter_estimates[["Convergence_check"]], "The model is definitely not converged")
    parameter_estimates[["time_for_sdreport"]] = Sys.time() - sd_time
    if (getHessian == TRUE) {
      parameter_estimates[["hessian"]] = h
    }
  }
  parameter_estimates[["time_for_run"]] = Sys.time() - start_time
  if (!is.null(savedir)) {
    save(parameter_estimates, file = file.path(savedir, "parameter_estimates.RData"))
    capture.output(parameter_estimates, file = file.path(savedir, "parameter_estimates.txt"))
  }
  if (quiet == FALSE & parameter_estimates[["Convergence_check"]] !=
      "There is no evidence that the model is not converged") {
    message("#########################")
    message(parameter_estimates[["Convergence_check"]])
    message("#########################")
  }
  return(parameter_estimates)
}


#####
## Andrew's VAST function edits
#####

# “Plot” function edits ---------------------------------------------------
plot_results_zz<- function (fit, settings = fit$settings, plot_set = 3, working_dir = paste0(OutFile, "/"), year_labels = fit$year_labels, years_to_plot = fit$years_to_plot, use_biascorr = TRUE, map_list, category_names, check_residuals = TRUE, projargs = "+proj=longlat", projargs2 = "+init=epsg:32619", zrange, n_samples = 100, ...) {
  
  if(FALSE){
    fit<- t1
    settings = fit$settings
    plot_set = c(3, 6, 7, 13, 14, 15, 16)
    working_dir = paste(shared.path(os.use = os_use, group = "Mills Lab", folder = "Projects"), "ForecastingChallenge/Temp Results/", sep = "/")
    year_labels = fit$year_labels
    years_to_plot = fit$years_to_plot
    use_biascorr =FALSE
    check_residuals = TRUE
    projargs = "+proj=longlat"
    projargs2 = "+init=epsg:32619"
    n_samples = 100
    category_names = "Atlantic cod"
    
    ... = NULL
  }
  
  if (is.null(fit$Report))
    stop("`fit$Report` is missing, please check inputs")
  if (missing(category_names))
    category_names = 1:fit$data_list$n_c
  message("\n### Making plots of data availability and knots")
  plot_data_zz(Extrapolation_List = fit$extrapolation_list, Spatial_List = fit$spatial_list, Lat_i = fit$data_frame[, "Lat_i"], Lon_i = fit$data_frame[, "Lon_i"], Year_i = fit$data_frame[, "t_i"], projargs = "+proj=longlat",  projargs2 = "+init=epsg:32619", PlotDir = working_dir, Year_Set = year_labels)
  if (missing(map_list)) {
    message("\n### Obtaining default settings for plotting maps")
    map_list = make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
  }
  message("\n### Making plot of encounter probability")
  Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = cbind(Catch_KG = fit$data_frame[, "b_i"]), DirName = working_dir)
  message("\n### Making plot of anisotropy")
  plot_anisotropy(FileName = paste0(working_dir, "Aniso.png"), Obj = fit, Report = fit$Report, TmbData = fit$data_list)
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
    Year_i = fit$data_frame[, "t_i"]
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

var_maps_zz<- function (Scenario, Season, Fit, Variable, Species, out_folder, ...) {
  if(FALSE){
   Scenario = pred_temp$Scenario[[1]]
   Season = pred_temp$Season[[1]]
   Fit = t1
   Variable = "Density"
   Species = c("Atlantic cod")
   out_folder = res_folder_use
  }
  
  # Get what we need from fitted object
  TmbData<- Fit$data_list
  Report<- Fit$Report
  Sdreport<- Fit$parameter_estimates$SD
  Year_Set = Fit$year_labels
  Years2Include = Fit$years_to_plot
  category_names = Species
  
  if ("D_gct" %in% names(Report)) {
    if (is.null(Year_Set))
      Year_Set = 1:dim(Report$D_gct)[3]
    if (is.null(Years2Include))
      Years2Include = 1:dim(Report$D_gct)[3]
    if (is.null(category_names))
      category_names = 1:dim(Report$D_gct)[2]
    Ncategories = dim(Report$D_gct)[2]
    Nyears = dim(Report$D_gct)[3]
  }
  if (Nyears != length(Year_Set)) {
    stop("Problem with `Year_Set`")
  }
  if (Ncategories != length(category_names)) {
    stop("Problem with `category_names`")
  }
  
  # Land info
  map_data = ne_countries(scale='medium', returnclass = 'sf') %>%
    filter(., admin %in% c("United States of America", "Canada"))
  xlim_use<- c(-77, -65)
  ylim_use<- c(35, 45)

  # Raster template
  # First get locations of extrapolation points
  dat_loc<- data.frame(Fit$spatial_list$latlon_g) %>%
    mutate(., "PtID" = seq(from = 1, to = nrow(.)))
  
  # Now, need an extent object from this...
  dat_loc_sf<- dat_loc %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) %>%
    st_transform(., crs = Fit$extrapolation_list$projargs)
  dat_loc_bbox<- st_bbox(dat_loc_sf)
  dat_rast_temp <- raster(xmn = dat_loc_bbox['xmin'], xmx = dat_loc_bbox['xmax'], ymx = dat_loc_bbox['ymax'], ymn = dat_loc_bbox['ymin'], res = Fit$settings$grid_size_km, crs = Fit$extrapolation_list$projargs)
  values(dat_rast_temp)<- NA
  
  for(i in seq_along(category_names)){
    spp_use<- category_names[i]
    
    dat_loc<- data.frame(Fit$spatial_list$latlon_g) %>%
      mutate(., "PtID" = seq(from = 1, to = nrow(.)))
    
    if(Variable == "Density"){
      dat_array<- log(Report$D_gct[,i,]+1)
      colnames(dat_array)<- Year_Set
      dat_df<- dat_array %>%
        as_tibble() %>%
        mutate(., "CommonName" = rep(spp_use, nrow(.)),
               "Variable" = rep(Variable, nrow(.)),
               "PtID" = seq(from = 1, to = nrow(.))) %>%
        gather(., "Year", "Value", -CommonName, -Variable, -PtID)
      
      dat_out<- dat_df %>%
        left_join(., dat_loc_sf, by = c("PtID" = "PtID")) %>%
        st_as_sf()
      
      rasts_out<- vector("list", length = length(unique(dat_out$Year)))
      rast_lims<- c(0, max(dat_array, na.rm = TRUE))
      
      # Make the plots
      for(j in seq_along(unique(dat_out$Year))){
        dat_plot<- dat_out %>%
          filter(., Year == unique(dat_out$Year)[j])
        dat_rast<- resample(rasterize(dat_plot, dat_rast_temp, field = "Value"), dat_rast_temp, na.rm = TRUE)
        # Re project back to regular lat/lon
        dat_rast_plot<- as.data.frame(projectRaster(dat_rast, crs = 4326), xy = TRUE)

        rasts_out[[j]]<- ggplot() +
          geom_tile(data = dat_rast_plot, aes(x = x, y = y, fill = layer)) +
          scale_fill_viridis_c(name = Variable, option = "viridis", na.value = "transparent", limits = rast_lims) +
          geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
          coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
          annotate("text", x = -68, y = 36, label = unique(dat_out$Year)[j]) +
          theme(panel.background = element_rect(fill = "white"), 
                panel.border = element_rect(fill = NA), 
                axis.text.x=element_blank(), 
                axis.text.y=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title = element_blank())
      }
      all_plot<- wrap_plots(rasts_out, ncol = 10, nrow = 4, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste(out_folder, "/", Scenario, "_", Variable, ".jpg", sep = ""), all_plot, width = 11, height = 8, units = "in", scale = 1.25)
      # Return the data...
      return(rasts_out)
    }
    
    if(Variable == "Density_CV"){
     
    }
    
    if(Variable == "Epsilon_1"){
      dat_array<- Report$Epsilon1_gct[,i,]
      colnames(dat_array)<- Year_Set
      dat_df<- dat_array %>%
        as_tibble() %>%
        mutate(., "CommonName" = rep(spp_use, nrow(.)),
               "Variable" = rep({{Variable}}, nrow(.)),
               "PtID" = seq(from = 1, to = nrow(.))) %>%
        gather(., "Year", "Value", -CommonName, -Variable, -PtID)
      
      dat_out<- dat_df %>%
        left_join(., dat_loc_sf, by = c("PtID" = "PtID")) %>%
        st_as_sf()
      
      rasts_out<- vector("list", length = length(unique(dat_out$Year)))
      rast_lims<- c(min(dat_array, na.rm = TRUE), max(dat_array, na.rm = TRUE))
      
      # Make the plots
      for(j in seq_along(unique(dat_out$Year))){
        dat_plot<- dat_out %>%
          filter(., Year == unique(dat_out$Year)[j])
        dat_rast<- resample(rasterize(dat_plot, dat_rast_temp, field = "Value"), dat_rast_temp, na.rm = TRUE)
        # Re project back to regular lat/lon
        dat_rast_plot<- as.data.frame(projectRaster(dat_rast, crs = 4326), xy = TRUE)
        
        rasts_out[[j]]<- ggplot() +
          geom_tile(data = dat_rast_plot, aes(x = x, y = y, fill = layer)) +
          scale_fill_gradient2(name = {{Variable}}, high = "#d73027", low = "#4575b4", midpoint = 0, mid = "#ffffbf", na.value = "transparent", limits = rast_lims) +
          geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
          coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
          annotate("text", x = -68, y = 36, label = unique(dat_out$Year)[j]) +
          theme(panel.background = element_rect(fill = "white"), 
                panel.border = element_rect(fill = NA), 
                axis.text.x=element_blank(), 
                axis.text.y=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title = element_blank())
      }
      all_plot<- wrap_plots(rasts_out, ncol = 10, nrow = 4, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste(out_folder, "/", Scenario, "_", Variable, ".jpg", sep = ""), all_plot, width = 11, height = 8, units = "in", scale = 1.25)
      # Return the data...
      return(rasts_out)
    }
    
    if(Variable == "Epsilon_2"){
      dat_array<- Report$Epsilon2_gct[,i,]
      colnames(dat_array)<- Year_Set
      dat_df<- dat_array %>%
        as_tibble() %>%
        mutate(., "CommonName" = rep(spp_use, nrow(.)),
               "Variable" = rep(Variable, nrow(.)),
               "PtID" = seq(from = 1, to = nrow(.))) %>%
        gather(., "Year", "Value", -CommonName, -Variable, -PtID)
      
      dat_out<- dat_df %>%
        left_join(., dat_loc_sf, by = c("PtID" = "PtID")) %>%
        st_as_sf()
      
      rasts_out<- vector("list", length = length(unique(dat_out$Year)))
      rast_lims<- c(min(dat_array, na.rm = TRUE), max(dat_array, na.rm = TRUE))
      
      # Make the plots
      for(j in seq_along(unique(dat_out$Year))){
        dat_plot<- dat_out %>%
          filter(., Year == unique(dat_out$Year)[j])
        dat_rast<- resample(rasterize(dat_plot, dat_rast_temp, field = "Value"), dat_rast_temp, na.rm = TRUE)
        # Re project back to regular lat/lon
        dat_rast_plot<- as.data.frame(projectRaster(dat_rast, crs = 4326), xy = TRUE)
        
        rasts_out[[j]]<- ggplot() +
          geom_tile(data = dat_rast_plot, aes(x = x, y = y, fill = layer)) +
          scale_fill_gradient2(name = Variable, high = "#d73027", low = "#4575b4", midpoint = 0, mid = "#ffffbf", na.value = "transparent", limits = rast_lims) +
          geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
          coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
          annotate("text", x = -68, y = 36, label = unique(dat_out$Year)[j]) +
          theme(panel.background = element_rect(fill = "white"), 
                panel.border = element_rect(fill = NA), 
                axis.text.x=element_blank(), 
                axis.text.y=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title = element_blank())
      }
      all_plot<- wrap_plots(rasts_out, ncol = 10, nrow = 4, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste(out_folder, "/", Scenario, "_", Variable, ".jpg", sep = ""), all_plot, width = 11, height = 8, units = "in", scale = 1.25)
      # Return the data...
      return(rasts_out)
    }
    
    if(Variable == "Omega_1"){
      dat_array<- Report$Omega1_gc[,]
      dat_df<- dat_array %>%
        as_tibble() %>%
        mutate(., "CommonName" = rep(spp_use, nrow(.)),
               "Variable" = rep(Variable, nrow(.)),
               "PtID" = seq(from = 1, to = nrow(.)))
      names(dat_df)[1]<- "Value"
      
      dat_out<- dat_df %>%
        left_join(., dat_loc_sf, by = c("PtID" = "PtID")) %>%
        st_as_sf()
      
      rast_lims<- c(min(dat_array, na.rm = TRUE), max(dat_array, na.rm = TRUE))
      
      # Make the plots
      dat_plot<- dat_out 
      dat_rast<- resample(rasterize(dat_plot, dat_rast_temp, field = "Value"), dat_rast_temp, na.rm = TRUE)
      # Re project back to regular lat/lon
      dat_rast_plot<- as.data.frame(projectRaster(dat_rast, crs = 4326), xy = TRUE)
      
      rast_out<- ggplot() +
        geom_tile(data = dat_rast_plot, aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(name = Variable, high = "#d73027", low = "#4575b4", midpoint = 0, mid = "#ffffbf", na.value = "transparent", limits = rast_lims) +
        geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
        coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
        annotate("text", x = -68, y = 36, label = Variable) +
        theme(panel.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill = NA), 
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title = element_blank())
      ggsave(filename = paste(out_folder, "/", Scenario, "_", Variable, ".jpg", sep = ""), rast_out, width = 11, height = 8, units = "in", scale = 1.25)
      # Return the data...
      return(rast_out)
    }
    
    if(Variable == "Omega_2"){
      dat_array<- Report$Omega2_gc[,]
      dat_df<- dat_array %>%
        as_tibble() %>%
        mutate(., "CommonName" = rep(spp_use, nrow(.)),
               "Variable" = rep(Variable, nrow(.)),
               "PtID" = seq(from = 1, to = nrow(.)))
      names(dat_df)[1]<- "Value"
      
      dat_out<- dat_df %>%
        left_join(., dat_loc_sf, by = c("PtID" = "PtID")) %>%
        st_as_sf()
      
      rast_lims<- c(min(dat_array, na.rm = TRUE), max(dat_array, na.rm = TRUE))
      
      # Make the plots
      dat_plot<- dat_out 
      dat_rast<- resample(rasterize(dat_plot, dat_rast_temp, field = "Value"), dat_rast_temp, na.rm = TRUE)
      # Re project back to regular lat/lon
      dat_rast_plot<- as.data.frame(projectRaster(dat_rast, crs = 4326), xy = TRUE)
      
      rast_out<- ggplot() +
        geom_tile(data = dat_rast_plot, aes(x = x, y = y, fill = layer)) +
        scale_fill_gradient2(name = Variable, high = "#d73027", low = "#4575b4", midpoint = 0, mid = "#ffffbf", na.value = "transparent", limits = rast_lims) +
        geom_sf(data = map_data, fill = land_color, lwd = 0.2) +
        coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
        annotate("text", x = -68, y = 36, label = Variable) +
        theme(panel.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill = NA), 
              axis.text.x=element_blank(), 
              axis.text.y=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title = element_blank())
      ggsave(filename = paste(out_folder, "/", Scenario, "_", Variable, ".jpg", sep = ""), rast_out, width = 11, height = 8, units = "in", scale = 1.25)
      # Return the data...
      return(rast_out)
    }
  }
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

biomass_index_zz<- function(Scenario, Season, Species, Fit, ...)
{
  if(FALSE){
    Scenario = fore_challenge_use
    Season = ifelse(grepl("SPRING", res_use), "Spring", "Fall")
    DirName = working_dir
    TmbData = fit$data_list
    Sdreport = fit$parameter_estimates$SD
    Year_Set = year_labels
    Years2Include = years_to_plot
    use_biascorr = use_biascorr
    category_names = category_names
    strata_names = colnames(fit$extrapolation_list$a_el)
    total_area_km2 = NULL
  }
  
  # Get what we need from fitted object
  TmbData<- Fit$data_list
  Sdreport<- Fit$parameter_estimates$SD
  Year_Set = Fit$year_labels
  Years2Include = Fit$years_to_plot
  use_biascorr = Fit$settings$bias.correct
  category_names = Species
  strata_names = colnames(Fit$extrapolation_list$a_el)
  total_area_km2 = NULL
  
  if (is.null(Sdreport))
    stop("Sdreport is NULL; please provide Sdreport")
  if (!is.null(category_names) && length(category_names) != TmbData$n_c)
    stop("`category_names` must have same length as `TmbData$n_c`")
  if (!is.null(Year_Set) && length(Year_Set) != TmbData$n_t)
    stop("`Year_Set` must have same length as `TmbData$n_t`")
  if (!is.null(strata_names) && length(strata_names) != TmbData$n_l)
    stop("`strata_names` must have same length as `TmbData$n_l`")
  if ("ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
    ParName = "Index_ctl"
  }
  if (!("t_iz" %in% names(TmbData))) {
    TmbData$t_iz = matrix(TmbData$t_i, ncol = 1)
  }
  if (!("t_yz" %in% names(TmbData))) {
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol = 1)
  }
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
  Array_ctl = Index_ctl
  names(dim(Array_ctl))<- c("Category", "Time", "Strata", "EstimateAndSD")
  log_Array_ctl = log_Index_ctl
  names(dim(log_Array_ctl))<- c("Category", "Time", "Strata", "EstimateAndSD")
  
  # Get index table...
  for (catI in 1:length(dim(Array_ctl)[1])) {
    temp_df_mean<- data.frame(Array_ctl[catI,,,1])
    colnames(temp_df_mean)<- strata_names
    temp_df_mean_l<- temp_df_mean %>%
      mutate(., "Year" = Year_Set) %>%
      gather(., "Strata", "Mean", -Year)
    
    temp_df_se<- data.frame(Array_ctl[catI,,,2])
    colnames(temp_df_se)<- strata_names
    temp_df_se_l<- temp_df_se %>%
      mutate(., "Year" = Year_Set) %>%
      gather(., "Strata", "SE", -Year)
    
    temp_df_out<- temp_df_mean_l %>%
      left_join(., temp_df_se_l) %>%
      mutate(., "Category" = rep(category_names[catI], nrow(.)),
             "Index" = rep("Biomass", nrow(.)),
             "Season" = rep(Season, nrow(.)),
             "Scenario" = rep(Scenario, nrow(.)))
    
    # Log array
    temp_df_log_mean<- data.frame(log_Array_ctl[catI,,,1])
    colnames(temp_df_log_mean)<- strata_names
    temp_df_log_mean_l<- temp_df_log_mean %>%
      mutate(., "Year" = Year_Set) %>%
      gather(., "Strata", "Mean", -Year)
    
    temp_df_log_se<- data.frame(log_Array_ctl[catI,,,2])
    colnames(temp_df_log_se)<- strata_names
    temp_df_log_se_l<- temp_df_log_se %>%
      mutate(., "Year" = Year_Set) %>%
      gather(., "Strata", "SE", -Year)
    
    temp_df_log_out<- temp_df_log_mean_l %>%
      left_join(., temp_df_log_se_l) %>%
      mutate(., "Category" = rep(category_names[catI], nrow(.)),
             "Index" = rep("Log_Biomass", nrow(.)),
             "Season" = rep(Season, nrow(.)),
             "Scenario" = rep(Scenario, nrow(.)))
    
    temp_df_out<- bind_rows(temp_df_out, temp_df_log_out)
    
    if(catI == 1){
      all_index_out<- temp_df_out
    } else {
      all_index_out<- bind_rows(all_index_out)
    }
  }
  
  # Return it
  return(all_index_out)
}


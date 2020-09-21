# Load packages
library(VAST)
library(splines)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )

# Define formula.  In this case I'm demonstrating how to use a basis-spline with
# three degrees of freedom to model a nonlinear effect of log-transformed bottom depth,
# based on example developed by Nicholas Ducharme-Barth.
formula = ~ bs( log(BOT_DEPTH), knots=3, intercept=FALSE)

# set Year = NA to treat all covariates as "static" (not changing among years)
# If using a mix of static and dynamic covariates, please email package author to add easy capability
example$covariate_data[,'Year'] = NA

# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 formula = formula,
                 covariate_data = example$covariate_data )


# Adjusting settings
fieldconfig_all1<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rhoconfig_intrw_epsar1<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 4, "Epsilon2" = 4)

fit = try(capture.output(fit_model( "settings" = settings,
                 "FieldConfig" = fieldconfig_all1,
                 "RhoConfig" = rhoconfig_intrw_epsar1,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 formula = formula,
                 covariate_data = example$covariate_data)))


fit = purrr::quietly(try(fit_model( "settings" = settings,
                                    "FieldConfig" = fieldconfig_all1,
                                    "RhoConfig" = rhoconfig_intrw_epsar1,
                                    Lat_i = example$sampling_data[,'Lat'],
                                    Lon_i = example$sampling_data[,'Lon'],
                                    t_i = example$sampling_data[,'Year'],
                                    b_i = example$sampling_data[,'Catch_KG'],
                                    a_i = example$sampling_data[,'AreaSwept_km2'],
                                    formula = formula,
                                    covariate_data = example$covariate_data)))
unlist(attributes(fit))
attr(fit, "try-error")
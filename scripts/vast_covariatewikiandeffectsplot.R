###### DEVELOPMENT VERSION OF VAST SHOULD BE USED
library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )

# Define formula.
# In this case I'm demonstrating how to use a basis-spline with
# three degrees of freedom to model a nonlinear effect of log-transformed bottom depth,
# based on example developed by Nicholas Ducharme-Barth.
X1_formula = ~ bs( log(BOT_DEPTH), knots=3, intercept=FALSE)
X2_formula = ~ bs( log(BOT_DEPTH), knots=3, intercept=FALSE)

# If all covariates as "static" (not changing among years),
#  then set Year = NA to cause values to be duplicated internally for all values of Year
# If using a mix of static and dynamic covariates,
#  then duplicate rows for static covariates for every value of Year
# Here, all covariates are static, so I'm using the first approach.
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
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data )


library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X2", 
                         xlevels = 100)
plot(pred)


# With mgcv???
gam_df<- data.frame("Response" = example$sampling_data[,'Catch_KG'], "BOT_DEPTH" = example$covariate_data[,'BOT_DEPTH'])
gam_fit<- gam(Response ~ s(BOT_DEPTH, bs = "bs", m = c(3,2)) - 1, data = gam_df)
plot(gam_fit)

plot(pred)


# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial
Partial = partial( object = fit,
                   pred.var = "CPE",
                   pred.fun = pred.fun,
                   train = fit$covariate_data )

Partial = partial( object = fit,
                   pred.var = "BOT_DEPTH",
                   pred.fun = pred.fun,
                   train = fit$covariate_data )

# Make plot using ggplot2
library(ggplot2)
autoplot(Partial)

# What about with an "old" model -- remove the effects
fit.old<- fit
fit.old$effects<- NULL
names(fit.old)

covariate_data_full<- fit.old$covariate_data
catchability_data_full<- fit.old$catchability_data


# Plot 1st linear predictor
dev.off()
pred = Effect.fit_model( fit.old,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1" )

# Error with prior weights...
# Now trying adjusted function that creates and adds to "effects" slot of model fitted object, as done in newer FishStatsUtils::fit_model functions...
new.func.path<- "~/GitHub/ForecastingChallenge/scripts/effects_aja.R"
source(new.func.path)
pred = Effect.fit_model_test( fit.old,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1" )

# Still an error, but a different one. After trying all I could think of, only other idea was seeing if effects need to be added first??
fit.old$effects
if(is.null(fit.old$effects)){
  fit.old$effects<- list()
  if(!is.null(fit.old$catchability_data)) {
    catchability_data_full = data.frame(fit.old$catchability_data, 
                                        linear_predictor = 0)
    Q1_formula_full = update.formula(fit.old$Q1_formula, linear_predictor ~ 
                                       . + 0)
    call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
    Q2_formula_full = update.formula(fit.old$Q2_formula, linear_predictor ~ 
                                       . + 0)
    call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
    fit.old$effects = c(Return$effects, list(call_Q1 = call_Q1, 
                                         call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
  }
  if (!is.null(fit.old$covariate_data)) {
    covariate_data_full = data.frame(fit.old$covariate_data, linear_predictor = 0)
    X1_formula_full = update.formula(fit.old$X1_formula, linear_predictor ~ 
                                       . + 0)
    call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
    X2_formula_full = update.formula(fit.old$X2_formula, linear_predictor ~ 
                                       . + 0)
    call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
    fit.old$effects = c(fit.old$effects, list(call_X1 = call_X1, 
                                      call_X2 = call_X2, covariate_data_full = covariate_data_full))
  }
}
pred = Effect.fit_model_test( fit.old,
                              focal.predictors = c("BOT_DEPTH"),
                              which_formula = "X1" )
plot(pred)
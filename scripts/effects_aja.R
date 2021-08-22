#' Calculate effects for plotting
#'
#' @title Adapts package \code{effects}
#'
#' @rawNamespace S3method(effects::Effect, fit_model)
#' @export
Effect.fit_model_test = function (focal.predictors, mod, which_formula="X1", ...) {
  
  if(FALSE){
    # Model with effects slot...
    mod = fit
    # Define formula.
    X1_formula = ~ bs( log(BOT_DEPTH), knots=10, intercept=FALSE) 
    X2_formula = ~ bs( log(BOT_DEPTH), knots=10, intercept=FALSE)
    # Must add data-frames to global environment (hope to fix in future)
    covariate_data_full<- mod$effects$covariate_data_full
    catchability_data_full<- mod$effects$catchability_data_full
    focal.predictors = c("BOT_DEPTH")
    which_formula = "X1"
    
    # Model without effects slot, fit with previous VAST version
    mod = t1
    # Define formula.
    X1_formula = mod$X1_formula 
    X2_formula = mod$X2_formula
    focal.predictors = c("AVGDEPTH")
    which_formula = "X1"
  }
  
  # Error checks
  if( mod$data_list$n_c>1 ){
    stop("`Effect.fit_model` is not currently designed for multivariate models")
  }
  if( !all(c("covariate_data_full","catchability_data_full") %in% ls(.GlobalEnv)) ){
    stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
  }
  if( !requireNamespace("effects") ){
    stop("please install the effects package")
  }
  
  # Check for effects slot in mod, if it isn't there, try to make it -- seems to happen if using older version of FishStatsUtils::fit_model
  if(is.null(mod$effects)){
    mod$effects<- list()
    if(!is.null(mod$catchability_data)) {
      catchability_data_full = data.frame(mod$catchability_data, 
                                          linear_predictor = 0)
      Q1_formula_full = update.formula(mod$Q1_formula, linear_predictor ~ 
                                         . + 0)
      call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
      Q2_formula_full = update.formula(mod$Q2_formula, linear_predictor ~ 
                                         . + 0)
      call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
      mod$effects = c(Return$effects, list(call_Q1 = call_Q1, 
                                           call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
    }
    if (!is.null(mod$covariate_data)) {
      covariate_data_full = data.frame(mod$covariate_data, linear_predictor = 0)
      X1_formula_full = update.formula(mod$X1_formula, linear_predictor ~ 
                                         . + 0)
      call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
      X2_formula_full = update.formula(mod$X2_formula, linear_predictor ~ 
                                         . + 0)
      call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
      mod$effects = c(mod$effects, list(call_X1 = call_X1, 
                                        call_X2 = call_X2, covariate_data_full = covariate_data_full))
    }
  }
  
  # Identify formula-specific stuff
  if( which_formula=="X1" ){
    formula_orig = mod$X1_formula
    parname = "gamma1_cp"
    mod$call = mod$effects$call_X1
  }else if( which_formula=="X2" ){
    formula_orig = mod$X2_formula
    parname = "gamma2_cp"
    mod$call = mod$effects$call_X2
  }else if( which_formula=="Q1" ){
    formula_orig = mod$Q1_formula
    parname = "lambda1_k"
    mod$call = mod$effects$call_Q1
  }else if( which_formula=="Q2" ){
    formula_orig = mod$Q2_formula
    parname = "lambda2_k"
    mod$call = mod$effects$call_Q2
  }else{
    stop("Check `which_formula` input")
  }
  
  # Extract parameters / covariance
  whichnum = which(names(mod$parameter_estimates$par)==parname)
  mod$parhat = mod$parameter_estimates$par[whichnum]
  mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum,whichnum,drop=FALSE]
  # Fill in values that are mapped off
  if( parname %in% names(mod$tmb_list$Obj$env$map) ){
    mod$parhat = mod$parhat[ mod$tmb_list$Obj$env$map[[parname]] ]
    mod$covhat = mod$covhat[ mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]], drop=FALSE ]
    mod$parhat = ifelse( is.na(mod$parhat), 0, mod$parhat)
    mod$covhat = ifelse( is.na(mod$covhat), 0, mod$covhat)
  }
  # add names
  names(mod$parhat)[] = parname
  rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)
  
  # Augment stuff
  formula_full = update.formula(formula_orig, linear_predictor~.+0)
  mod$coefficients = mod$parhat
  mod$vcov = mod$covhat
  mod$formula = formula_full
  mod$family = gaussian(link = "identity")
  
  # Functions for package
  family.fit_model = function(x,...) x$family
  vcov.fit_model = function(x,...) x$vcov
  
  # dummy functions to make Effect.default work
  dummyfuns = list(variance = function(mu) mu,
                   initialize = expression(mustart = y + 0.1),
                   dev.resids = function(...) poisson()$dev.res(...) )
  
  # Replace family (for reasons I don't really understand)
  fam = mod$family
  for( i in names(dummyfuns) ){
    if( is.null(fam[[i]]) ) fam[[i]] = dummyfuns[[i]]
  }
  
  # allow calculation of effects ...
  if (length(formals(fam$variance))>1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance = dummyfuns$variance
  }
  
  # Bundle arguments
  args = list(call = mod$call,
              coefficients = mod$coefficients,
              vcov = mod$vcov,
              family = fam,
              formula = formula_full)

  # Do call
  effects::Effect.default(focal.predictors,
                          mod,
                          ...,
                          sources = args)
}

Effect.default.test<- function (focal.predictors, mod, ..., sources) {
  if(FALSE){
    focal.predictors = focal.predictors
    mod = mod
    sources = args
  }
  sources <- if (missing(sources)) 
    effSources(mod)
  else sources
  formula <- if (is.null(sources$formula)) 
    insight::find_formula(mod)$conditional
  else sources$formula
  if (is.null(focal.predictors)) 
    return(formula)
  cl <- if (is.null(sources$call)) {
    if (isS4(mod)) 
      mod@call
    else mod$call
  }
  else sources$call
  cl$formula <- formula
  type <- if (is.null(sources$type)) 
    "glm"
  else sources$type
  fam <- try(family(mod), silent = TRUE)
  if (inherits(fam, "try-error")) 
    fam <- NULL
  if (!is.null(sources$family)) {
    fam <- sources$family
  }
  if (!is.null(fam)) {
    fam$aic <- function(...) NULL
    if (!is.null(fam$variance)) {
      if (length(formals(fam$variance)) > 1) 
        stop("Effect plots are not implemented for families with more than\n             one parameter in the variance function (e.g., negitave binomials).")
    }
  }
  cl$family <- fam
  coefficients <- if (is.null(sources$coefficients)) 
    effCoef(mod)
  else sources$coefficients
  vcov <- if (is.null(sources$vcov)) 
    as.matrix(vcov(mod, complete = TRUE))
  else sources$vcov
  zeta <- if (is.null(sources$zeta)) 
    NULL
  else sources$zeta
  cl$control <- switch(type, glm = glm.control(epsilon = Inf, 
                                               maxit = 1), polr = list(maxit = 1), multinom = c(maxit = 1))
  cl$method <- sources$method
  .m <- switch(type, glm = match(c("formula", "data", "family", 
                                   "contrasts", "subset", "control", "offset"), names(cl), 
                                 0L), polr = match(c("formula", "data", "family", "contrasts", 
                                                     "subset", "control", "method"), names(cl), 0L), multinom = match(c("formula", 
                                                                                                                        "data", "family", "contrasts", "subset", "family", "maxit", 
                                                                                                                        "offset"), names(cl), 0L))
  cl <- cl[c(1L, .m)]
  cl[[1L]] <- as.name(type)
  mod2 <- eval(cl)
  mod2$coefficients <- coefficients
  mod2$vcov <- vcov
  if (!is.null(zeta)) 
    mod2$zeta <- zeta
  if (type == "glm") {
    mod2$weights <- as.vector(with(mod2, prior.weights * 
                                     (family$mu.eta(linear.predictors)^2/family$variance(fitted.values))))
  }
  class(mod2) <- c("fakeeffmod", class(mod2))
  Effect(focal.predictors, mod2, ...)
}

Effect.test<- function (focal.predictors, mod, ...) 
{
  if(FALSE){
    focal.predictors = focal.predictors
    mod = mod
  }
  if (!checkFormula(mod)) 
    stop("model formula should not contain calls to", "\n  factor(), as.factor(), ordered(), as.ordered(),", 
         " as.numeric(), or as.integer();", "\n  see 'Warnings and Limitations' in ?Effect")
  UseMethod("Effect", mod)
}
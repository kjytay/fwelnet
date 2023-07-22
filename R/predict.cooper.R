#' Survival prediction for COoper models.
#'
#' @param object Object of class `cooper`, containing two cause-specific Cox model fits.
#' @param xnew `data.frame` with new observations for predictions.
#' @param type Either `"risk"` for linear predictors, or `absrisk` (default) for CIF.
#' @param event Event of interest, either `1` (default) or `2`.
#' @param times Event times to evaluate at.
#' @param use_initial_fit (`FALSE`) If `TRUE`, uses the initial `glmnet` fit of the cooper algorithm.
#' @param reference **NYI** Either `"zero"` or `"sample"`.
#'
#' @return For `type = "risk"`, a data.frame with columns "id", "event", "lp", "eXb".
#' For `type = "absrisk"`, a data.frame with columns "id" "event" "horizon" "absolute_risk".
#' @export
#'
#' @examples
#' 
#' data(pbc, package = "survival")
#' pbc <- na.omit(pbc)
#' xtrain <- pbc[1:270, -1]
#' xtest <- pbc[270:276, -1]
#' 
#' set.seed(12)
#' cooperfit <- fwelnet_mt_cox(
#'   xtrain, mt_max_iter = 3,
#'   stratify_by_status = TRUE,
#'   standardize = FALSE # Currently required to ensure correct results
#' )
#' 
#' predict(cooperfit, xnew = xtest, type = "risk", event = 1)
#' predict(cooperfit, xnew = xtest, times = c(1180.5, 1788), type = "absrisk", event = 1)
predict.cooper <- function(
    object, xnew, 
    type = c("absrisk", "risk"), 
    event = 1,
    times = NULL,
    use_initial_fit = FALSE,
    reference = c("zero", "sample") # NYI
) {
  
  type <- match.arg(type)
  reference <- match.arg(reference)
  #browser()
  # Accessing either initial glmnet fit or last fwelnet fits
  # Names of objects in cooper return value list
  which_model <- ifelse(use_initial_fit, "initial_fit", "fwelfits")
  
  # subset to ensure only predictors are contained in new data
  # if new data is a data.table this would need with = FALSE
  if (inherits(xnew, "data.table")) xnew <- as.data.frame(xnew)
  xnew[, object$predictors, drop = FALSE]

  if (type == "risk") {
    # If event is length 2, lapply
    linpred <- do.call(rbind, lapply(event, function(e) {
      lp <- predict(object[[which_model]][[e]], xnew = x_modelmat(object, xnew), type = "link", s = "lambda.min")
      lp <- as.vector(lp)
      
      data.frame(
        id = seq_along(lp),
        event = e,
        lp = lp,
        eXb = exp(lp)
      )
    }))
   
    return(linpred)
  }
  
  if (is.null(times)) stop("`times` needs to be supplied for absolute risk predictions.")
  
  time_grid <- expand.grid(event = event, horizon = times)

  # I miss purrr.
  do.call(rbind, Map(
    function(event, horizon) {
      get_abs_risk(object, xnew, event, horizon, use_initial_fit)
    },
    event = time_grid$event,
    horizon = time_grid$horizon
  ))
  
}

# Simple utility to convert a data.frame input of new data to model matrix,
# which takes care of dummy encoding etc. Note the intercept is dropped explicitly rather than removed
# via formula to ensure correct recoding.
# Also ensures only predictors used during model fit are in model matrix
x_modelmat <- function(object, xnew) {
  model.matrix(reformulate(object$predictors), data = xnew)[, -1]
}

get_abs_risk <- function(object, xnew, event, horizon, use_initial_fit = FALSE) {
  
  event_times <- get_event_times(object$y[[event]], horizon = horizon)
  #browser()
  base_haz <- baseline_hazard(object, event = event, use_initial_fit = use_initial_fit)
  base_haz <- base_haz[base_haz$time %in% event_times, "hazard", drop = FALSE]
  # need is as an 1 x times matrix
  base_haz <- t(base_haz)
  
  # Get cumhazards H_e for both event, subset to event times, to get marginal hazard, survival
  cumhazards <- lapply(1:2, function(e) {
    cumulative_hazard(object, xnew, event = e, times = event_times, use_initial_fit = use_initial_fit)
  })
  
  # Get exp(eta), vector of length nrow(xnew)
  risks <- predict(object, event = event, xnew = xnew, type = "risk")[["eXb"]]

  marginal_chazard <- Reduce(`+`, cumhazards)


  marginal_surv <- exp(-marginal_chazard)
  # (exp(eta) * h0(t)) * S(t)  -> he(t) * S(t)
  #abs_risk <- as.vector(risks * as.vector(base_haz) * marginal_surv)
  abs_risk <- as.vector(risks * base_haz %*% marginal_surv)
  
  data.frame(
    id = seq_along(abs_risk),
    event = event,
    horizon = horizon,
    absolute_risk = abs_risk
  )
}

# Get event times (i.e. status != 0) from Surv() object
# Since these models have status in 0/1, the event param is kind of superfluous and only 1 is used.
get_event_times <- function(surv, horizon = NULL) {
  if (is.null(horizon)) horizon <- max(surv[, 1])
  
  idx <- which((surv[, 1] <= horizon) & (surv[, 2] == 1))
  surv[idx, 1]
}

#' predictRisk for Cooper models
#'
#' @inheritParams predict.cooper
#' @param newdata New data to predict on.
#' @param times Event times to evaluate at
#' @param cause Event of interest (`1` or `2`).
#' @param ... Unused
#'
#' @return Matrix of dimension `nrow(newdata)` x `length(times)`
# @export
#' @exportS3Method riskRegression::predictRisk
#' @examples
#' if (requireNamespace("riskRegression")) {
#' library(riskRegression)
#' 
#' data(pbc, package = "survival")
#' pbc <- na.omit(pbc)
#' xtrain <- pbc[1:270, -1]
#' xtest <- pbc[270:276, -1]
#' 
#' set.seed(12)
#' cooperfit <- fwelnet_mt_cox(
#'   xtrain, mt_max_iter = 3,
#'   stratify_by_status = TRUE,
#'   standardize = FALSE # Currently required to ensure correct results
#' )
#' 
#' predictRisk(cooperfit, xtest, cause = 1, times = c(1180.5, 1788))
#' 
#' }
predictRisk.cooper <- function(object, newdata, times, cause, use_initial_fit = FALSE, ...) {
  checkmate::assert_int(cause, lower = 1, upper = 2)
  
  res <- do.call(cbind, lapply(times, function(time) {
    risk <- predict(
      object, 
      xnew = newdata, 
      type = "absrisk", 
      event = cause, 
      times = time, 
      use_initial_fit = use_initial_fit
    )[["absolute_risk"]]
    as.matrix(risk)
  }))
  
  colnames(res) <- times
  res
}


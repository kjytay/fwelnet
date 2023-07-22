
baseline_hazard <- function(object, event = 1, use_initial_fit = FALSE, times = NULL) {
  
  # get baseline hazards using x = 0 as reference, possibly bad idea according to ?basehaz
  newdata_zero <- object$x[1, ] * 0
  
  # Accesses the glmnet fit returned with the fwelnet fit with the cv.fwelnet object (naming is bad here yes)
  
  if (use_initial_fit) {
    model_fit <- object$initial_fits[[event]]
    lambda_min <- model_fit$lambda.min
  } else {
    model_fit <- object$fwelfits[[event]]$glmfit$glmfit
    lambda_min <- object$fwelfits[[event]]$lambda.min
  }
  
  
  bhraw <- survival::survfit(
    model_fit, s = lambda_min, 
    x = object$x,
    y = object$y[[event]],
    newx = newdata_zero,
    se.fit = FALSE
  )
  
  ret <- data.frame(
    event = event,
    time = bhraw$time,
    n.risk = bhraw$n.risk,
    surv = bhraw$surv,
    cumhaz = bhraw$cumhaz,
    hazard = diff(c(0, bhraw$cumhaz))
  )
  #browser()
  if (!is.null(times)) {
    ret <- ret[which(ret$time %in% times), ]
  }
  
  ret
}


cumulative_hazard <- function(object, xnew, event = 1, times = NULL, use_initial_fit = FALSE) {
  
  if (use_initial_fit) {
    # cv.glmnet object with lambda.min
    model_fit <- object$initial_fits[[event]]
    lambda_min <- model_fit$lambda.min
  } else {
    # glmnet fit object within fwelnet object within cv.fwelnet object
    model_fit <- object$fwelfits[[event]]$glmfit$glmfit
    lambda_min <- object$fwelfits[[event]]$lambda.min
  }
  # browser()
  bfit <- survival::survfit(
    model_fit, 
    x = object$x,
    y = object$y[[event]],
    newx = x_modelmat(object, xnew), s = lambda_min
  )

  ret <- bfit$cumhaz
  
  if (!is.null(times)) {
    # Make sure to subset by [rows, columns, drop = FALSE]
    # Don't forget the second , or you will have a bad time
    # Ask me how I know
    ret <- ret[which(bfit$time %in% times), , drop = FALSE]
  }
  
  ret
}

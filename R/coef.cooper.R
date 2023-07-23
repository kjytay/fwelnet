#' @export
#' @importFrom stats setNames
coef.cooper <- function(object, event = 1, s = "lambda.min", use_initial_fit = FALSE,...) {
  

  if (!use_initial_fit) {
    model_which <- "fwelfits"
    
    if (s == "lambda.min") {
      s <- which(object[[model_which]][[event]]$glmfit$lambda == object[[model_which]][[event]]$lambda.min)
    }
    
    coefs <- object[[model_which]][[event]]$glmfit$beta[, s]
  } else {
    model_which <- "initial_fits"
    
    if (s == "lambda.min") {
      s <- which(object[[model_which]][[event]]$lambda == object[[model_which]][[event]]$lambda.min)
    }
    
    coefs <- object[[model_which]][[event]]$glmnet.fit$beta[, s]
  }
  
  stats::setNames(
    object = coefs,
    nm = colnames(object$x)
  )
  
}


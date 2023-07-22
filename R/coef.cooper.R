#' @export
#' @importFrom stats setNames
coef.cooper <- function(object, event = 1, s = "lambda.min", use_initial_fit = FALSE,...) {
  

  # We may want to either get predictions/coefs from the initital glmnet fit (for comparison)
  # or the final fit. We're not storing intermediate models/lambdas, only first and last.
  # Hence, "iteration" 
  
    
  if (!use_initial_fit) {
    
    model_which <- "fwelfits"
    
    if (s == "lambda.min") {
      s <- which(object[[model_which]]$glmfit$lambda == object[[model_which]]$lambda.min)
    }
    
    coefs <- object[[model_which]]$glmfit$beta[, s]
  } else {
    model_which <- "initial_fit"
    
    if (s == "lambda.min") {
      s <- which(object[[model_which]]$lambda == object[[model_which]]$lambda.min)
    }
    
    coefs <- object[[model_which]]$glmnet.fit$beta[, s]
  }
  
  stats::setNames(
    object = coefs,
    nm = colnames(object$x)
  )
  
}



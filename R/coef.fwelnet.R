#' @importFrom Matrix Matrix
#' @export
coef.fwelnet <- function(object, ...) {

    if (object$family == "cox") {
      out <- Matrix::Matrix(object$beta, sparse = TRUE)
      row.names(out) <- paste0("V", 1:nrow(object$beta))
    } else {
      out <- Matrix::Matrix(rbind(object$a0, object$beta), sparse = TRUE)
      row.names(out) <- c("(Intercept)", paste0("V", 1:nrow(object$beta)))
    }

    return(out)
}

#' @export
coef.cooper <- function(object, event = 1, s = "lambda.min", ...) {
  
  fwfit_which <- paste0("fwfit", event)
  
  if (s == "lambda.min") {
    s <- which(object[[fwfit_which]]$glmfit$lambda == object[[fwfit_which]]$lambda.min)
  }
  
  setNames(
    object = object[[fwfit_which]]$glmfit$beta[, s],
    nm = colnames(cooperfit$x)
  )
}



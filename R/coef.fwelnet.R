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

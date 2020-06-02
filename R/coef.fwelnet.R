#' @importFrom Matrix Matrix
#' @export
coef.fwelnet <- function(object, ...) {
    out <- Matrix::Matrix(rbind(object$a0, object$beta), sparse = TRUE)
    row.names(out) <- c("(Intercept)", paste0("V", 1:nrow(object$beta)))
    return(out)
}

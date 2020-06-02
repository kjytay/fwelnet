#' Make predictions from a "cv.fwelnet" object
#'
#' This function returns the predictions for a new data matrix from a
#' cross-validated \code{fwelnet} model by using the stored "\code{glmfit}"
#' object and the optimal value chosen for \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object Fitted "\code{cv.fwelnet}" object.
#' @param xnew Matrix of new values for \code{x} at which predictions are to
#' be made.
#' @param s Value of the penalty parameter \code{lambda} at which predictions are
#' required. Default is the value \code{s="lambda.1se"} stored in the CV
#' \code{fit}. Alternatively, \code{s="lambda.min"} can be used.
#' @param ... Other arguments to be passed to \code{predict.fwelnet}.
#'
#' @return Predictions which the cross-validated model makes for \code{xnew} at
#' the optimal value of \code{lambda}. Note that the default is the "lambda.1se" for lambda,
#' to make this function consistent with \code{cv.glmnet} in the \code{glmnet}
#' package.
#'
#' @seealso \code{\link{cv.fwelnet}} and \code{\link{predict.fwelnet}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' z <- cbind(1, abs(beta) + rnorm(p))
#'
#' cvfit <- cv.fwelnet(x, y, z)
#' predict(cvfit, xnew = x[1:5, ])
#' predict(cvfit, xnew = x[1:5, ], s = "lambda.min")
#'
#' @importFrom stats predict
#' @export
predict.cv.fwelnet <- function(object, xnew, s = c("lambda.1se", "lambda.min"),
                               ...) {
    s <- match.arg(s)
    lambda <- object[[s]]
    idx <- which(object$lambda == lambda)

    predictions <- predict(object$glmfit, xnew, ...)
    return(predictions[, idx, drop = FALSE])
}

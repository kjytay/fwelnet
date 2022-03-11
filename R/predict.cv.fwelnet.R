#' Make predictions from a "cv.fwelnet" object
#'
#' This function returns the predictions for a new data matrix from a
#' cross-validated `fwelnet` model by using the stored "`glmfit`"
#' object and the optimal value chosen for `lambda`.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object Fitted "`cv.fwelnet`" object.
#' @param xnew Matrix of new values for `x` at which predictions are to
#' be made.
#' @param s Value of the penalty parameter `lambda` at which predictions are
#' required. Default is the value `s="lambda.1se"` stored in the CV
#' `fit`. Alternatively, `s="lambda.min"` can be used.
#' @param ... Other arguments to be passed to `predict.fwelnet`.
#'
#' @return Predictions which the cross-validated model makes for `xnew` at
#' the optimal value of `lambda`. Note that the default is the "lambda.1se" for lambda,
#' to make this function consistent with `cv.glmnet` in the `glmnet`
#' package.
#'
#' @seealso [cv.fwelnet()] and [predict.fwelnet()].
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

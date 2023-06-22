#' Make predictions from a "fwelnet" object
#'
#' This function returns the predictions from a "`fwelnet`" object
#' for a new data matrix.
#'
#' @param object Fitted "`fwelnet`" object.
#' @param xnew Matrix of new values for `x` at which predictions are to
#' be made.
#' @param type Type of prediction required. Type `"link"` (default) gives the
#' linear predictors. Type `"response"` gives the linear predictor for "gaussian"
#' and fitted probabilities for "binomial" family.
#' `"coefficients"` will return the beta vector or matrix.
#' @param ... Potentially other arguments to be passed to and from methods;
#' currently not in use.
#'
#' @return Predictions which the model `object` makes at `xnew`.
#'
#' @seealso [fwelnet()].
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' z <- cbind(1, abs(beta) + rnorm(p))
#'
#' fit <- fwelnet(x, y, z)
#' predict(fit, xnew = x[1:5, ])
#'
#' biny <- ifelse(y > 0, 1, 0)
#' fit <- fwelnet(x, biny, z, family = "binomial")
#' # linear predictor values
#' predict(fit, xnew = x[1:5, ])
#' # fitted probabilities
#' predict(fit, xnew = x[1:5, ], type = "response")
#'
#' @export
#' @importFrom stats model.matrix coef
predict.fwelnet <- function(object, xnew, type = c("link", "response", "coefficients"),
                            ...) {
    type <- match.arg(type)

    # Allow data.frame/categorical input for survival without intercept
    if (is.data.frame(xnew)) {
      if (object$family == "cox") {
        #xnew <- model.matrix(~ . -1, xnew)
        xnew <- model.matrix(~ ., xnew)[, -1, drop = FALSE]
      } else {
        xnew <- model.matrix(~ ., xnew)
      }
    }

    # Coefs only as glmnet does, no lambda interpolation for now
    if (type == "coefficients") {
      out <- object$beta
      return(out)
    }

    if (object$family == "cox") {
        # a0 is numeric(0) or NULL for family = "cox"
        out <- xnew %*% object$beta
    } else {
        out <- t(object$a0 + t(xnew %*% object$beta))
    }

    if (type == "response" && object$family == "binomial") {
        out <- 1 / (1 + exp(-out))
    }

    return(out)
}

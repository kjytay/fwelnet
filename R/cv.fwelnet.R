#' Cross-validation for fwelnet
#'
#' Does `k`-fold cross-validation for `fwelnet`.
#'
#' This function runs `fwelnet nfolds+1` times: the first to get the
#' `lambda` sequence, and the remaining `nfolds` times to compute the
#' fit with each of the folds omitted. The error is accumulated, and the mean
#' error and standard deviation over the folds is computed. Note that
#' `cv.pcLasso` does NOT search for values of `alpha`. A specific
#' value of `alpha` should be supplied.
#'
#' @param x `x` matrix as in `fwelnet`.
#' @param y `y` matrix as in `fwelnet`.
#' @param z `z` matrix as in `fwelnet`.
#' @param family Response type. Either `"gaussian"` (default) for linear
#' regression or `"binomial"` for logistic regression.
#' @param lambda A user supplied `lambda` sequence. Typical usage is to
#' have the program compute its own `lambda` sequence; supplying a value of
#' lambda overrides this.
#' @param type.measure Loss to use for cross-validation. Currently five options,
#' not all available for all models. The default is type.measure="deviance",
#' which uses squared-error for gaussian models (a.k.a type.measure="mse" there)
#' and deviance for logistic and cox regression.
#' `type.measure = "class"` applies to binomial logistic regression only,
#' and gives misclassification error. `type.measure = "auc"` is for two-class
#' logistic regression only, and gives area under the ROC curve.
#' `type.measure = "mse"` or type.measure="mae" (mean absolute error) can be used by
#' all models except cox.
#' @param nfolds Number of folds for CV (default is 10). Although `nfolds`
#' can be as large as the sample size (leave-one-out CV), it is not recommended
#' for large datasets. Smallest value allowable is `nfolds = 3`.
#' @param foldid An optional vector of values between 1 and `nfolds`
#' identifying what fold each observation is in. If supplied, `nfolds` can
#' be missing.
#' @param keep If `keep = TRUE`, a prevalidated array is returned
#' containing fitted values for each observation at each value of lambda. This
#' means these fits are computed with this observation and the rest of its fold
#' omitted. Default is `FALSE`.
#' @param verbose Print information as model is being fit? Default is FALSE.
#' @param ... Other arguments that can be passed to `fwelnet`.
#' @inheritParams fwelnet
#'
#' @return An object of class `"cv.fwelnet"`, which is a list with the
#' ingredients of the cross-validation fit.
#' \item{glmfit}{A fitted `fwelnet` object for the full data.}
#' \item{lambda}{The values of `lambda` used in the fits.}
#' \item{nzero}{The number of non-zero coefficients in the model `glmfit`.}
#' \item{fit.preval}{If `keep=TRUE`, this is the array of prevalidated
#'   fits.}
#' \item{cvm}{The mean cross-validated error: a vector of length
#'   `length(lambda)`.}
#' \item{cvsd}{Estimate of standard error of `cvm`.}
#' \item{cvlo}{Lower curve = `cvm - cvsd`.}
#' \item{cvup}{Upper curve = `cvm + cvsd`.}
#' \item{lambda.min}{The value of `lambda` that gives minimum
#'   `cvm`.}
#' \item{lambda.1se}{The largest value of `lambda` such that the CV
#'   error is within one standard error of the minimum.}
#' \item{foldid}{If `keep=TRUE`, the fold assignments used.}
#' \item{name}{Name of error measurement used for CV.}
#' \item{call}{The call that produced this object.}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' z <- cbind(1, abs(beta) + rnorm(p))
#'
#' cvfit1 <- cv.fwelnet(x, y, z)
#'
#' # change no. of CV folds
#' cvfit2 <- cv.fwelnet(x, y, z, nfolds = 5)
#' # specify which observations are in each fold
#' foldid <- sample(rep(seq(5), length = length(y)))
#' cvfit3 <- cv.fwelnet(x, y, z, foldid = foldid)
#' # keep=TRUE to have pre-validated fits and foldid returned
#' cvfit4 <- cv.fwelnet(x, y, z, keep = TRUE)
#'
#' @importFrom stats predict
#' @importFrom glmnet glmnet coxnet.deviance
#' @export
cv.fwelnet <- function(x, y, z, family = c("gaussian", "binomial", "cox"), lambda = NULL,
                       type.measure = c("mse", "deviance", "class", "auc", "mae", "nll"),
                       nfolds = 10, foldid = NULL, keep = FALSE, verbose = FALSE, t = 1, a = 0.5,
                       ...) {
    this.call <- match.call()

    n <- nrow(x); p <- ncol(x); K <- ncol(z)
    family <- match.arg(family)
    if (family != "cox") {
      # Previous default behavior unsuitable for y if Surv() object
      y <- as.vector(y)
    }

    # get type.measure
    if (missing(type.measure)) {
      # Previous behaviour
      # type.measure <- ifelse(family == "gaussian", "mse", "deviance")

      # Setting cox -> deviance, and deviance for binomial as previous default
        type.measure <- switch (family,
          "gaussian" = "mse",
          "binomial" = "deviance",
          "cox" = "deviance"
        )
    } else {
        type.measure <- match.arg(type.measure)
    }

    # if fold IDs not given, randomly create them
    # if given, count the number of folds
    if (is.null(foldid)) {
        foldid <- sample(rep(seq(nfolds), length = n))
    } else {
        nfolds <- max(foldid)
    }
    if (nfolds < 3) {
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }

    # get fwelnet fit for all of the data
    fit0 <- fwelnet(x, y, z, lambda = lambda, family = family, t = t, a = a, ...)
    if (verbose) {
        cat("Initial fit done", fill = TRUE)
    }
    nz <- colSums(abs(fit0$beta) > 0)

    # fit fwelnet on folds
    fits <- vector("list", nfolds)
    for (ii in 1:nfolds) {
        if (verbose) {
            cat(c("Fold = ", ii), fill = TRUE)
        }
        oo <- foldid == ii
        xc <- x[!oo, , drop = FALSE]
        yy <- y[!oo]
        fits[[ii]] <- fwelnet(xc, yy, z, lambda = fit0$lambda, family = family,
                               verbose = verbose, t = t, a = a, ...)
    }

    # get predictions
    predmat <- matrix(NA, n, length(fit0$lambda))
    for (ii in 1:nfolds) {
      oo <- foldid == ii

      if (family != "cox") {
        out <- predict(fits[[ii]], x[oo, , drop = F])
        predmat[oo, 1:ncol(out)] <- out
      }

      if (family == "cox") {
        # If cox, add nll here?
        # fit on test data
        fits[[ii]]$nll <- glmnet::coxnet.deviance(
          x = x[oo, , drop = F],
          y = y[oo, , drop = F],
          beta = fits[[ii]]$beta
        )

        # Get number of events in fold, if yy is Surv() obj. second column
        # is status with 1 == event
        # fits[[ii]]$n_events <- sum(y[oo, "status", drop = F])
      }
    }

    # compute CV stuff (machinery borrowed from glmnet)
    weights <- rep(1, n); grouped = TRUE  # for compatibility w glmnet
    if (family == "gaussian") {
        cvstuff <- cv.elnet(predmat, y, type.measure, weights, foldid, grouped)
        name <- ifelse(type.measure == "mae", "Mean Absolute Error",
                       "Mean-Squared Error")
    } else if (family == "binomial") {
        cvstuff <- cv.lognet(predmat, y, type.measure, weights, foldid, grouped)
        name <- switch(type.measure,
                       deviance = "Binomial Deviance",
                       class = "Misclassification Error",
                       auc = "AUC")
    } else if (family == "cox") {
      # nll calc done above during fit step so access x, y in folds
      # or better to use glmnet::coxnet.deviance probably
      # nll per fold is fine for now unless it isn't?
      cvstuff <- list(
        # usually N x nlam matrix, but we only have nfold x nlam, so yeah.
        # This is not a nice way to do it but it works, I guess
        cvraw = t(vapply(fits, function(x) x$nll, FUN.VALUE = double(length(fit0$lambda)))),
        weights = weights, # See above
        # Number of events in fold maybe rather than raw N?
        # N = vapply(fits, function(x) x$n_events, FUN.VALUE = double(1)),
        type.measure = type.measure,
        grouped = TRUE # for compatibility? see above
      )
      name <- "Cox Deviance"
    }

    # Hacky special handling for the somewhat incomplete cox case
    # FIXME: All of this should just use glmnet machinery but time is finite
    if (family == "cox") {
      # average deviance per lambda across folds (rows = folds, cols = lambdas)
      mean_nll_per_lambda <- apply(cvstuff$cvraw, 2, mean)
      # index of lambda with smallest nll
      best_lambda_i <- which.min(mean_nll_per_lambda)
      lambda.min <- fit0$lambda[best_lambda_i]

      out <- list(
        glmfit = fit0,
        lambda = fit0$lambda,
        nzero = fit0$nzero,
        # fit.preval=predmat.preval,
        cvm = colMeans(cvstuff$cvraw),
        # cvsd=cvsd,
        # cvlo=cvlo,
        # cvup=cvup,
        lambda.min = lambda.min,
        # lambda.1se=lambda.1se,
        # foldid=foldid_copy,
        name = name,
        call = this.call
      )
      class(out) <- "cv.fwelnet"

      return(out)
    }

    cvout <- cvstats(cvstuff, foldid, nfolds, lambda, nz, cvstuff$grouped)
    cvm <- cvout$cvm
    cvsd <- cvout$cvsd
    cvup <- cvout$cvup
    cvlo <- cvout$cvlo

    predmat.preval <- NULL
    foldid_copy <- NULL
    if (keep) {
        predmat.preval <- predmat
        foldid_copy <- foldid
    }

    # get best and 1se lambda
    if (type.measure == "auc") {
        imin <- which.max(cvm)
        lambda.min <- fit0$lambda[imin]
        imin.1se <- which(cvm > cvm[imin] - cvsd[imin])[1]
        lambda.1se <- fit0$lambda[imin.1se]
    } else {
        imin <- which.min(cvm)
        lambda.min <- fit0$lambda[imin]
        imin.1se <- which(cvm < cvm[imin] + cvsd[imin])[1]
        lambda.1se <- fit0$lambda[imin.1se]
    }

    out <- list(glmfit=fit0, lambda=fit0$lambda, nzero=fit0$nzero,
                fit.preval=predmat.preval, cvm=cvm, cvsd=cvsd, cvlo=cvlo,
                cvup=cvup, lambda.min=lambda.min, lambda.1se=lambda.1se,
                foldid=foldid_copy, name=name, call=this.call)
    class(out) <- "cv.fwelnet"
    return(out)
}

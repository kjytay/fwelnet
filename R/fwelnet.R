#' Fit feature-weighted elastic net
#'
#' Fit a model with feature-weighted elastic net for a path of lambda values.
#' Fits linear and logistic regression models.
#'
#' \code{fwelnet} always mean centers the columns of the \code{x} matrix. If
#' \code{standardize=TRUE}, \code{fwelnet} will also scale the columns to have
#' standard deviation 1. In all cases, the \code{beta} coefficients returned are
#' for the original \code{x} values (i.e. uncentered and unscaled).
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is
#' an observation vector.
#' @param y Response variable. Quantitative for \code{family = "gaussian"}. For
#' \code{family="binomial"}, should be a numeric vector consisting of 0s and 1s.
#' @param z Feature of features matrix, with dimension \code{nvars x nfeaturevars}.
#' @param lambda A user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence; supplying a value of
#' lambda overrides this.
#' @param family Response type. Either \code{"gaussian"} (default) for linear
#' regression or \code{"binomial"} for logistic regression.
#' @param alpha The elastic net mixing hyperparameter, a real value number
#' between 0 and 1 (inclusive). Default value is 1.
#' @param standardize If \code{TRUE}, the columns of the input matrix are
#' standardized before the algorithm is run. Default is \code{TRUE}.
#' @param max_iter The number of iterations for the optimization. Default is 1.
#' @param ave_mode If equal to 1 (default), the gradient descent direction for
#' \code{theta} is the mean gradient across the lambda values. If equal to 2,
#' it is the component-wise median gradient across the lambda values.
#' @param thresh_mode If equal to 1 (default), backtracking line search for
#' \code{theta} is done so that the mean objective function (across lambda
#' values) decreases. If equal to 2, it is done so that the median objective
#' function decreases.
#' @param t The initial step size for \code{theta} backtracking line search
#' (default value is 1).
#' @param a The factor by which step size is decreased in \code{theta}
#' backtracking line search (default value is 0.5).
#' @param thresh If the mean/median objective function does not decrease by at
#' least this factor, we terminate the optimization early. Default is 1e-4.
#' @param verbose If \code{TRUE}, prints information to console as model is
#' being fit. Default is \code{FALSE}.
#'
#' @return An object of class \code{"fwelnet"}.
#' \item{beta}{A \code{p x length(lambda)} matrix of coefficients.}
#' \item{theta}{Theta value, a \code{nfeaturevars x 1} matrix.}
#' \item{a0}{Intercept sequence of length \code{length(lambda)}.}
#' \item{lambda}{The actual sequence of \code{lambda} values used.}
#' \item{nzero}{The number of non-zero coefficients for each value of
#' \code{lambda}.}
#' \item{family}{Response type.}
#' \item{call}{The call that produced this object.}
#' \item{obj}{A \code{max_iter + 1} by \code{length(lambda)} matrix of
#' objective function values. (Number of rows could be fewer if the optimization
#' stopped early.)}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' z <- cbind(1, abs(beta) + rnorm(p))
#'
#' fwelnet(x, y, z)
#' fwelnet(x, y, z, ave_mode = 2)
#' fwelnet(x, y, z, ave_mode = 2, thresh_mode = 2)
#'
#' @importFrom stats median sd
#' @export
fwelnet <- function(x, y, z, lambda = NULL, family = c("gaussian", "binomial", "cox"),
                    alpha = 1, standardize = TRUE, max_iter = 1, ave_mode = 1,
                    thresh_mode = 1, t = 1, a = 0.5, thresh = 1e-4,
                    verbose = FALSE) {
    this.call <- match.call()

    if (alpha > 1 || alpha < 0) {
        stop("alpha must be between 0 and 1 (inclusive)")
    }

    n <- nrow(x); p <- ncol(x); K <- ncol(z)
    if (!inherits(y, "Surv")) {
      # Need to pass survival outcome via Surv or as matrix, this breaks it
      y <- as.vector(y)
    }
    family <- match.arg(family)
    if (family == "binomial" && any(!(unique(y) %in% c(0, 1)))) {
        stop("If family is binomial, y can only contain 0s and 1s")
    }

    # center y and columns of x
    mx <- colMeans(x)
    if (standardize) {
        sx <- apply(x, 2, sd) * sqrt((n-1) / n)
    } else {
        sx <- rep(1, ncol(x))
    }
    x <- scale(x, mx, sx)

    y_mean <- NA
    if (family == "gaussian") {
        y_mean <- mean(y)
        y <- y - y_mean
    }

    # fit elastic net; if lambda sequence not provided, use the glmnet default
    # come up with lambda values: use the lambda sequence from glmnet
    if (is.null(lambda)) {
        glmfit <- glmnet::glmnet(x, y, alpha = alpha, family = family,
                                 standardize = FALSE)
        lambda <- glmfit$lambda
    } else {
        glmfit <- glmnet::glmnet(x, y, alpha = alpha, family = family,
                                 standardize = FALSE, lambda = lambda)
    }

    # function for averaging objective
    # get_ave_obj <- ifelse(thresh_mode == 1, mean, median)
    # ifelse throws an error in interactive testing otherwise
    get_ave_obj <- switch(thresh_mode,
        `1` = mean,
        `2` = median
    )

    ave_fn_name <- ifelse(thresh_mode == 1, "Mean", "Median")

    # initialize beta, a0 at elastic net solution, theta at 0
    beta <- matrix(glmfit$beta, nrow = p)
    a0 <- glmfit$a0
    theta <- matrix(0, nrow = K, ncol = 1)
    iter <- 0    # no. of beta/theta minimizations

    # store obj values across iterations (including original elastic net soln)
    obj_store <- matrix(NA, nrow = max_iter + 1, ncol = length(lambda))
    # FIXME:
    # 1) a0 is NULL in family = "cox" case
    # 2) objective_fn shoddily cox'd, w/o a0 though
    # browser()
    obj_store[1, ] <- objective_fn(x, y, z, beta, a0, theta, lambda, alpha, family)
    prev_iter_m_obj_value <- get_ave_obj(obj_store[1, ])
    if (verbose) cat(paste(ave_fn_name,
        "objective function value of elastic net solution:",
        format(round(prev_iter_m_obj_value, 3), nsmall = 3)),
        fill = TRUE)
    while (iter < max_iter) {
        iter <- iter + 1
        if (verbose) cat("Optimization iteration", iter, fill = TRUE)

        # OPTIMIZATION FOR THETA: ONE GRADIENT BACKTRACKING STEP
        # compute direction of gradient descent
        grad_theta <- gradient_fn(x, y, z, beta, theta, lambda, alpha)
        if (ave_mode == 1) {
            ave_grad_theta <- matrix(rowMeans(grad_theta), ncol = 1)
        } else {
            ave_grad_theta <- matrix(apply(grad_theta, 1, median), ncol = 1)
        }

        # backtracking for theta
        # compute current objective
        nll_value <- nll_fn(x, y, beta, a0, family)
        obj_value <- nll_value + penalty_fn(z, beta, theta, lambda, alpha)
        m_obj_value <- get_ave_obj(obj_value)

        while (TRUE) {
            newtheta <- theta - t * ave_grad_theta

            # compute new objective
            new_obj_value <- nll_value + penalty_fn(z, beta, newtheta, lambda, alpha)
            new_m_obj_value <- get_ave_obj(new_obj_value)

            if (!is.finite(new_m_obj_value) || new_m_obj_value > m_obj_value) {
                t <- a * t
            } else {
                theta <- newtheta
                break
            }
        }
        if (verbose) cat(paste(ave_fn_name,
            "objective function value after theta minimization:",
            format(round(new_m_obj_value, 3), nsmall = 3)),
            fill = TRUE)

        # OPTIMIZATION FOR BETA: ONE GLMNET STEP
        pen.weights = mean(exp(z %*% theta)) * exp(- z %*% theta)
        fit <- glmnet::glmnet(x, y, family = family, standardize = FALSE,
                              alpha = alpha, lambda = lambda * mean(pen.weights),
                              penalty.factor = pen.weights)
        beta <- matrix(fit$beta, nrow = p)
        a0 <- fit$a0

        # compute new objective and store it
        new_obj_value <- objective_fn(x, y, z, beta, a0, theta, lambda, alpha, family)
        new_m_obj_value <- get_ave_obj(new_obj_value)
        obj_store[iter + 1, ] <- new_obj_value
        if (verbose) cat(paste(ave_fn_name,
            "objective function value after beta minimization:",
            format(round(new_m_obj_value, 3), nsmall = 3)),
            fill = TRUE)

        # if mean/median obj function value not decreased enough, terminate early
        if ((prev_iter_m_obj_value - new_m_obj_value) /
            prev_iter_m_obj_value < thresh) {
            obj_store <- obj_store[1:(iter + 1), ]
            if (verbose) cat(paste(ave_fn_name,
                "objective function value not decreased enough, stopping early"))
            break
        } else {
            prev_iter_m_obj_value <- new_m_obj_value
        }
    }

    # rescale intercept & beta values (beta will only change when standardize = TRUE)
    a0 <- a0 - colSums(beta * mx / sx)
    if (family == "gaussian") a0 <- a0 + y_mean
    beta <- beta / sx

    # count no. of non-zero coefficients for each model
    nzero <- colSums(beta != 0)

    out <- list(beta = beta, theta = theta, a0 = a0, lambda = lambda, nzero = nzero,
                family = family, call = this.call, obj = obj_store)
    class(out) <- "fwelnet"
    return(out)
}

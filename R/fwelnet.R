#' Fit feature-weighted elastic net
#'
#' Fit a model with feature-weighted elastic net for a path of lambda values.
#' Fits linear and logistic regression models.
#'
#' `fwelnet` always mean centers the columns of the `x` matrix. If
#' `standardize=TRUE`, `fwelnet` will also scale the columns to have
#' standard deviation 1. In all cases, the `beta` coefficients returned are
#' for the original `x` values (i.e. uncentered and unscaled).
#'
#' @param x Input matrix, of dimension `nobs x nvars`; each row is
#' an observation vector.
#' @param y Response variable. Quantitative for `family = "gaussian"`. For
#' `family="binomial"`, should be a numeric vector consisting of 0s and 1s.
#' For `family = "cox"`, should be a [`survival::Surv`] object.
#' @param z Feature of features matrix, with dimension `nvars x nfeaturevars`.
#' @param lambda A user supplied `lambda` sequence. Typical usage is to
#' have the program compute its own `lambda` sequence; supplying a value of
#' lambda overrides this.
#' @param family Response type. Either `"gaussian"` (default) for linear
#' regression, `"binomial"` for logistic regression, or `"cox"` for cox
#' regression.
#' @param alpha The elastic net mixing hyperparameter, a real value number
#' between 0 and 1 (inclusive). Default value is 1.
#' @param standardize If `TRUE`, the columns of the input matrix are
#' standardized before the algorithm is run. Default is `TRUE`.
#' @param max_iter The number of iterations for the optimization.
#' @param ave_mode If equal to 1 (default), the gradient descent direction for
#' `theta` is the mean gradient across the lambda values. If equal to 2,
#' it is the component-wise median gradient across the lambda values.
#' @param thresh_mode If equal to 1 (default), backtracking line search for
#' `theta` is done so that the mean objective function (across lambda
#' values) decreases. If equal to 2, it is done so that the median objective
#' function decreases.
#' @param weight_fun A function relating `z` and `theta` resulting in penalization
#'  weights. Must be either `"original"` for the function implemented as in
#'  Tay et. al. (2020), or a function with parameters `z` and `theta`.
#'  
#'  `"original"` corresponds to
#'  
#'  `function(z, theta) mean(exp(z %*% theta)) * exp(- z %*% theta)`
#'@param theta `[NULL]` For development purposes only: If set to a matrix with
#'  dimensions ncol(z) x 1, no theta optimization step will take place and
#'  theta will be used for `weight_fun()` as is. In regular settings, this
#'  is not advisable.
#' @param t The initial step size for `theta` backtracking line search
#' (default value is 1).
#' @param a The factor by which step size is decreased in `theta`
#' backtracking line search (default value is 0.5).
#' @param thresh If the mean/median objective function does not decrease by at
#' least this factor, we terminate the optimization early. Default is 1e-4.
#' @param verbose If `TRUE`, prints information to console as model is
#' being fit. Default is `FALSE`.
#'
#' @return An object of class `"fwelnet"`.
#' \item{beta}{A `p x length(lambda)` matrix of coefficients.}
#' \item{theta}{Theta value, a `nfeaturevars x 1` matrix.}
#' \item{a0}{Intercept sequence of length `length(lambda)`.}
#' \item{lambda}{The actual sequence of `lambda` values used.}
#' \item{nzero}{The number of non-zero coefficients for each value of
#' `lambda`.}
#' \item{family}{Response type.}
#' \item{call}{The call that produced this object.}
#' \item{obj}{A `max_iter + 1` by `length(lambda)` matrix of
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
#' @importFrom stats median sd model.matrix
#' @importFrom glmnet glmnet
#' @export
fwelnet <- function(x, y, z, lambda = NULL, family = c("gaussian", "binomial", "cox"),
                    alpha = 1, standardize = TRUE, max_iter = 20, ave_mode = 1,
                    thresh_mode = 1, t = 1, a = 0.5, thresh = 1e-3,
                    weight_fun = "original",
                    theta = NULL,
                    verbose = FALSE) {
    this.call <- match.call()

    assert_number(alpha, lower = 1, upper = 1)
    assert_choice(family, choices = c("gaussian", "binomial", "cox"))

    # If we allow weight_fun to be an actual function, we can't rely on == 
    # for character comparison. Maybe pre-defining weight_funs and using
    # switch() for selection might be simpler.
    if (identical(weight_fun, "original")) {
      weight_fun <- function(z, theta) mean(exp(z %*% theta)) * exp(- z %*% theta)
    }

    # Allow data.frame/categorical input for survival without intercept
    # This is problematic though since it changes the number of columns in X
    # and z has to be pre-specified based on the "effective" ncol(x)
    if (is.data.frame(x)) {
      if (inherits(y, "Surv")) {
        x <- model.matrix(~ . -1, x)
      } else {
        x <- model.matrix(~ ., x)
      }
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
     glmfit <- glmnet::glmnet(x, y,
       alpha = alpha, family = family,
       standardize = FALSE
     )
     lambda <- glmfit$lambda
   } else {
     # Ensure lambda values are sorted in accordance with glmnet
     lambda <- as.double(rev(sort(lambda)))
     glmfit <- glmnet::glmnet(x, y,
       alpha = alpha, family = family,
       standardize = FALSE, lambda = lambda
     )
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
    iter <- 0    # no. of beta/theta minimizations

    # MT debug:
    if (!is.null(theta)) {
      # If theta is fixed, we still need to step through the optimization
      # machinery exactly once and selectively suppress optimization steps
      # Since fixing theta is not intended originally, this is rather clunky.
      # First, ensure we only iterate once
      # iter <- max_iter - 1
      max_iter <- iter + 1
      
      # Set flag for easier disabling of steps below
      skip_theta_optim <- TRUE
      
      # check if provided theta is even valid
      if (!is.matrix(theta)) theta <- as.matrix(theta) # in case of scalar input
      stopifnot(identical(dim(theta), c(as.integer(K), 1L)))
    } else {
      # Otherwise theta is initialized with zeros normally
      theta <- matrix(0, nrow = K, ncol = 1)
      skip_theta_optim <- FALSE
    }

    # store obj values across iterations (including original elastic net soln)
    obj_store <- matrix(NA, nrow = max_iter + 1, ncol = length(lambda))
    obj_store[1, ] <- objective_fn(x, y, z, beta, a0, theta, lambda, alpha, family)
    prev_iter_m_obj_value <- get_ave_obj(obj_store[1, ])
    
    # MT debug:: theta is scalar in our case, store it in a matrix per iteration step
    theta_store <- NULL
    theta_store <- matrix(NA_real_, nrow = max_iter + 1, ncol = K)
    theta_store[1, ] <- 0
  
    
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

        if (!skip_theta_optim) {
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
        }
        # OPTIMIZATION FOR BETA: ONE GLMNET STEP
        # pen.weights = mean(exp(z %*% theta)) * exp(- z %*% theta)
        pen.weights <- weight_fun(z, theta)
        fit <- glmnet::glmnet(x, y, family = family, standardize = FALSE,
                              alpha = alpha, lambda = lambda * mean(pen.weights),
                              penalty.factor = pen.weights)
        beta <- matrix(fit$beta, nrow = p)
        a0 <- fit$a0
        
        # MT debug: If we fix theta, we break early here since we don't need the rest
        if (skip_theta_optim) break

        # compute new objective and store it
        new_obj_value <- objective_fn(x, y, z, beta, a0, theta, lambda, alpha, family)
        new_m_obj_value <- get_ave_obj(new_obj_value)
        obj_store[iter + 1, ] <- new_obj_value
        
        # MT debug: Store theta if it's a scalar only
        theta_store[iter + 1, ] <- t(theta)

        
        if (verbose) cat(paste(ave_fn_name,
            "objective function value after beta minimization:",
            format(round(new_m_obj_value, 3), nsmall = 3)),
            fill = TRUE)

        # if mean/median obj function value not decreased enough, terminate early
        if ((prev_iter_m_obj_value - new_m_obj_value) /
            prev_iter_m_obj_value < thresh) {
            obj_store <- obj_store[1:(iter + 1), ]
            
            # MT debug: Trim theta if it's a scalar only
            theta_store <- theta_store[1:(iter + 1), ]
            
            
            if (verbose) cat(paste(ave_fn_name,
                "objective function value not decreased enough, stopping early\n"))
            break
        } else {
            prev_iter_m_obj_value <- new_m_obj_value
        }
    }

    # rescale intercept & beta values (beta will only change when standardize = TRUE)
    # is NULL for "cox" anyway
    if(!is.null(a0)) a0 <- a0 - colSums(beta * mx / sx)
    if (family == "gaussian") a0 <- a0 + y_mean
    beta <- beta / sx

    # count no. of non-zero coefficients for each model
    nzero <- colSums(beta != 0)

    out <- list(beta = beta, theta = theta, a0 = a0, lambda = lambda, nzero = nzero,
                family = family, call = this.call, obj = obj_store, theta_store = theta_store)
    class(out) <- "fwelnet"
    return(out)
}

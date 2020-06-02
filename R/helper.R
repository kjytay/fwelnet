msefun <- function(yhat,y) {
    (y - yhat)^2
}

binfun <- function(yhat, y) {
    - y * log(yhat) - (1 - y) * log(1 - yhat)
}


#' @importFrom graphics segments
# helper function for plotting CV curve: draws the error bars
error.bars <- function(x, upper, lower, width = 0.02, ...) {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

# gradient of the objective function w.r.t. theta
# beta: p x nlam matrix
# theta: K x 1 matrix
# lambda: vector of length nlam
# returns K x nlam matrix
gradient_fn <- function(x, y, z, beta, theta, lambda, alpha = 1) {
    p <- ncol(x)
    zT <- t(z)  # K x p
    exp_term <- exp(z %*% theta)  # p x 1
    beta_exp_term <- (alpha * abs(beta) + (1 - alpha) * beta^2 / 2) /
        drop(exp_term) # p x nlam
    scale(-zT %*% beta_exp_term * sum(exp_term) +
              zT %*% exp_term %*% matrix(colSums(beta_exp_term), nrow = 1),
          center = FALSE,
          scale = p / lambda)
}

# gradient of the objective function w.r.t. theta
# old version (superseded by gradient_fn): only works for single lambda value
# beta and theta should be matrices of column 1
# lambda should be a single real value
old_gradient_fn <- function(x, y, z, beta, theta, lambda, alpha = 1) {
    p <- ncol(x)
    zT <- t(z)
    exp_term <- exp(z %*% theta)
    beta_exp_term <- (alpha * abs(beta) + (1 - alpha) * beta^2 / 2) / exp_term
    grad_theta <- lambda / p * ( sum(beta_exp_term) * zT %*% exp_term -
                                     sum(exp_term) * zT %*% beta_exp_term )
    grad_theta
}

# penalty function
# beta: p x nlam matrix
# theta: K x 1 matrix
# lambda: vector of length nlam
penalty_fn <- function(z, beta, theta, lambda, alpha = 1) {
    p <- nrow(beta)
    exp_term <- drop(exp(z %*% theta))
    sum(exp_term) / p * lambda *
        colSums((alpha * abs(beta) + (1 - alpha) * beta^2 / 2 ) / exp_term)
}

# NLL function
# beta: p x nlam matrix
# a0: vector of length nlam
nll_fn <- function(x, y, beta, a0, family = "gaussian") {
    y <- drop(y)
    pred <- scale(x %*% beta, center = -a0, scale = FALSE)
    if (family == "gaussian") {
        nll <- colMeans((y - pred)^2) / 2
    } else if (family == "binomial") {
        nll <- colMeans(y * log(1 + exp(-pred)) + (1 - y) * log(1 + exp(pred)))
    } else {
        stop("Invalid value for family option")
    }
}

# objective function
# beta: p x nlam matrix
# theta: K x 1 matrix
# a0: vector of length nlam
# lambda: vector of length nlam
objective_fn <- function(x, y, z, beta, a0, theta, lambda, alpha = 1, family = "gaussian") {
    nll_fn(x, y, beta, a0, family) +
        penalty_fn(z, beta, theta, lambda, alpha)
}

# Reduced version of Binder et al (2009), basically just cox exponential model with normal predictors and equal event and censoring props
# Uses DGP as per Bender et al (2005)
#' Simple competing risk simulation
#' @param n,p Number of observations and number of features.
#' @param beta1,beta2 Coefficient vectors of the Cox-exponential model.
#' @param lambda1,lambda2 Baseline hazard constants of the models for event 1 and 2.
#' @param lambda_c Baseline hazard constant for the uninformative censoring distribution.
#' @param r Correlation between features
#' @export
#' 
#' @examples
#' sim_cr_coxexp()
sim_cr_coxexp <- function(n = 100, p = 4, beta1 = runif(p, -1, 1), beta2 = beta1, lambda1 = 0.1, lambda2 = lambda1, lambda_c = 0.1, r = 0) {
  
  checkmate::assert_int(n, lower = 10)
  checkmate::assert_int(p, lower = 2)
  checkmate::assert_numeric(beta1, any.missing = FALSE, len = p)
  checkmate::assert_numeric(beta2, any.missing = FALSE, len = p)
  lapply(list(lambda1, lambda2, lambda_c), \(x) {
    checkmate::assert_number(x, lower = 1e-8)
  })
  checkmate::assert_number(r, lower = 0, upper = 1)
  
  if (r == 0) {
    X <- matrix(rnorm(n * p), ncol = p, nrow = n)
  } else {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop("Package mvtnorm is required for simulating correlated data.")
    }
    X <- mvtnorm::rmvnorm(n = n, sigma = toeplitz(r^(0:(p - 1))))
  }
  
  # Linear predictors
  lp1 <- X %*% beta1
  lp2 <- X %*% beta2
  
  # Survival and censoring times
  Ti1 <- -log(runif(n)) / (lambda1 * exp(lp1))
  Ti2 <- -log(runif(n)) / (lambda2 * exp(lp2))
  # Default lambda_c = 0.1 should yield roughly 36% censored events.
  Ci <- -log(runif(n)) / lambda_c
  
  ti <- pmin(Ti1, Ti2, Ci)
  # Censored where neither event occurs before censoring time
  # if event before censoring time -> TRUE -> integer 1
  di <- as.integer(Ti1 <= Ci | Ti2 <= Ci)
  # Set status == 2 if Ti2 is observed
  di[which(Ti2 <= Ti1 & Ti2 <= Ci)] <- 2
  
  colnames(X) <- paste0("x", seq_len(p))
  xdat <- data.frame(time = ti, status = di, X)
  
  
  # separate event indicators for cause-specific fitting
  status_c1 <- status_c2 <- di
  status_c1[which(di == 2)] <- 0
  
  status_c2[which(di == 1)] <- 0
  status_c2[which(di == 2)] <- 1
  
  structure(list(
    data = xdat,
    Xmat = X,
    time = ti,
    status_c1 = status_c1,
    status_c2 = status_c2,
    beta1 = beta1,
    beta2 = beta2,
    lp1 = lp1, 
    lp2 = lp2),
    class = "sim_cr"
  )
}

print.sim_cr <- function(x, ...) {
  cat(sprintf("Simulated dataset with N = %d and p = %d\n", nrow(x$Xmat), ncol(x$Xmat)))
  cat("beta1:", sprintf("%.3f", x$beta1), "\n")
  cat("beta2:", sprintf("%.3f", x$beta2))
}

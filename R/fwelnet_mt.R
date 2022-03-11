#' Fit fwelnet for Multi-Task Learning
#'
#' Survival-focused implementation of the fwelnet-multitask algorithm described
#' in Algorithm 2 on page 13 (Tay et. al. 2020)
#'
#' @param data A data.frame or matrix holding predictors and outcome,
#' with outcome variables assumed to be named `"time"` and `"status"`.
#' @param causes Integer vector indicating causes, e.g. `1:2` for two causes.
#' @param mt_max_iter `[5]` number of mt-iterations to perform. Will break early
#' if no change in per-cause beta vector between iterations is detected.
#' If set to `0`, no `fwelnet` iteration will be performed and the returned
#' coefficients will be the result of fitting a cause-specific `glmnet`.
#' @param z_scale `[1]` Scalar for `z = abs(beta)` step.
#' @param z_method `["original"]` Either assign `z1` to be informed by `beta2`
#' of the current iteration (default behavior, as described in Algorithm 2
#' in Tay et. al. 2020), or `"aligned"` to have both `z1` and `z2` be informed
#' by `beta2` and `beta1` from the previous iteration step respectively.
#' @param alpha `[1]` Passed to [`glmnet()`] and [`fwelnet()`].
#' @param verbose Display informative message on the state of the mt fit.
#' @param ... Passed to [`fwelnet()`]
#'
#' @export
#' @importFrom survival Surv
#' @importFrom glmnet cv.glmnet
#' @return A `list` containing per-cause beta matrices for each iteration step
#'
fwelnet_mt_cox <- function(data, causes = 1:2,
                           mt_max_iter = 5,
                           z_scale = 1,
                           z_method = c("original", "aligned"),
                           alpha = 1, # pass to glmnet and fwelnet
                           verbose = FALSE, ...) {

  z_method <- match.arg(z_method)

  # data prep ---------------------------------------------------------------
  rowidx <- list()
  for (i in causes) {
    # hold row ids for events of each cause
    rowidx[[i]] <- which(data$status == i)
  }

  # Slow copy of data :(
  status_c1 <- data$status
  status_c2 <- data$status

  # Set events of cause 2 to censored for cause 1 dataset and vice versa
  status_c1[rowidx[[2]]] <- 0
  status_c2[rowidx[[1]]] <- 0
  # Status indicator for cause 2 should be 1 instead of 2
  status_c2[rowidx[[2]]] <- 1

  # Sanity check ----
  orig_status <- table(data$status)
  c1_status <- table(status_c1)
  c2_status <- table(status_c2)

  # All have same N
  stopifnot(all.equal(sum(orig_status), sum(c1_status), sum(c2_status)))

  # Events check out
  stopifnot(orig_status[["1"]] == c1_status[["1"]])
  stopifnot(orig_status[["2"]] == c2_status[["1"]])

  # Predictor matrix shared for all causes
  X <- data[, !(names(data) %in% c("time", "status"))]
  # Treatment-encode factors, ensure numeric design matrix without intercept
  X <- model.matrix(~ . -1, X)
  # time is shared between causes
  y1 <- survival::Surv(data$time, event = status_c1)
  y2 <- survival::Surv(data$time, event = status_c2)


  # fwelnet multi-task: Algorithm 2,  p. 13 ---------------------------------

  # Alg step 1) Initialize b1_0, b2_0 at lambda.min glmnet solution
  # for y_1, y_2 respectively
  gl1 <- glmnet::cv.glmnet(X, y1, family = "cox", alpha = alpha)
  gl2 <- glmnet::cv.glmnet(X, y2, family = "cox", alpha = alpha)

  b1_0 <- gl1$glmnet.fit$beta[, which(gl1$lambda == gl1$lambda.min)]
  b2_0 <- gl2$glmnet.fit$beta[, which(gl2$lambda == gl2$lambda.min)]

  # Hold betas in matrix, columns = iter+1, rows = p
  # beta1[j, k] holds cause-1 beta_j at iteration k+1 with k = 1 being glmnet solution
  beta1 <- matrix(NA_real_, ncol = mt_max_iter + 1, nrow = ncol(X), dimnames = list(dimnames(X)[[2]]))
  beta2 <- beta1

  # Initialize first col with glmnet fit
  beta1[, 1] <- as.numeric(b1_0)
  beta2[, 1] <- as.numeric(b2_0)

  # Alg step 2) init k = 0, while k < mt_max_iter:
  # go from 1 to max_iter+1 so we're doing max_iter steps and 1 can be 0.
  k <- 1

  while (k < mt_max_iter + 1) {

    if (verbose) message("k = ", k)

    # Alg step 2a)
    z2 <- z_scale * abs(beta1[, k, drop = FALSE])

    # Get betas from fwelnet fit, requires finding lambda.min via cv first
    fw2 <- cv.fwelnet(X, y2, z2, family = "cox", alpha = alpha, ...)
    beta2[, k + 1] <- fw2$glmfit$beta[, which(fw2$lambda == fw2$lambda.min)]

    # Alg step 2b)
    z1 <- switch (z_method,
      "original" = z_scale * abs(beta2[, k + 1, drop = FALSE]),
      "aligned"  = z_scale * abs(beta2[, k, drop = FALSE])
    )

    fw1 <- cv.fwelnet(X, y1, z1, family = "cox", alpha = alpha, ...)
    beta1[, k + 1] <- fw1$glmfit$beta[, which(fw1$lambda == fw1$lambda.min)]

    # Check beta differences, break if differences are 0
    beta1_diff <- beta1[, k] - beta1[, k + 1]
    beta2_diff <- beta2[, k] - beta2[, k + 1]

    beta1_stagnant <- isTRUE(all.equal(beta1_diff, rep(0, length(beta1_diff))))
    beta2_stagnant <- isTRUE(all.equal(beta2_diff, rep(0, length(beta2_diff))))

    if (verbose) {
      if (beta1_stagnant) message("No change in beta1 after k = ", k)
      if (beta2_stagnant) message("No change in beta2 after k = ", k)
    }

    if (beta1_stagnant & beta2_stagnant) {
      message("No change in beta{1,2} at k = ", k)
      break
    }

    if (verbose) {
      browser()

      message("Change in beta{1,2}: ")
      print(rbind(beta1_diff, beta2_diff), digits = 3)
    }

    # Increment k at the very last step
    k <- k + 1
  }

  # Returns betas and anything else that might be interesting
  list(
    beta1 = beta1,
    beta2 = beta2,
    # Check mt iterations later to assess reasonable values
    mt_iter = k - 1, # Adjust since k can't start at 0
    mt_max_iter = mt_max_iter,
    converged = (k < mt_max_iter),
    z_scale = z_scale,
    z_method = z_method
  )
}


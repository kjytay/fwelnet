#' Fit fwelnet for Multi-Task Learning
#'
#' Survival-focused implementation of the fwelnet-multitask algorithm described
#' in Algorithm 2 on page 13 (Tay et. al. 2020)
#'
#' @param data A data.frame or matrix holding predictors and outcome,
#' with outcome variables assumed to be named `"time"` and `"status"`.
#' @param causes (Unused for now) Integer vector indicating causes, e.g. `1:2` for two causes.
#' @param mt_max_iter `[5]` number of mt-iterations to perform. Will break early
#' if no change in per-cause beta vector between iterations is detected.
#' If set to `0`, no `fwelnet` iteration will be performed and the returned
#' coefficients will be the result of fitting a cause-specific `glmnet`.
#' @param z_method `("original")` Either assign `z1` to be informed by `beta2`
#' of the current iteration (default behavior, as described in Algorithm 2
#' in Tay et. al. 2020), or `"aligned"` to have both `z1` and `z2` be informed
#' by `beta2` and `beta1` from the previous iteration step respectively.
#' @param alpha `(1)` Passed to [`glmnet()`] and [`fwelnet()`].
#' @param stratify_by_status (`FALSE`) If `TRUE`, the internal cross-validation folds are sampled by the `status`
#'  variable, leading to folds with roughly equal proportions of events as in `data`.  
#'  This is helpful in cases where many censored observations in folds lead to errors during model fits.
#' @param standardize Passed to [`cv.fwelnet`] and [`cv.glmnet`] respectively.
#' @param verbose Display informative message on the state of the mt fit.
#' @param nfolds (`10`) Passed to [`cv.fwelnet`] and [`cv.glmnet`] respectively.
#' @param include_mt_beta_history `[FALSE]` If `TRUE`, output `list` includes
#'   components `beta1` and `beta2`, matrices of dimensions `p` x `mt_iter_max + 1`
#'   containing coefficient vectors for causes 1 and 2 for each multi-task iteration,
#'   with the first column corresponding to the original `cv.glmnet` solution.
#' @param ... Passed to [`fwelnet()`]
#' @inheritParams fwelnet
#'
#' @export
#' @importFrom survival Surv
#' @importFrom glmnet cv.glmnet
#' @return An object of class `cooper` with `cv.fwelnet` objects of the final iteration.
#' If `include_mt_beta_history = TRUE`, contains per-cause beta matrices for each iteration step.
#'
fwelnet_mt_cox <- function(data, 
                           causes = 1:2, # Mostly unused for now
                           mt_max_iter = 5,
                           z_method = "original",
                           stratify_by_status = FALSE,
                           alpha = 1, # pass to glmnet and fwelnet
                           standardize = TRUE,
                           verbose = FALSE, t = 1, a = 0.5, 
                           thresh = 1e-3,
                           nfolds = 10,
                           include_mt_beta_history = FALSE,
                           ...) {
  
  # Convert to data.frame just in case it's a data.table
  data <- as.data.frame(data)
  assert_data_frame(data, any.missing = FALSE)
  assert_integer(causes, any.missing = FALSE, min.len = 2)
  assert_int(mt_max_iter, lower = 1)
  assert_choice(z_method, choices = c("original", "aligned"))
  assert_numeric(alpha, any.missing = FALSE, lower = 0, upper = 1, len = 1)
  assert_logical(verbose, len = 1)
  # FIXME: Rename t variable to avoid clash with base::t
  assert_numeric(t, any.missing = FALSE, len = 1)
  assert_numeric(alpha, any.missing = FALSE, lower = 0, upper = 1, len = 1)
  assert_numeric(thresh, any.missing = FALSE, lower = 0, len = 1)
  assert_int(nfolds, lower = 3)
  assert_logical(include_mt_beta_history, len = 1)

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
  X <- model.matrix(~ ., X)
  # Dropping intercept also drops attributes, not sure of this is a problem
  X <- X[, -1, drop = FALSE] 
  
  # time is shared between causes
  y_list <- list(
    survival::Surv(data$time, event = status_c1),
    survival::Surv(data$time, event = status_c2)
  )
  
  # Stratify by status, generating per-observation fold IDs and passing those to initial
  # cv.glmnet and cv.fwelnet later on
  y1foldids <- NULL
  y2foldids <- NULL
  
  if (stratify_by_status) {
    y1df <- data.frame(status = status_c1)
    y2df <- data.frame(status = status_c2)
    
    y1foldids <- stratified_cv_folds(y1df, nfolds = nfolds)[["fold"]]
    y2foldids <- stratified_cv_folds(y2df, nfolds = nfolds)[["fold"]]
  }

  # fwelnet multi-task: Algorithm 2,  p. 13 ---------------------------------

  # Alg step 1) Initialize b1_0, b2_0 at lambda.min glmnet solution
  # for y_1, y_2 respectively
  gl1 <- glmnet::cv.glmnet(X, y_list[[1]], family = "cox", alpha = alpha, foldid = y1foldids, standardize = standardize)
  gl2 <- glmnet::cv.glmnet(X, y_list[[2]], family = "cox", alpha = alpha, foldid = y2foldids, standardize = standardize)

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
    
    fw_cv_list <- vector(mode = "list", length = 2)
    
    # Stratify by status, generating per-observation fold IDs and passing those cv.fwelnet fits
    # using different fold assignments for each iteration, maybe should assign folds once initially and reuse here?
    if (stratify_by_status) {
      y1foldids <- stratified_cv_folds(y1df, nfolds = nfolds)[["fold"]]
      y2foldids <- stratified_cv_folds(y2df, nfolds = nfolds)[["fold"]]
    }

    # Alg step 2a)
    z2 <- abs(beta1[, k, drop = FALSE])

    # Get betas from fwelnet fit, requires finding lambda.min via cv first
    fw_cv_list[[2]] <- cv.fwelnet(X, y_list[[2]], z2, family = "cox", alpha = alpha, t = t, a = a, thresh = thresh, foldid = y2foldids, standardize = standardize, ...)
    beta2[, k + 1] <- fw_cv_list[[2]]$glmfit$beta[, which(fw_cv_list[[2]]$lambda == fw_cv_list[[2]]$lambda.min)]

    # Alg step 2b)
    z1 <- switch (z_method,
      "original" = abs(beta2[, k + 1, drop = FALSE]),
      "aligned"  = abs(beta2[, k, drop = FALSE]),
      stop("z_method ", z_method, " not known")
    )

    fw_cv_list[[1]] <- cv.fwelnet(X, y_list[[1]], z1, family = "cox", alpha = alpha, t = t, a = a, thresh = thresh, foldid = y1foldids, standardize = standardize, ...)
    beta1[, k + 1] <- fw_cv_list[[1]]$glmfit$beta[, which(fw_cv_list[[1]]$lambda == fw_cv_list[[1]]$lambda.min)]

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
      message("Change in beta{1,2}: ")
      print(rbind(beta1_diff, beta2_diff), digits = 3)
    }

    # Increment k at the very last step
    k <- k + 1
  }
  
  # get baseline hazards using x = 0 as reference, possibly bad idea according to ?basehaz
  newdata_zero <- X[1, ] * 0
  
  basehazards <- lapply(causes, \(e) {
    # Accesses the glmnet fit returned with the fwelnet fit with the cv.fwelnet object (naming is bad here yes)
    bhraw <- survival::survfit(
      fw_cv_list[[e]]$glmfit$glmfit, s = fw_cv_list[[e]]$lambda.min, 
      x = X,
      y = y_list[[e]],
      newx = newdata_zero,
      se.fit = FALSE
    )
    
    data.frame(
      event = e,
      time = bhraw$time,
      n.risk = bhraw$n.risk,
      surv = bhraw$surv,
      cumhaz = bhraw$cumhaz,
      hazard = diff(c(0, bhraw$cumhaz))
    )
  })
  

  # Returns betas and anything else that might be interesting
  ret <- list(
    # final fwelnet fit objects for later predictions
    fwfit1 = fw_cv_list[[1]], # To be removed
    fwfit2 = fw_cv_list[[2]], # To be removed
    # 0-th step glmnet solutions
    glmfit1 = gl1, # To be removed
    glmfit2 = gl2, # To be removed
    # Collected versions of same objects
    initial_fits = list(gl1, gl2),
    fwelfits = fw_cv_list,
    # model matrix
    x = X, y = y_list,
    predictors = setdiff(colnames(data), c("time", "status")),
    # Baseline hazards
    basehazards = basehazards,
    # Check mt iterations later to assess reasonable values
    mt_iter = k - 1, # Adjust since k can't start at 0
    mt_max_iter = mt_max_iter,
    converged = ((k - 1) < mt_max_iter)
  )
  
  if (include_mt_beta_history) {
    #  History of betas
    ret$beta1 <- beta1
    ret$beta2 <- beta2
  }
  
  class(ret) <- c("cooper")
  ret
}

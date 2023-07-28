library(riskRegression)
library(survival)
library(fwelnet)
library(ggplot2)

stopifnot(
  "Update fwelnet: remotes::install_github('jemus42/fwelnet')" = packageVersion("fwelnet") >= package_version("0.1.0.9006")
)

# Reduced version of Bindert et al dgp, basically just cox exponential model mit normal predictors and equal event and censoring props
simple_cr <- function(n = 100, p = 4, beta1 = rep(0.5, p), beta2 = beta1, lambda1 = 0.1, lambda2 = lambda1, lambda_c = 0.1, r = 0) {
  
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
    if (!requireNamespace("mvtnorm")) {
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
  
  list(
    data = xdat,
    Xmat = X,
    time = ti,
    status_c1 = status_c1,
    status_c2 = status_c2,
    beta1 = beta1,
    beta2 = beta2,
    lp1 = lp1, 
    lp2 = lp2
  )
}

# Simple cleanup for riskRegression::Score() output for ggplot usage
cleanup_score <- function(scores) {
  result <- scores$AUC$score[scores$Brier$score, on = c("model", "times")]
  result <- data.table::melt(result , id.vars = c("model", "times"),
                             value.name = "score", variable.name = "metric",
                             meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
  # Exclude some superfluous output
  result <- result[!(is.na(score) & model == "Null model"), ]
  result <- result[!(metric == "IPA" & model == "Null model"), ]
  result
}


# Simulate data ---------------------------------------------------------------------------------------------------
set.seed(123)
instance <- simple_cr(
  n = 1000, p = 4, 
  beta1 = c(0.4, 0.3, -0.2, 0.5), 
  beta2 = c(0.5, -0.3, -0.7, -0.4),
  # Make event 1 a bit rarer
  lambda1 = 0.05, lambda2 = 0.1,
  lambda_c = 0.1
)
xtrain <- instance$data
table(xtrain$status)

survfit(Surv(instance$time, instance$status_c1) ~ 1) |>
  plot()

# Models: full model sees all features, reduced model only sees x1 and x2
# Focusing on event 1

# Fit Cooper ------------------------------------------------------------------------------------------------------
cooper_full <- cooper(
  data = xtrain, mt_max_iter = 5, stratify_by_status = TRUE, alpha = 1, verbose = TRUE, 
  #t = 100, thresh = 1e-7,
  epsilon1 = 1e-15, epsilon2 = 1e-15
)

cbind(
  cooper = coef(cooper_full, event = 1),
  truth = instance$beta1
)

cooper_reduced <- cooper(
  data = xtrain[, c("time", "status", "x1", "x2")], mt_max_iter = 5, stratify_by_status = TRUE, verbose = TRUE, epsilon1 = 1e-15, epsilon2 = 1e-15
)

# Fit vanilla glmnet ----------------------------------------------------------------------------------------------
rr_glmnet_full <- GLMnet(
  formula = reformulate(cooper_full$predictors, response = "Surv(time, status)"), 
  data = within(xtrain, status[status == 2] <- 0),
  lambda = cooper_full$initial_fits[[1]]$lambda.min,
  cv = FALSE,
  standardize = TRUE, alpha = 1
)

rr_glmnet_reduced <- GLMnet(
  formula = reformulate(cooper_reduced$predictors, response = "Surv(time, status)"), 
  data = within(xtrain, status[status == 2] <- 0),
  lambda = cooper_reduced$initial_fits[[1]]$lambda.min,
  cv = FALSE,
  standardize = TRUE, alpha = 1
)


# Fit regular cause-specific cox ----------------------------------------------------------------------------------
csc_full <- CSC(reformulate(cooper_full$predictors, response = "Hist(time, status)"), data = xtrain)
csc_reduced <- CSC(reformulate(cooper_reduced$predictors, response = "Hist(time, status)"), data = xtrain)


# Evaluate --------------------------------------------------------------------------------------------------------
# Use 10% through 80% observed event 1 times in 10% steps
eval_times <- quantile(xtrain$time[xtrain$status == 1], probs = seq(0.1, 0.8, .1), type = 2, names = FALSE)

scores <- Score(
  list(cooper_full = cooper_full, cooper_reduced = cooper_reduced, 
       rr_glmnet_full = rr_glmnet_full, rr_glmnet_reduced = rr_glmnet_reduced, 
       csc_full = csc_full, csc_reduced = csc_reduced),
  formula = Hist(time, status) ~ 1,
  data = xtrain,
  cause = 1,
  metrics = c("Brier", "AUC"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
) |>
  cleanup_score()


# Plot facetted to compare full and reduced
scores |>
  dplyr::filter(metric %in% c("Brier", "AUC")) |>
  dplyr::mutate(
    model_spec = stringr::str_extract(model, "full|reduced"),
    model = stringr::str_remove(model, "_(full|reduced)")
  ) |>
  ggplot(aes(x = times, y = score, color = model_spec, fill = model_spec)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), cols = vars(model), scales = "free", switch = "y") +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Full vs reduced models",
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Plot facetted to compare cooper and glmnet
scores |>
  dplyr::filter(metric %in% c("Brier", "AUC")) |>
  dplyr::mutate(
    model_spec = stringr::str_extract(model, "full|reduced"),
    model = stringr::str_remove(model, "_(full|reduced)")
  ) |>
  dplyr::filter(model_spec %in% c("full", NA)) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free", switch = "y") +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  #scale_color_brewer(palette = "Dark2", breaks = c("cooper", "rr_glmnet", "csc", "Null model")) +
  scale_color_discrete(breaks = c("cooper", "rr_glmnet", "csc", "Null model")) +
  labs(
    title = "Full (correctly specified) models",
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


# Plot all in one, not very easy to compare
# scores |>
#   dplyr::filter(metric %in% c("Brier", "AUC")) |>
#   ggplot(aes(x = times, y = score, color = model, fill = model)) +
#   #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
#   facet_grid(rows = vars(metric), scales = "free", switch = "y") +
#   geom_line() +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_brewer(palette = "Dark2") +
#   labs(
#     x = "Event Time (t)", y = "Metric",
#     color = NULL, fill = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")

# Compare coeffcients ---------------------------------------------------------------------------------------------

# for reduced models, coefs will get recycled, because
# coef(cooper_full, event = 1) and coef(cooper_reduced, event = 1) have different lengths obviously
# We fix that inelegantly
coef_df <- data.frame(
  x = cooper_full$predictors,
  cooper_full = coef(cooper_full, event = 1),
  cooper_reduced = c(coef(cooper_reduced, event = 1), 0, 0),
  rr_glmnet_full = as.vector(coef(rr_glmnet_full$fit, event = 1)),
  rr_glmnet_reduced = c(as.vector(coef(rr_glmnet_reduced$fit, event = 1)), 0, 0),
  csc_full = coef(csc_full$models$`Cause 1`),
  csc_reduced = c(coef(csc_reduced$models$`Cause 1`), 0, 0),
  truth_full = instance$beta1,
  truth_reduced = instance$beta1,
  row.names = NULL
)

coef_df

coef_df |>
  tidyr::pivot_longer(cols = -1, names_to = "model", values_to = "coef") |>
  dplyr::mutate(
    model_spec = stringr::str_extract(model, "full|reduced"),
    #model_spec = ifelse(model == "truth", "", model_spec),
    model = stringr::str_remove(model, "_(full|reduced)"),
  ) |>
  tidyr::pivot_wider(id_cols = c("model_spec", "x"), names_from = "model", values_from = "coef")
  # tidyr::pivot_wider(id_cols = c("model_spec", "model"), names_from = "x", values_from = "coef")
  # tidyr::pivot_wider(id_cols = c("model", "model_spec"), names_from = "x", values_from = "coef")


X_full <- instance$Xmat
X_reduced <- instance$Xmat[, c("x1", "x2")]

linpred_df <- data.frame(
  truth = instance$lp1,
  cooper_full = as.numeric(X_full %*% coef(cooper_full, event = 1)),
  cooper_reduced = as.numeric(X_reduced %*% coef(cooper_reduced, event = 1)),
  rr_glmnet_full = as.numeric(X_full %*% as.vector(coef(rr_glmnet_full$fit, event = 1))),
  rr_glmnet_reduced = as.numeric(X_reduced %*% as.vector(coef(rr_glmnet_reduced$fit, event = 1))),
  csc_full = X_full %*% coef(csc_full$models$`Cause 1`),
  csc_reduced = X_reduced %*% coef(csc_reduced$models$`Cause 1`),
  row.names = NULL
)

linpred_df |>
  dplyr::mutate(dplyr::across(-"truth", \(x) truth - x)) |>
  dplyr::select(-"truth") |>
  tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "model", values_to = "lp_error") |>
  dplyr::mutate(
    model_spec = stringr::str_extract(model, "full|reduced"),
    #model_spec = ifelse(model == "truth", "", model_spec),
    model = stringr::str_remove(model, "_(full|reduced)"),
  ) |>
  ggplot(aes(x = lp_error, color = model_spec, fill = model_spec)) +
  facet_wrap(vars(model), ncol = 1) +
  geom_density(alpha = 1/3) +
  labs(
    title = "Linear Predictor Estimates",
    subtitle = "Error as [(true lp) - (estimated lp)] of currect (full) and incomplete (reduced) models",
    x = expression(eta-hat(eta))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# Fit glmnet / fwelnet directly on single-cause subset ----------------------------------------------------------------

p <- ncol(instance$Xmat)
# Constant z will result in no weighting, equal to glmnet (no matter if all 1 or all 2 or all 0...)
z_ones <- matrix(rep(1, p))
# Some signal, doesn't matter which direction, just that it's not constant to see if it has an effect
z_signal <- matrix(abs(rnorm(p, sd = 5)))

set.seed(3)
fwelnet_no_z <- cv.fwelnet(instance$Xmat, y = Surv(instance$time, instance$status_c1), z = z_ones, family = "cox")
# Save lambda sequence to ensure identical lambdas across models
lambda_seq <- fwelnet_no_z$lambda

set.seed(3)
fwelnet_signal_z <- cv.fwelnet(instance$Xmat, y = Surv(instance$time, instance$status_c1), z = z_signal, family = "cox", lambda = lambda_seq)

set.seed(3)
# Fit cv.glmnet with same lambda as 
glmnet_fit <- cv.glmnet(instance$Xmat, y = Surv(instance$time, instance$status_c1), family = "cox", lambda = lambda_seq)

# Expecting the to be idential
identical(fwelnet_no_z$lambda, glmnet_fit$lambda)
waldo::compare(fwelnet_no_z$lambda.min, glmnet_fit$lambda.min)

# Comparing stored glmnet object to "vanilla" glmnet, expecting no meaningful differences
waldo::compare(fwelnet_no_z$glmfit$glmfit, glmnet_fit$glmnet.fit)

# Coefficient matrices differ though, but not by much
waldo::compare(Matrix::Matrix(fwelnet_no_z$glmfit$beta, sparse = TRUE), glmnet_fit$glmnet.fit$beta, max_diffs = 6)
# Difference negligible (1e-16)
as.numeric(Matrix::Matrix(fwelnet_no_z$glmfit$beta, sparse = TRUE) - glmnet_fit$glmnet.fit$beta) |> 
  summary()

# Fwelnet with and without informative z should differ
identical(fwelnet_no_z$lambda, fwelnet_signal_z$lambda)
waldo::compare(fwelnet_no_z$lambda.min, fwelnet_signal_z$lambda.min)

# Comparing stored glmnet object to "vanilla" glmnet, expecting no meaningful differences
waldo::compare(fwelnet_no_z$glmfit$glmfit, fwelnet_signal_z$glmfit$glmfit)

# Coefficient matrices differ enough to assume z has influence, not a big one but not the point here
waldo::compare(fwelnet_no_z$glmfit$beta, fwelnet_signal_z$glmfit$beta, max_diffs = 1)

as.numeric(fwelnet_no_z$glmfit$beta - fwelnet_signal_z$glmfit$beta) |> 
  summary()

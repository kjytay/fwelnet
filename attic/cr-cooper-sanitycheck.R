library(riskRegression)
library(survival)
library(fwelnet)
library(ggplot2)

stopifnot(
  "Update fwelnet: remotes::install_github('jemus42/fwelnet')" = packageVersion("fwelnet") >= package_version("0.1.0.9005")
)

# Reduced version of Bindert et al dgp, basically just cox exponential model mit normal predictors and equal event and censoring props
simple_cr <- function(n = 100, p = 4, beta1 = rep(0.5, p), beta2 = beta1, lambda1 = 0.1, lambda2 = lambda1, lambda_c = 0.1) {
  X <- matrix(rnorm(n * p), ncol = p, nrow = n)
  
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
  
  list(
    data = xdat,
    Xmat = X,
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
instance <- simple_cr(n = 1000, p = 4)
xtrain <- instance$data
table(xtrain$status)


# Fit Cooper ------------------------------------------------------------------------------------------------------
fitfull <- cooper(
  data = xtrain, mt_max_iter = 3, stratify_by_status = TRUE, alpha = 1, standardize = TRUE
)

fitreduced <- cooper(
  data = xtrain[, c("time", "status", "x1", "x2")], mt_max_iter = 3, stratify_by_status = TRUE, alpha = 1, standardize = TRUE
)


# Fit vanilla glmnet ----------------------------------------------------------------------------------------------
rr_glmnet_full <- GLMnet(
  formula = reformulate(fitfull$predictors, response = "Surv(time, status)"), 
  data = within(xtrain, status[status == 2] <- 0),
  lambda = fitfull$initial_fits[[1]]$lambda.min,
  cv = FALSE,
  standardize = TRUE, alpha = 1
)

rr_glmnet_reduced <- GLMnet(
  formula = reformulate(fitreduced$predictors, response = "Surv(time, status)"), 
  data = within(xtrain, status[status == 2] <- 0),
  lambda = fitreduced$initial_fits[[1]]$lambda.min,
  cv = FALSE,
  standardize = TRUE, alpha = 1
)


# Fit regular cause-specific cox ----------------------------------------------------------------------------------
csc_full <- CSC(reformulate(fitfull$predictors, response = "Hist(time, status)"), data = xtrain)
csc_reduced <- CSC(reformulate(fitreduced$predictors, response = "Hist(time, status)"), data = xtrain)


# Evaluate --------------------------------------------------------------------------------------------------------
# Use 10% through 80% observed event 1 times in 10% steps
eval_times <- quantile(xtrain$time[xtrain$status == 1], probs = seq(0.1, 0.8, .1), type = 2, names = FALSE)

scores <- Score(
  list(cooper_full = fitfull, cooper_reduced = fitreduced, 
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

# Plot all in one, not very easy to compare
scores |>
  dplyr::filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free", switch = "y") +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Plot facetted
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
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


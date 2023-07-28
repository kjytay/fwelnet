library(survival)
library(riskRegression)
library(fwelnet)

data(pbc, package = "survival")
pbc <- na.omit(pbc)
xtrain <- pbc[1:270, -1]
xtest <- pbc[271:276, -1]

eval_times <- quantile(xtrain$time, probs = seq(0.1, 0.6, .1), type = 2, names = FALSE)

set.seed(12)
cooperfit <- fwelnet_mt_cox(
  xtrain, mt_max_iter = 3,
  stratify_by_status = TRUE,
  standardize = TRUE # FALSE currently required to ensure correct results
)

csc <- CSC(Hist(time, status) ~ ., data = xtrain, cause = 1)

predictRisk(csc, xtest, cause = 1, times = eval_times)

xtest_c1 <- within(xtrain, status[status == 2] <- 0)
rr_glmnet <- GLMnet(
  formula = reformulate(cooperfit$predictors, response = "Surv(time, status)"), 
  data = xtest_c1,
  lambda = cooperfit$initial_fits[[1]]$lambda.min,
  cv = FALSE,
  standardize = FALSE, alpha = 1
)

scores <- Score(
  list(cooper = cooperfit, csc = csc, rr_glmnet = rr_glmnet),
  formula = Hist(time, status) ~ 1,
  data = xtrain,
  cause = 1,
  metrics = c("Brier", "AUC"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
)
scores$Brier$score
scores$AUC$score

debugonce(riskRegression:::Score.list)
debugonce(riskRegression:::getResponse)
debugonce(fwelnet:::get_abs_risk)




result <- scores$AUC$score[scores$Brier$score, on = c("model", "times")]
result <- data.table::melt(result , id.vars = c("model", "times"),
                           value.name = "score", variable.name = "metric",
                           meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
# Exclude some superfluous output
result <- result[!(is.na(score) & model == "Null model"), ]
result <- result[!(metric == "IPA" & model == "Null model"), ]
result

library(ggplot2)

result |>
  subset(metric == "Brier") |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

result |>
  subset(metric == "AUC") |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



# inital glmnet? --------------------------------------------------------------------------------------------------

rr_glmnet$fit$beta
rr_glmnet$fit$lambda

cooperfit$initial_fits[[1]]$glmnet.fit$beta[cooperfit$initial_fits[[1]]$lambda == cooperfit$initial_fits[[1]]$lambda.min]
waldo::compare(
  coef(cooperfit$initial_fits[[1]], s = "lambda.min"),
  coef(rr_glmnet$fit)
)



# compare coefs ---------------------------------------------------------------------------------------------------

library(survival)
library(riskRegression)
library(fwelnet)

data(pbc, package = "survival")
pbc <- na.omit(pbc)
xtrain <- pbc[1:270, -1]
xtest <- pbc[271:276, -1]

eval_times <- quantile(xtrain$time, probs = seq(0.1, 0.6, .1), type = 2, names = FALSE)

set.seed(12)
cooperfit <- fwelnet_mt_cox(
  xtrain, mt_max_iter = 3,
  stratify_by_status = TRUE,
  standardize = TRUE # FALSE currently required to ensure correct results
)

cvfwelnet <- cooperfit$fwelfits[[1]]

lambda_idx <- which(cvfwelnet$lambda.min == cvfwelnet$lambda)
cvfwelnet$glmfit$beta[, lambda_idx]

fwelnet_beta <- cvfwelnet$glmfit$beta[, lambda_idx]
glmnet_beta <- cvfwelnet$glmfit$glmfit$beta[, lambda_idx]

cbind(fwel = fwelnet_beta, glmn = glmnet_beta) |> round(3)

sx <- apply(cooperfit$x, 2, sd) * sqrt((nrow(cooperfit$x) - 1) / nrow(cooperfit$x))
mx <- colMeans(cooperfit$x)


glmnet_beta / sx



eval_times <- quantile(xtrain$time, probs = seq(0.1, 0.6, .1), type = 2, names = FALSE)

set.seed(12)
cooperfit <- fwelnet_mt_cox(
  xtrain, mt_max_iter = 3,
  stratify_by_status = TRUE,
  standardize = TRUE # FALSE currently required to ensure correct results
)

set.seed(12)
cooperfit_nonst <- fwelnet_mt_cox(
  xtrain, mt_max_iter = 3,
  stratify_by_status = TRUE,
  standardize = FALSE # FALSE currently required to ensure correct results
)


lambda_idx <- which(cooperfit$fwelfits[[1]]$lambda.min == cooperfit$fwelfits[[1]]$lambda)

xn <- dimnames(cooperfit$fwelfits[[1]]$glmfit$glmfit$beta)
xn2 <- dimnames(cooperfit$fwelfits[[2]]$glmfit$glmfit$beta)

identical(xn, xn2)

# fwelbeta <- Matrix::Matrix(cooperfit$fwelfits[[1]]$glmfit$beta, dimnames = xn)
# 
# fwelbeta[, lambda_idx]
# cooperfit$fwelfits[[1]]$glmfit$glmfit$beta[, lambda_idx]

cooperfit$fwelfits[[1]]$glmfit$glmfit$beta <- Matrix::Matrix(cooperfit$fwelfits[[1]]$glmfit$beta, dimnames = xn)
cooperfit$fwelfits[[2]]$glmfit$glmfit$beta <- Matrix::Matrix(cooperfit$fwelfits[[2]]$glmfit$beta, dimnames = xn2)

scores <- Score(
  list(cooper = cooperfit, cooper_nonst = cooperfit_nonst),
  formula = Hist(time, status) ~ 1,
  data = xtrain,
  cause = 1,
  metrics = c("Brier", "AUC"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
)

cleanup_score(scores) |>
  dplyr::filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

debugonce(fwelnet:::predict.cooper)
debugonce(fwelnet:::get_abs_risk)


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

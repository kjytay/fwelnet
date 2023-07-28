library(riskRegression)
# data(Melanoma)
# 
# melanoma <- subset(Melanoma, select = -3)
# testdat <- melanoma[1:5, -c(1:2)]
# #testdat <- model.matrix(~ ., data = melanoma[1:5, -c(1:2)])[, -1]
# set.seed(12)
# cooperfit <- fwelnet_mt_cox(
#   melanoma, mt_max_iter = 3,
#   standardize = FALSE # Currently required to ensure correct results
# )

library(fwelnet)
data(pbc, package = "survival")
pbc <- na.omit(pbc)
xtrain <- pbc[1:270, -1]
xtest <- pbc[271:276, -1]

set.seed(12)
cooperfit <- fwelnet_mt_cox(
  xtrain, mt_max_iter = 3,
  stratify_by_status = TRUE,
  standardize = FALSE # Currently required to ensure correct results
)

predict(cooperfit, xnew = xtest, type = "risk", event = 1)
predict(cooperfit, xnew = xtest, times = c(1180.5, 1788), type = "absrisk", event = 1)

predict(cooperfit, xnew = xtest, times = c(1180.5, 1788), type = "absrisk", event = 1, use_initial_fit = TRUE)

if (requireNamespace("riskRegression")) {
  library(riskRegression)
  predictRisk(cooperfit, xtest, cause = 1, times = c(1180.5, 1788))
}

survfit(Surv(time, status) ~ 1, data = within(pbc, status[status == 2] <- 0)) |>
  plot()

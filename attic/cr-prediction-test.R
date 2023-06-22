# Prediction from fwelnet
library(survival)
library(fwelnet)
library(riskRegression)
library(data.table)
library(ggplot2)

if (!requireNamespace("pem.xgb")) {
  pak::pak("adibender/pem.xgb")
}

set.seed(13)
# Simulate some CR data -------------------------------------------------------------------------------------------
n <- 500 
# create data set with covariates
xdf <- tibble::tibble(
  x0 = sample(c(-1,1), n, .3),
  x1 = runif(n, -3, 3),
  x2 = runif(n, -3, 3),
  x3 = runif(n, -3, 3))
xdf2 <- mvtnorm::rmvnorm(n = nrow(xdf), mean = rep(0, 10))
# noise variables
colnames(xdf2) <- paste0("x", 4:(ncol(xdf2) + 3))
xdf <- cbind(xdf, xdf2)

sim_wrapper_cr <- function(formula, data, time_grid  = seq(0, 10, by = 0.1)) {
  # baseline hazard
  instance <- pem.xgb::sim_pexp_cr(
    formula = formula,
    data    = data,
    cut     = time_grid
  )
  instance$status <- instance$type
  instance$type <- NULL
  # add censoring
  cens_times <- runif(nrow(instance), 0, 20)
  cens <- instance$time > cens_times
  instance$time <- pmin(instance$time, cens_times)
  instance$status[cens] <- 0
  instance$status[instance$time == 10] <- 0
  instance <- as.data.frame(instance)
  instance$id <- instance$hazard1 <- instance$hazard2 <- NULL
  instance$x0 <- as.factor(instance$x0)
  
  instance
  
}

sim1 <- sim_wrapper_cr(
  formula = ~ -3.5 + 2*dgamma(t, 8, 2) + 0.75 * x1 + 0.6 * x2| -3.5 + 2 * dgamma(t, 8, 2) + 0.5 * x1 + 0.5 * x2 + 0.25 * x3,
  data = xdf
)

# from ?basehaz_cpp: WARNING stoptimes status eXb and strata must be sorted by strata, stoptimes, and status
data.table::setorder(sim1, time, status)

eval_times <- quantile(sim1$time[sim1$status == 1], prob = seq(.1, .7, .05), type = 2, names = FALSE)
event_times_c1 <- sim1$time[sim1$status == 1]

# Cause-specific cox via survival::coxph --------------------------------------------------------------------------

sim1_c1 <- sim1
sim1_c2 <- sim1

sim1_c1$status[sim1_c1$status == 2] <- 0
sim1_c2$status[sim1_c2$status == 1] <- 0
sim1_c2$status[sim1_c2$status == 2] <- 1

coxph_c1 <- survival::coxph(Surv(time, status) ~ ., data = sim1_c1, x = TRUE)
coxph_c2 <- survival::coxph(Surv(time, status) ~ ., data = sim1_c2, x = TRUE)

coxph_coefs <- cbind(c1 = coef(coxph_c1), c2 = coef(coxph_c2))

# riskRegression CSC for comparison
rr_csc <- CSC(Hist(time, status) ~ ., data = sim1)

csc_coefs <- cbind(
  csc1 = rr_csc$models$`Cause 1`$coefficients,
  csc2 = rr_csc$models$`Cause 2`$coefficients
)

cbind(coxph_coefs, csc_coefs) |>
  round(3)

scores <- Score(
  list(rr_csc = rr_csc, coxph = coxph_c1), 
  formula = Hist(time, status) ~ 1, 
  data = sim1, 
  times = eval_times
)

scores$Brier$score |>
  ggplot(aes(x = times, y = Brier, color = model, fill = model)) +
  geom_line()


# Baseline survs ---------------------------------------------------------------------------------------------------

# centered = FALSE does not use mean of x as reference, makes more sense for binary vars
basehaz_survival <- basehaz(coxph_c1, centered = FALSE)

basehaz_rr <- baseHaz_cpp(
  starttimes = rep(0, nrow(sim1)),
  stoptimes = sim1$time,
  status = sim1$status,
  eXb = exp(coxph_c1$linear.predictors),
  predtimes = eval_times,
  nPatients = nrow(sim1),
  emaxtimes = max(sim1$time),
  cause = 1,
  strata = rep(1, nrow(sim1)),
  nStrata = 1,
  Efron = FALSE
)
basehaz_rr$strata <- NULL

basehaz_survival <- as.data.table(basehaz_survival[basehaz_survival$time %in% eval_times, ])
setnames(basehaz_survival, new = c("cumhazard_survival", "time"))

basehaz_rr <- as.data.table(basehaz_rr)
setnames(basehaz_rr, new = c("time", "hazard_rr", "cumhazard_rr"))

# basehaz does not provide same time points
cumhazards <- merge(basehaz_rr, basehaz_survival, by = "time")

# Calculate baseline hazard via survfit (what basehaz does)
newdata_zero <- matrix(0, ncol = 14, nrow = 1) |>
  as.data.frame() |>
  setNames(paste0("x", 0:13)) |>
  mutate(x0 = factor(-1))

basehaz_survfit0 <- survfit(coxph_c1, newdata = newdata_zero)

basehaz_survfit <- data.table(time = basehaz_survfit0$time, cumhazard_survfit = basehaz_survfit0$cumhaz, s0_survfit_orig = basehaz_survfit0$surv)
basehaz_survfit <- basehaz_survfit[basehaz_survfit$time %in% eval_times, ]
cumhazards <- merge(cumhazards, basehaz_survfit, by = "time")

# See if it checks out
cumhazards |>
  mutate(
    s0_rr = exp(-cumhazard_rr),
    s0_survfit = exp(-cumhazard_survfit)
  )

# Manual basehazard? ----------------------------------------------------------------------------------------------

xtemp1 <- sim1 |>
  select(time, status) |>
  mutate(
    elp1 = exp(coxph_c1$x %*% coxph_c1$coefficients),
    elp2 = exp(coxph_c2$x %*% coxph_c2$coefficients)
  )

# FALSE since $linear.predictors are centered, not sure if better to take those?
# Difference is constant
(coxph_c1$x %*% coxph_c1$coefficients) == coxph_c1$linear.predictors


# event times are already ordered
event_times1 <- unique(xtemp1$time[xtemp1$status == 1])
event_times2 <- unique(xtemp1$time[xtemp1$status == 2])

event_elp1 <- unique(xtemp1$elp1[xtemp1$status == 1])
event_elp2 <- unique(xtemp1$elp2[xtemp1$status == 2])

# cause 1 baseline hazard
h0c1 <- vapply(event_times1, \(ti) {
  # subset of observations with time >= current time are at risk
  riskset <- xtemp1[xtemp1$time >= ti, ]
  # 1 / (exp(xT beta)) => h_0_i
  1/sum(riskset$elp1)
}, FUN.VALUE = numeric(1))

h0c2 <- vapply(event_times2, \(ti) {
  # subset of observations with time >= current time are at risk
  riskset <- xtemp1[xtemp1$time >= ti, ]
  # 1 / (exp(xT beta)) => h_0_i
  1/sum(riskset$elp2)
}, FUN.VALUE = numeric(1))

# cause 1 cumulative baseline hazard
H0c1 <- cumsum(h0c1)
H0c2 <- cumsum(h0c2)

# Sanity check should look more or less identical to cause1 cumulative baseline hazard from riskRegression
plot(x = event_times1[event_times1 %in% eval_times], y = H0c1[event_times1 %in% eval_times])
plot(x = cumhazards$time, y = cumhazards$cumhazard_rr)

S0c1 <- exp(-H0c1)
S0c2 <- exp(-H0c2)

# have H0 per cause and event times per cause, need to add up at each time point, but have different event times


# Overall survival would be
# exp(-H(t)), H(t) = H1(t) + H2(t)

# CIF would be
# cumsum(hc1(t) * S(t))


# fwelnet_mt / cooper fit -----------------------------------------------------------------------------------------

cooperfit <- fwelnet_mt_cox(
  data = sim1, causes = 1:2, mt_max_iter = 2, stratify_by_status = TRUE,
  alpha = 1, t = 10, thresh = 1e-4,
  verbose = TRUE, include_mt_beta_history = TRUE
)

cooperfit$fwfit1$glmfit$beta[, which(cooperfit$fwfit1$glmfit$lambda == cooperfit$fwfit1$lambda.min)]
cooperfit$beta1[, ncol(cooperfit$beta1)]

# getting linear predictor (type = "link" is default)
predict(cooperfit$fwfit1, xnew = sim1[,1:14], s = "lambda.min", type = "link")


# pammtools -------------------------------------------------------------------------------------------------------
library(pammtools)
library(mgcv)

sim1_ped <- sim1 %>%
  as_ped(Surv(time, status) ~ ., id = "id") |>
  mutate(cause = factor(cause))

pam_formula <- reformulate(c("s(tend, by = cause)", paste0("x", 1:13)), "ped_status")

pam_csh <- pamm(pam_formula, data = sim1_ped)
summary(pam_csh)
# pam_csh <- lapply(sim1_ped, \(x) pamm(pam_formula, data = x))
# lapply(pam_csh, summary)

sim1_cif_pam <- sim1_ped |>
  make_newdata(tend = unique(tend), cause = 1:2) |>
  group_by(cause) |>
  add_cif(pam_csh) |>
  mutate(cause = factor(cause))

ggplot(sim1_cif_pam, aes(x = tend, y = cif, color = cause, fill = cause)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), alpha = 1/5)


library(fwelnet)
library(ggplot2)

xp <- c("time", "status", "sex", "age", "trt", "logbili", "logprotime", "protimegrp", "stage")

train <- riskRegression::simPBC(500)[xp]
test <- riskRegression::simPBC(300)[xp]

fit <- fwelnet_mt_cox(
  train,
  standardize = TRUE, nfolds = 10, stratify_by_status = TRUE,
  thresh = 1e-7, t = 100,
  mt_max_iter = 10, 
  include_mt_beta_history = TRUE
)

fit$fwelfits[[1]]$glmfit$theta_store

betadf <- extract_beta_history(fit)

betadf[iteration <= 10, coef := round(coef, 3)] |>
  ggplot(aes(x = iteration, y = coef, color = x, fill = x)) +
  facet_grid(vars(x), vars(event), scales = "free_y") +
  geom_path() +
  theme_minimal()

betadf[, list(mean = mean(coef), sd = sd(coef)), by = "x"]

# Initial estimate (glmnet) and boxplots of subsequent betas by cooper
betadf |>
  ggplot(aes(x = x, y = coef, color = factor(event), fill = factor(event))) +
  facet_wrap(vars(event), ncol = 1, scales = "free_x") +
  geom_boxplot(alpha = 1/3) +
  geom_point(data = betadf[iteration == 1], size = 3, stroke = 0.3, color = "black", shape = 21, position = position_dodge(width = 2/3)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "top")

xp <- c("time", "status", "sex", "age", "trt", "logbili", "logprotime", "protimegrp", "stage")

train <- riskRegression::simPBC(300)[xp]

trainmat <- do.call(cbind, lapply(train, as.numeric))
trainmat <- trainmat[, -match(names(train), c("time", "status"), nomatch = 0)]

train$status[train$status == 2] <- 0

table(train$status)
z <- matrix(1, nrow = ncol(trainmat))

set.seed(13)
z <- matrix(rnorm(ncol(trainmat), sd = 4), ncol = 1)

set.seed(2)
cvg <- cv.glmnet(trainmat, y = survival::Surv(train$time, train$status), family = "cox", nfolds = 3, keep = TRUE)
set.seed(2)
cvf <- cv.fwelnet(trainmat, y = survival::Surv(train$time, train$status), z = z, family = "cox", lambda = cvg$lambda, foldid = cvg$foldid)

debugonce(glmnet:::cv.glmnet.raw)
debugonce(glmnet:::getOptcv.glmnet)

buildPredmat.fwelnetlist <- function()
  
View(glmnet::cv.glmnet)
View(glmnet:::cv.glmnet.raw)
View(glmnet:::buildPredmat.coxnetlist)
debugonce(glmnet:::buildPredmat.coxnetlist)


cvg$cvm
cvf$cvm

cvf$lambda
cvg$lambda

cvg$lambda.min
cvf$lambda.min

plot(cvg)
plot(cvf)

idx_g <- cvg$index[1]
idx_f <- which(cvf$lambda == cvf$lambda.min)

cbind(
  glmnet = cvg$glmnet.fit$beta[,idx_g],
  fwelnet = cvf$glmfit$beta[,idx_f],
  glmnet_in_fwelnet = cvf$glmfit$glmfit$beta[,idx_f]
)




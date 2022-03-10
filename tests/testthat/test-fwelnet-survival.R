library(survival)

x <- as.matrix(veteran[, 5:8]) # c(-3, -4)
y <- Surv(veteran$time, veteran$status)

test_that("Parity with glmnet", {
  # if z = 1, should be same as glmnet
  z <- matrix(1, ncol = 1, nrow = ncol(x))

  fwfit <- fwelnet(x, y, z, family = "cox")
  gfit <- glmnet::glmnet(x, y, family = "cox")

  expect_equal(fwfit$beta, as.matrix(gfit$beta), ignore_attr = TRUE, tolerance = 1e-5)
  expect_equal(gfit$lambda, fwfit$lambda)
})

lambda_seq <- seq(1e-4, .5, length.out = 20)

test_that("CV parity with glmnet", {
  # if z = 1, should be same as glmnet
  z <- matrix(1, ncol = 1, nrow = ncol(x))

  fwfit <- cv.fwelnet(x, y, z, family = "cox")
  gfit <- glmnet::cv.glmnet(x, y, family = "cox")
  expect_equal(gfit$lambda, fwfit$lambda)

  set.seed(1234)
  fwfit <- cv.fwelnet(x, y, z, family = "cox", lambda = lambda_seq)
  gfit <- glmnet::cv.glmnet(x, y, family = "cox", lambda = lambda_seq)
  expect_equal(gfit$lambda, fwfit$lambda)
  expect_equal(gfit$lambda.min, fwfit$lambda.min)

  set.seed(1234)
  fwfit <- cv.fwelnet(x, y, z, family = "cox", lambda = rev(lambda_seq))
  gfit <- glmnet::cv.glmnet(x, y, family = "cox", lambda = rev(lambda_seq))
  expect_equal(gfit$lambda, fwfit$lambda)
  expect_equal(gfit$lambda.min, fwfit$lambda.min)
  expect_equal(
    fwfit$glmfit$beta,
    as.matrix(gfit$glmnet.fit$beta),
    ignore_attr = TRUE, tolerance = 1e-5
  )
})

test_that("CV for cox is not much slower than glmnet", {
  z <- matrix(1, ncol = 1, nrow = ncol(x))
  t1 <- system.time(cv.fwelnet(x, y, z, family = "cox"))
  t2 <- system.time(glmnet::cv.glmnet(x, y, family = "cox"))

  expect_equal(mean(t1), mean(t2), tolerance = 0.1)
})



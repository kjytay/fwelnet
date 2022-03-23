n <- 100
p <- 20
x <- matrix(rnorm(n * p), n, p)
beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
y <- x %*% beta + rnorm(n)

test_that("Parity with glmnet", {
  # if z = 1, should be same as glmnet
  z <- matrix(1, ncol = 1, nrow = p)

  fwfit <- fwelnet(x, y, z, family = "gaussian")
  gfit <- glmnet::glmnet(x, y, family = "gaussian")

  expect_equal(fwfit$beta, as.matrix(gfit$beta), ignore_attr = TRUE, tolerance = 1e-5)
  expect_equal(gfit$lambda, fwfit$lambda)
})


test_that("Not-parity with glmnet", {
  # if z != 1, should be not-same as glmnet
  set.seed(5)
  z <- cbind(1, abs(beta) + rnorm(p, sd = 4), abs(beta) + rnorm(p, sd = 4))

  fwfit <- fwelnet(x, y, z, family = "gaussian")
  gfit <- glmnet::glmnet(x, y, family = "gaussian")

  expect_failure(
    expect_equal(
      fwfit$beta[1,], as.matrix(gfit$beta)[1,],
      ignore_attr = TRUE, tolerance = 1e-5
    )
  )
  expect_equal(gfit$lambda, fwfit$lambda)
})


# CV ----------------------------------------------------------------------

lambda_seq <- seq(1e-4, .5, length.out = 20)

test_that("CV parity with glmnet", {
  # if z = 1, should be same as glmnet
  z <- matrix(1, ncol = 1, nrow = p)

  # Even setting seeds would not ensure equality of lambda.min due to
  # differing numbers of rng calls resulting in slightly different folds?
  # Results only roughly equiv, not testable I think, requires foldid + nfolds
  fwfit <- cv.fwelnet(x, y, z, family = "gaussian")
  gfit <- glmnet::cv.glmnet(x, y, family = "gaussian")
  expect_equal(gfit$lambda, fwfit$lambda)

  fwfit <- cv.fwelnet(x, y, z, family = "gaussian", lambda = lambda_seq)
  gfit <- glmnet::cv.glmnet(x, y, family = "gaussian", lambda = lambda_seq)
  expect_equal(gfit$lambda, fwfit$lambda)


  fwfit <- cv.fwelnet(x, y, z, family = "gaussian", lambda = rev(lambda_seq))
  gfit <- glmnet::cv.glmnet(x, y, family = "gaussian", lambda = rev(lambda_seq))
  expect_equal(gfit$lambda, fwfit$lambda)
})

test_that("Complete parity with glmnet", {
  # Fully identical setup
  z <- matrix(1, ncol = 1, nrow = p)
  set.seed(4)
  nfolds <- 10
  foldid <- sample(rep(seq(nfolds), length = n))

  fwfit <- cv.fwelnet(x, y, z, family = "gaussian", lambda = lambda_seq, foldid = foldid, nfolds = nfolds)
  gfit <- glmnet::cv.glmnet(x, y, family = "gaussian", lambda = lambda_seq, foldid = foldid, nfolds = nfolds)

  expect_equal(gfit$lambda.min, fwfit$lambda.min)
  expect_equal(gfit$cvm, fwfit$cvm)
})

test_that("CV non-parity with glmnet", {
  # Fully non-identical setup
  set.seed(45)
  z <- cbind(1, abs(beta) + rnorm(p, sd = 6))
  nfolds <- 10
  foldid <- rep(seq(nfolds), length = n)

  fwfit <- cv.fwelnet(x, y, z, family = "gaussian", lambda = lambda_seq, foldid = foldid, nfolds = nfolds)
  gfit <- glmnet::cv.glmnet(x, y, family = "gaussian", lambda = lambda_seq, foldid = foldid, nfolds = nfolds)

  expect_equal(gfit$lambda, fwfit$lambda)
  # expect_failure(expect_equal(gfit$lambda.min, fwfit$lambda.min))
  expect_failure(expect_equal(gfit$cvm, fwfit$cvm))
})


# Fixing theta ------------------------------------------------------------

test_that("Fixing theta works", {
  z <- matrix(rnorm(p, sd = 2), ncol = 1, nrow = p)
  
  set.seed(3)
  fwfit <- fwelnet(x, y, z, family = "gaussian")

  fwfit_fix <- fwelnet(x, y, z, family = "gaussian", theta = fwfit$theta, lambda = fwfit$lambda)
  
  expect_identical(fwfit$theta, fwfit_fix$theta)
  expect_identical(fwfit$beta, fwfit_fix$beta)
  expect_identical(fwfit$nzero, fwfit_fix$nzero)
  expect_identical(fwfit$lambda, fwfit_fix$lambda)
})

test_that("Fixing theta works, reducing to glmnet solution", {
  z <- matrix(1, ncol = 1, nrow = p)
  
  set.seed(3)
  glmfit <- glmnet::glmnet(x, y, family = "gaussian")
  fwfit_fix <- fwelnet(x, y, z, family = "gaussian", theta = 1, lambda = glmfit$lambda)
  
  expect_identical(fwfit_fix$theta, as.matrix(1))
  expect_equal(fwfit_fix$beta, as.matrix(glmfit$beta), ignore_attr = TRUE, tolerance = 1e-5)
})

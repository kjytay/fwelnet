<!-- README.md is generated from README.Rmd. Please edit that file -->
Feature-Weighted Elastic Net (fwelnet)
======================================

`fwelnet` is a package that fits the ***feature-weighted elastic net
(fwelnet)***, a variant of the elastic net which has feature-specific
penalties. These penalties are based on additional information that the
user has on the features. `fwelnet` works with continuous and binary
responses. For details, please see the preprint (to come soon). For a
short tutorial on how to use the package, please see the vignette in the
`vignettes/` folder.

An example
----------

Here is a simple example to illustrate how to use this package. First,
let’s generate some data. In this example, we assume that we have 40
features, and that these features come in 4 groups of 10. The response
is a linear combination of the features from the first 2 groups with
additive Gaussian noise.

``` r
set.seed(1)
n <- 100
p <- 40
groups <- list(1:10, 11:20, 21:30, 31:40)  # which features belong to which group
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- matrix(rep(1:0, each = 20), ncol = 1)
y <- x %*% beta + rnorm(n)
```

In order to fit a fwelnet model, we have to define a feature information
matrix. In our example, we have **Z** ∈ ℝ<sup>40 × 4</sup>, with
*z*<sub>*j**k*</sub> = 1{feature *j* belongs to group *k*}.

``` r
# generate Z matrix
z <- matrix(0, nrow = p, ncol = length(groups))
for (i in 1:length(groups)) {
    z[groups[[i]], i] <- 1
}
```

Once `z` is specified, we can fit the fwelnet model with `fwelnet()`.

``` r
library(fwelnet)
fit <- fwelnet(x, y, z)
```

“fwelnet” objects are equipped with `predict` and `coef` methods which
allow the user to make predictions on new data and to view the model
coefficients. By default predictions and coefficients are returned for
the whole `lambda` path.

``` r
# predictions for first 5 observations at 20th lambda value
predict(fit, x[1:5, ])[, 20]
#  [1]  1.0342118  5.9991002 -1.8476885 -1.1713328 -0.5891373

# coefficients at the 20th lambda value (including intercept)
as.numeric(coef(fit)[, 20])
#   [1] -0.1724001  0.5835508  0.4862924  0.4525343  0.7083713  0.9203059
#   [7]  0.6250100  0.7257856  1.1662045  0.9173184  0.6505217  0.6570609
#  [13]  0.6398280  0.9016092  0.6402918  0.7126133  0.8255739  1.2130462
#  [19]  0.4603750  1.0811643  0.6380584  0.0000000  0.0000000  0.0000000
#  [25]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#  [31]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#  [37]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
```

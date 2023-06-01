#' Stratified CV folds for `cv.fwelnet`
#'
#' @param xdf An input dataset with at least `"status"` variable.
#' @param nfolds (`10`) Integer number of CV folds, default is equal to that of [`cv.glmnet`].
#'
#' @return A `data.frame` with variables `folds` for fold assignments and `status` for checking proportions.
#' @export
#' @importFrom data.table as.data.table :=
#'
#' @examples
#' xdf <- data.frame(
#'   row_id = 1:50, 
#'   status = c(rep(0, 20), rep(1, 18), rep(2, 12)) 
#' )
#' 
#' resdf <- stratified_cv_folds(xdf)
#' table(resdf$status, resdf$fold)
stratified_cv_folds <- function(xdf, nfolds = 10L) {
  if (!checkmate::test_data_table(xdf)) xdf <- data.table::as.data.table(xdf)
  checkmate::assert_data_table(xdf, any.missing = FALSE)
  checkmate::assert("status" %in% names(xdf))
  checkmate::assert_int(nfolds, lower = 2)
  cleanup <- TRUE
  if (min(floor(table(xdf[["status"]]) / nfolds)) == 0) {
    stop("nfolds likely too large, expecting empty cv folds")
  }
  
  if (!("row_id" %in% names(xdf))) {
    xdf[, row_id := .I]
  }
  
  # Stratified sampling with DT and helpers ripped straight from mlr3
  shuffle <- function(x, n = length(x), ...) {
    x[sample.int(length(x), size = n, ...)]
  }
  
  seq_along0 <- function(x) {
    n <- length(x)
    if (n >= 1L)  {
      0L:(n - 1L)
    } else {
      integer(0L)
    }
  }
  
  xdf[, fold := shuffle(seq_along0(row_id) %% as.integer(nfolds) + 1L), by = "status"]
  as.data.frame(xdf[, c("row_id", "fold", "status")])
}

globalVariables(c(".I", "row_id", "fold"))

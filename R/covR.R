#' @title A covariance calculator using R 
#' @description A function for computing the covaiance of two set of samples using R
#' @param X the first set of samples
#' @param Y the second set of samples with the same size of X
#' @return the covariance of the data
#' @examples
#' \dontrun{
#' X <- runif(100)
#' Y <- rnorm(100)
#' cov <- covR(X,Y)
#' }
#' @export
covR <- function(X, Y) {
  n <- length(X)
  cov <- sum((X-mean(X))*(Y-mean(Y)))/(n-1)
  return(cov)
}
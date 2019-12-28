#' @title A correlation coefficient calculator using R
#' @description A function for computing the correlation coefficient of two set of samples using R
#' @param X the first set of samples
#' @param Y the second set of samples with the same size of X
#' @return the correlation coefficient of the data
#' @examples
#' \dontrun{
#' X <- runif(100)
#' Y <- rnorm(100)
#' corr <- corrR(X,Y)
#' }
#' @export
corrR <- function(X, Y) {
  n <- length(X)
  varX <- sum((X-mean(X))^2)/(n-1)
  varY <- sum((Y-mean(Y))^2)/(n-1)
  cov <- sum((X-mean(X))*(Y-mean(Y)))/(n-1)
  corr <- cov/sqrt(varX*varY)
  return(corr)
}
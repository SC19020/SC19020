#' @title A kernel density estimator using R
#' @description A function for estimating the density function of a given set of data using R
#' @param X the given set of data of one-dimension
#' @param h the bandwidth
#' @param x the estimated density function value where you want to know
#' @importFrom stats dnorm
#' @return the estimated density function value
#' @examples
#' \dontrun{
#' X <- rnorm(100)
#' x <- seq(-1, 1, by = 0.01)
#' h <- 0.1
#' for(i in 1:length(x)){
#'   y[i] <- KDER(X, h, x[i])
#'   }
#' plot(x, y, type="l")
#' }
#' @export
KDER <- function(X, h, x) {
  n <- length(X)
  out <- (1/h)*mean(dnorm((x-X)/h))
  return(out)
}
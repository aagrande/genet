#' @title The Power Law Distribution
#'
#' @description These functions provide information about a power law
#' distribution on \code{seq(x_min, x_max)} with exponent \code{alpha}.
#' \code{dpowerlaw} gives the pmf, \code{qpowerlaw} gives the
#' quantile function, and \code{rpowerlaw} generates random deviates.
#'
#' @name power_law
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param alpha power law's exponent. Must be positive.
#' @param x_min,x_max lower and upper bounds of the distribution.
#' Must be finite.
#' @param step number; size of increment of the sequence defining the support.
#'
#' @details
#' The pmf of the distribution reads:
#' \deqn{p(x;\alpha) = x^{-\alpha} / C,}
#' for \eqn{x = x_min, x_min + 1, \dots, x_max}, \eqn{\alpha>0}, and where
#' \eqn{C>0} is a normalizing constant.
#'
#' @return \code{dpowerlaw} gives the pmf, \code{qpowerlaw} gives the
#' quantile function, and \code{rpowerlaw} generates random deviates.
#'
#' @export
#'
dpower_law <- function(x, alpha, x_min, x_max, step = 1) {
  if (x_min >= x_max) stop("`x_min` should be < than `x_max`.")
  if (any(c(is.infinite(x_min), is.infinite(x_max)))) {
    stop("`x_min` and `x_max` must be finite.")
  }
  x_min <- x_min
  x_max <- x_max
  range <- seq(x_min, x_max, step)
  in_range <- which(x >= x_min & x <= x_max)
  log_c <- -log(sum(range ^ (- alpha)))
  p <- rep(0, length(x))
  p[in_range] = exp(log_c - alpha * log(x[in_range]))
  return(p)
}

#' @rdname power_law
#' @export
rpower_law <- function(n, alpha, x_min, x_max, step = 1) {
  if (x_min >= x_max) stop("`x_min` should be < than `x_max`.")
  if (any(c(is.infinite(x_min), is.infinite(x_max)))) {
    stop("`x_min` and `x_max` must be finite.")
  }
  range <- seq(x_min, x_max, step)
  x <- sample(x = range, size = n, replace = T,
              prob = dpower_law(range, alpha, x_min, x_max, step))
  return(x)
}

#' @rdname power_law
#' @export
qpower_law <- function(p, alpha, x_min, x_max, step = 1){
  if (any(x_min >= x_max)) stop("`x_min` should be < than `x_max`.")
  if (any(c(is.infinite(x_min), is.infinite(x_max)))) {
    stop("`x_min` and `x_max` must be finite.")
  }
  if (any(c(p < 0, p > 1))) stop("p must be in [0, 1].")
  range <- seq(x_min, x_max, step)
  cdf <- cumsum(dpower_law(range, alpha, x_min, x_max, step))
  idx <- sapply(p, function(x) min(which(x < cdf)))
  x <- range[idx[!is.infinite(idx)]]
  x <- c(x, rep(x_max, sum(is.infinite(idx))))
  if (length(x) > 1) {
    for (i in 2:length(x)){
      if (x[i] <= x[i-1]) x[i] = x[i-1] + step
    }
  }
  return(x)
}

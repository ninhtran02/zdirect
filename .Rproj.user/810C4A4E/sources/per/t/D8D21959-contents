#' Noncentral t-statistics normalization
#'
#' An approximate normalization method for noncentral t-distributed z-values. The resulting normalized z-values will be approximately normally distributed with standard deviation 1.
#'
#' @param zv A numeric vector of noncentral t-statistics.
#' @param v A numeric vector of the degrees of freedom.
#'
#' @return Normalized noncentral t-statistics.
#' @export
#'
#' @examples
noncentral_t_normalization <- function(zv,v){
  a <- sqrt(v/(v-2))
  b <- sqrt(2 * gamma(0.5*v)^2 / ((v-2)*gamma(0.5*v - 0.5)^2)  - 1)
  ALPHA <- 1/b
  BETA <- b/a
  return(ALPHA * asinh(BETA * zv))
}

dStorey <- function(p,alpha, B = 1000, lambdas = seq(0.01,0.99,0.01)){

  m <- length(p)

  pi_vec <- c()
  for(lambda in lambdas){
    sorted.p <- sort(p)
    W_lambda <- sum(p > lambda)
    pi_lambda <- (W_lambda + 1)/( m * (1 - lambda) )
    pi_vec <- c(pi_vec, pi_lambda)
  }
  min_pi <- min(pi_vec)

  MSE_vec <- c()
  for(lambda in lambdas){
    pi.b.vec <- c()
    for(b in 1:B){
      p.b <- sample(x = p, size = m, replace = TRUE)
      sorted.p.b <- sort(p.b)
      W_lambda <- sum(p.b > lambda)
      pi_lambda <- (W_lambda + 1)/( m * (1 - lambda) )
      pi.b.vec <- c(pi.b.vec,pi_lambda)
    }
    MSE <- (1/B)*sum(pi.b.vec - min_pi)^2
    MSE_vec <- c(MSE_vec, MSE)
  }
  chosen.lambda <- lambdas[which.min(MSE_vec)[1]]

  sorted.p <- sort(p)
  W_lambda <- sum(p > chosen.lambda)
  pi_lambda <- (W_lambda + 1)/( m * (1 - chosen.lambda) )
  dFDR <- (pi_lambda*sorted.p)/( (1:m) / m )
  dFDR[sorted.p > chosen.lambda] <- 1
  ref <- sum(dFDR <= alpha)

  return( ref )
}

#' Directional Storey, Taylor and Siegmund procedure.
#'
#' An adaptive procedure that controls the directional false discovery rate by estimating the null proportion \eqn{\pi_0}.
#'
#' @param pv A numeric of p-values.
#' @param alpha The directional false discovery rate target.
#' @param B Number of bootstrap samples.
#' @param lambdas A numeric of candidate \eqn{\lambda} values for the automatic directional Storey, Taylor and Siegmund procedure. If a single value is provided, then the standard directional Storey, Taylor and Siegmund procedure will run.
#'
#' @return The rejection indices.
#' @export
#'
#' @examples
STS <- function(pv, alpha, B = 1000, lambdas = 0.5){
  # Storey rejection set
  dStorey.obj <- dStorey(pv,
                         alpha,
                         B = B,
                         lambdas = lambdas)
  if(dStorey.obj == 0){
    rej_set_Storey = integer(0)
  } else {
    rej_set_Storey = which(pv <= sort(pv)[dStorey.obj])
  }

  return(rej_set_Storey)
}

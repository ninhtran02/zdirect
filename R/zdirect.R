#' zdirect
#'
#' A procedure that controls the directional false discovery rate via "data masking" for z-values of the form \eqn{z_i \sim \text{N}(\theta,\sigma_i)}. If the z-values available are noncentral t-distributed, \eqn{\texttt{noncentral\_t\_normalization}} should be used to approximately normalize them first before inputting them into \eqn{\texttt{zdirect}}.
#'
#' @param betahat A numeric of z-values \eqn{(z_1,\dots,z_m)}.
#' @param sebetahat A numeric of the standard deviations \eqn{(\sigma_1,\dots,\sigma_m)}.
#' @param s_l The left bound for \eqn{\mathcal{M}_1}. See below.
#' @param s_r The right bound for \eqn{\mathcal{M}_1} so that \eqn{ \mathcal{M}_1 = \{ i: u'_i \leq \texttt{s\_l} \text{ or } u'_i \geq \texttt{s\_r} \} }.
#' @param alpha The directional false discovery rate target.
#' @param mixcompdist The approximation model \eqn{\hat{G}(\cdot)} for \eqn{G(\cdot)}. If \eqn{\texttt{"halfuniform"}} is chosen, then \eqn{\hat{G}(\cdot)} will be a mixture of uniform distributions. If \eqn{\texttt{"halfnormal"}} is chosen, then \eqn{\hat{G}(\cdot)} will be a mixture of mixture of truncated normal distributions.
#' @param prior A string, or numeric vector indicating Dirichlet prior on mixture proportions: \eqn{\texttt{"nullbiased"}}, \eqn{\texttt{c(nullweight,altweight,...,altweight)}}, puts more weight on first component; \eqn{\texttt{"uniform"}} is \eqn{\texttt{c(1,1...,1)}}; \eqn{\texttt{"unit"}} is \eqn{\texttt{c(1/K,...,1/K)}}.
#' @param nullweight A scalar, the weight put on the prior for nulls under "nullbiased" specification, see \eqn{\texttt{prior}}.
#' @param altweight A scalar, the weight put on the prior for non-nulls under "nullbiased" specification, see \eqn{\texttt{prior}}.
#' @param nfits A scalar that ensures that \eqn{\hat{g}_t(\cdot)} is only re-estimated every \eqn{\lceil m / \texttt{nfits} \rceil} steps where \eqn{m} is the length of \eqn{\texttt{betahat}}.
#' @param epsilon A small scalar which widens the width of the non-masked data set at each iteration.
#'
#' @return The rejection indices.
#' @export
#'
#' @examples
zdirect = function(betahat, sebetahat,
                 s_l = NULL,
                 s_r =  NULL,
                 alpha = 0.1,
                 mixcompdist = "halfuniform",
                 prior = "nullbiased",
                 nullweight = 0.8,
                 altweight = 1e+08,
                 nfits = 20,
                 epsilon = 1e-10
                 ){
  # initialize the true beta in the same way as REBayes::KWDual() function
  # eps <- 1e-4
  # beta <- seq(min(betahat) - eps, max(betahat) + eps, length = K)


  n = length(betahat)
  if (is.null(s_l)) s_l = rep(0.2, n)
  if (is.null(s_r)) s_r = rep(0.8, n)


  parts_init = initparts(betahat, sebetahat) # get the initialized parts for the zdirect procedure
  r_set = parts_init$r_set
  l_set = parts_init$l_set
  betahat_mirror = parts_init$betahat_mirror
  betahat_outer = parts_init$betahat_outer
  U = parts_init$U
  U_outer = parts_init$U_outer
  unmask_set= which((s_l <= U & U<=  (0.5 - s_l))| ((1.5 -   s_r) <= U & U <= s_r ))
  mask_set = setdiff(1:n, unmask_set)
  rej_set = which( (U <  s_l |   U > s_r) )
  acp_set <- which(U>  (0.5 - s_l) &   U  < (1.5 -   s_r))
  R = length(rej_set)
  A = length(acp_set)
  FDPest = (1 + A)/max(R, 1)


  count = 0
  while(FDPest > alpha & R > 0 ) {
    if (count%%(n%/%nfits) == 0){
    # print(paste("THIS IS THE ", count%/%(n%/%nfits) +1, " TIME OF IP UPDATE. ","FDP_est = ", FDPest ))
      get_lfsr_mask_obj = get_lfsr_mask(
                                betahat,
                                betahat_mirror,
                                betahat_outer,
                                sebetahat,
                                mask_set,
                                mixcompdist,
                                prior,
                                nullweight,
                                altweight)

      lfsr_mask = get_lfsr_mask_obj$lfsr_mask
      matrix_lik_outer_mask = get_lfsr_mask_obj$matrix_lik_outer_mask

      # print(paste0("length of mask_set is ", length(mask_set) ))

    # print(paste0("dim of matrix_lik_outer_mask is ", dim(matrix_lik_outer_mask) ))

    }

    max_position = which.max(lfsr_mask)
    new_s = U_outer[mask_set[max_position]]

    # print(paste0("max position is ", max_position))
    # print(paste0("new_s is ", new_s))
    if (new_s <= 0.5){
      s_l[mask_set[max_position]] = new_s - epsilon
    }else{
      s_r[mask_set[max_position]] = new_s + epsilon
    }


    #
    unmask_set = which((s_l <= U & U<=  (0.5 - s_l))| ((1.5 -   s_r) <= U & U <= s_r ))
    mask_set = setdiff(1:n, unmask_set)
    lfsr_mask = lfsr_mask[- max_position]
    # print(paste0("length of lfsr_mask is ", length(lfsr_mask) ))
    # print(paste0("length of mask_set is ", length(mask_set) ))

    ## update FDP estimate, rej_set, and count variable
    rej_set <- which( (U <  s_l |   U > s_r) )
    acp_set <- which(U>  (0.5 - s_l) &   U  < (1.5 -   s_r))
    R = length(rej_set)
    A = length( acp_set )
    FDPest = (1 + A)/max(R, 1)
   # print(paste0("min of s_l is ", min(s_l) ))
   # print(paste0("max of s_r is ", max(s_r)))
    count = count +1
  }
# print("zdirect finite done!")
  return(list(rej_set = rej_set,
              acp_set = acp_set
              # ,
              # pihat = get_lfsr_mask_obj$pihat,
              # a = get_lfsr_mask_obj$a,
              # b = get_lfsr_mask_obj$b
              )
         )

}



## intialize parts for the zdirect procedure
initparts  = function(betahat, sebetahat){
  n = length(betahat)
  r_set = which(betahat > 0)
  l_set = setdiff( 1:n, r_set)
  betahat_mirror = numeric(n)
  betahat_outer = numeric(n)
  betahat_mirror[r_set] = qnorm(0.5 + pnorm(q = betahat[r_set],
                                            lower.tail = F, sd = sebetahat[r_set]), sd = sebetahat[r_set])
  betahat_mirror[l_set] = qnorm(0.5 - pnorm(q = betahat[l_set],
                                            lower.tail = T, sd = sebetahat[l_set]), sd = sebetahat[l_set])
  betahat_outer [r_set] = pmax( betahat[r_set] , betahat_mirror[r_set] )
  betahat_outer [l_set] = pmin( betahat[l_set] , betahat_mirror[l_set] )
  U = pnorm(q = betahat,sd = sebetahat)
  U_outer = pnorm(q =betahat_outer, sd = sebetahat)

  return(
    list(
      r_set = r_set,
      l_set = l_set,
      betahat_mirror = betahat_mirror,
      betahat_outer = betahat_outer,
      U = U,
      U_outer = U_outer
    )
  )

}





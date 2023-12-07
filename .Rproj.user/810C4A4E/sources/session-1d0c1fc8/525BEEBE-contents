

# compute lfsr quantities based on partially masked data
get_lfsr_mask = function(betahat,
                         betahat_mirror, 
                         betahat_outer, 
                         sebetahat, 
                         mask_set, 
                         mixcompdist = c( "halfuniform","halfnormal"),
                         prior =  prior, 
                         nullweight = 10, 
                         altweight = 1,
                         mode = 0,
                         mixsd = NULL,
                         gridmult = sqrt(2),
                         grange = c(-Inf, Inf),
                         lik = NULL,
                         df = NULL,
                         pointmass = TRUE
                         ){
  n = length(betahat) 
  unmask_set = setdiff(1:n, mask_set)
  
  # determine the likelihood used
  if(is.null(lik)){
    if(is.null(df)){
      lik = lik_normal()
    } else {lik = lik_t(df)}
  }
  
  # determine mixsd
  if (is.null(mixsd)) {
    betahat_data = betahat_outer
    betahat_data[unmask_set] = betahat[unmask_set]
    mixsd = 
      autoselect.mixsd( ashr::set_data(betahat_data, sebetahat, lik), 
      gridmult, mode, grange, mixcompdist)
  }
  if (pointmass) {
    mixsd = c(0, mixsd)
  }
  k  = length(mixsd)
  prior = setprior(prior,k,nullweight, altweight, which.min(mixsd)) # parameters of Dirichlet prior on mixing proportions. First coordinate corresponds to null if pointmass is true

  ## construct the component likelihood matrix, for g1 components
  if (mixcompdist == "halfuniform"){
    if (min(mixsd) > 0) { #meaning no point mass
      # Simply reflect the components.
      a = c(mode - mixsd, rep(mode, k))  # the left endpts of the half uniform g1 components
      b = c(rep(mode, k), mode + mixsd)  # the right endpts of the half uniform g1 components
      prior = c(prior, prior)
    } else {
      # Define two sets of components, but don't duplicate null component. 
      null.comp = which.min(mixsd)  # should be 1
      a = c(mode - mixsd, rep(mode, k-1)) # if mode = 0, then null followed by -ve of the mixsd followed by zeroes
      b = c(rep(mode, k), (mode + mixsd)[-null.comp]) # if mode = 0, then null followed by zeroes followed by mixsd
      prior = c(prior, prior[-null.comp])
    }
    matrix_lik = matrix_lik_halfuniform(betahat, betahat_mirror,
                                        sebetahat, a, b,  unmask_set, mask_set, pointmass)
  } else if (mixcompdist == "halfnormal"){
    if (min(mixsd) > 0) { #no point mass
      mixsd = c(mixsd, mixsd)
      prior = c(prior, prior)
    } else {
      null.comp = which.min(mixsd)
      mixsd = c(mixsd, mixsd[-null.comp])
      prior = c(prior, prior[-null.comp])
    }
    
    matrix_lik = matrix_lik_halfnormal(betahat, betahat_mirror, sebetahat, 
                                       mixsd, unmask_set, mask_set, pointmass)
  }
    
    pihat = ashr::mixIP(matrix_lik= matrix_lik, prior = prior)$pihat 

  # print(paste0("est. prob of point mass " , pihat[1]))
  # print(paste0("est. prob of left mass " ,sum(pihat[2:k])))
  # print(paste0("est. prob of right mass ", sum(pihat[(k+1):(2*k - 1)])))
    
    if (mixcompdist == "halfuniform"){
      matrix_lik_outer_mask = matrix_lik_halfuniform_outer(betahat_outer[mask_set], sebetahat[mask_set], a, b,  pointmass)
    } else if (mixcompdist == "halfnormal"){
      matrix_lik_outer_mask = matrix_lik_halfnormal_outer(betahat_outer[mask_set], sebetahat[mask_set], mixsd,  pointmass)
    }
    
    lik_outer_mask = as.numeric(matrix_lik_outer_mask%*%pihat)  # this step should only use betahat_outer
    
    if (min(mixsd) > 0) {  # no null point
      
      neg_lik_outer_mask = as.numeric(matrix_lik_outer_mask[, 1:k]%*%pihat[1:k])
      
      lfsr_mask = ashr::compute_lfsr(NegativeProb = neg_lik_outer_mask/lik_outer_mask,
                                     ZeroProb= rep(0, length(mask_set))
                                     )
    }else{
    
      neg_lik_outer_mask = as.numeric(matrix_lik_outer_mask[, 2:k]%*%pihat[2:k])
      zero_lik_outer_mask = as.numeric(matrix_lik_outer_mask[ , 1]*pihat[1])

      lfsr_mask = ashr::compute_lfsr(NegativeProb = neg_lik_outer_mask/lik_outer_mask,
                                     ZeroProb = zero_lik_outer_mask/lik_outer_mask)
    }
  
  

  
  return(list( lfsr_mask = lfsr_mask,
               matrix_lik_outer_mask = matrix_lik_outer_mask,
               pihat = pihat))
}

# computing the L matrix for null and half-uniform components
matrix_lik_halfuniform = function(betahat, betahat_mirror, sebetahat, 
                                  a, b, unmask_set, mask_set, pointmass){
  K = length(a)
  n = length(betahat)
  n_mask = length(mask_set)
  n_unmask = length(unmask_set)  
  matrix_lik = matrix(0, nr = n, nc = K)
  
  unmask_entries = 
    pnorm(q = rep(betahat[unmask_set], each= K),
          mean  = rep(a, n_unmask), 
          sd = rep(sebetahat[unmask_set], each = K))-
    pnorm(q = rep(betahat[unmask_set], each= K) ,
          mean  = rep(b, n_unmask), 
          sd = rep(sebetahat[unmask_set], each = K))
  
  mask_entries = 
    (pnorm(q = rep(betahat[mask_set], each= K) ,
          mean  = rep(a, n_mask), 
          sd = rep(sebetahat[mask_set], each = K)) +
    pnorm(q = rep(betahat_mirror[mask_set], each= K) ,
          mean  = rep(a, n_mask),
          sd = rep(sebetahat[mask_set], each = K))
    )-
    (pnorm(q = rep(betahat[mask_set], each= K) ,
          mean  = rep(b, n_mask), 
          sd = rep(sebetahat[mask_set], each = K)) +
    pnorm(q = rep(betahat_mirror[mask_set], each= K) ,
          mean  = rep(b, n_mask),
          sd = rep(sebetahat[mask_set], each = K))

    )
  
  unmask_entries = unmask_entries/(rep(b, n_unmask) - rep(a, n_unmask))
  mask_entries = mask_entries/(rep(b, n_mask) - rep(a, n_mask))
  
  matrix_lik[unmask_set, ] = matrix(unmask_entries, byrow = T, nc = K)
  matrix_lik[mask_set, ] = matrix(mask_entries, byrow = T, nc = K)
  
  if (pointmass) {
    matrix_lik[unmask_set, 1] = 
      dnorm(x = betahat[unmask_set], mean = 0, sd = sebetahat[unmask_set])
    
    matrix_lik[mask_set, 1] = 
      dnorm(x = betahat[mask_set], mean = 0, sd = sebetahat[mask_set])+
      dnorm(x = betahat_mirror[mask_set], mean = 0, sd = sebetahat[mask_set])
  }
  return(matrix_lik)

}
  


# this should return a length(betahat_outer) by length(a) class marginal likelihood matrix, for the "outer" betahat
matrix_lik_halfuniform_outer = function(betahat_outer, sebetahat, a, b,  pointmass){
  K = length(a)
  n = length(betahat_outer)
  lik_entries = 
    pnorm(q = rep(betahat_outer, each= K),
          mean  = rep(a, n), 
          sd = rep(sebetahat, each = K))  - 
    pnorm(q = rep(betahat_outer, each= K) ,
          mean  = rep(b, n), 
          sd = rep(sebetahat, each = K))
  
  lik_entries = lik_entries/(rep(b, n) - rep(a, n))
  matrix_lik = matrix(lik_entries, byrow = T, nc = K)
  if (pointmass)  matrix_lik[, 1] = dnorm(x = betahat_outer, mean = 0, sd = sebetahat)
  return(matrix_lik)
}

# computing the L matrix for null and half-normal components
matrix_lik_halfnormal = function(betahat, betahat_mirror, sebetahat, 
                                 mixsd, unmask_set, mask_set, pointmass){
  K = length(mixsd)
  H = K%/%2
  n = length(betahat)
  n_mask = length(mask_set)
  n_unmask = length(unmask_set)  
  matrix_lik = matrix(0, nr = n, nc = K)  # the matrix will be n by K
  
  if (pointmass == T){
    left_trunc_pts = c(-Inf, rep(-Inf,H), rep(0, H) )
    right_trunc_pts = c(Inf, rep(0,H), rep(Inf, H))
  } else {
    left_trunc_pts = c( rep(-Inf,H), rep(0, H) )
    right_trunc_pts = c(rep(0,H), rep(Inf, H))
  }
  
  
  unmask_entries = tnn_conv(v = rep(betahat[unmask_set], each= K), 
                            s = rep(sebetahat[unmask_set], each = K), 
                            mu = 0, 
                            sigma = rep(mixsd, n_unmask),
                            a = rep(left_trunc_pts, n_unmask),
                            b = rep(right_trunc_pts, n_unmask))
  
  mask_entries = 
    tnn_conv(v = rep(betahat[mask_set], each= K), 
             s = rep(sebetahat[mask_set], each = K), 
             mu = 0, 
             sigma = rep(mixsd, n_mask),
             a = rep(left_trunc_pts, n_mask),
             b = rep(right_trunc_pts, n_mask)) +
    tnn_conv(v = rep(betahat_mirror[mask_set], each= K),
             s = rep(sebetahat[mask_set], each = K),
             mu = 0,
             sigma = rep(mixsd, n_mask),
             a = rep(left_trunc_pts, n_mask),
             b = rep(right_trunc_pts, n_mask))

  matrix_lik[unmask_set, ] = matrix(unmask_entries, byrow = T, nc = K)
  matrix_lik[mask_set, ] = matrix(mask_entries, byrow = T, nc = K)
  
  if (pointmass) {
    matrix_lik[unmask_set, 1] = 
      dnorm(x = betahat[unmask_set], mean = 0, sd = sebetahat[unmask_set])
    
    matrix_lik[mask_set, 1] = 
      dnorm(x = betahat[mask_set], mean = 0, sd = sebetahat[mask_set]) +
      dnorm(x = betahat_mirror[mask_set], mean = 0, sd = sebetahat[mask_set])
  }
  return(matrix_lik)
  
}




matrix_lik_halfnormal_outer = function(betahat_outer, sebetahat, mixsd,  pointmass){
  K = length(mixsd)
  H = K%/%2
  n = length(betahat_outer)
  
  if (pointmass == T){
    left_trunc_pts = c(-Inf, rep(-Inf,H), rep(0, H) )
    right_trunc_pts = c(Inf, rep(0,H), rep(Inf, H))
  } else {
    left_trunc_pts = c( rep(-Inf,H), rep(0, H) )
    right_trunc_pts = c(rep(0,H), rep(Inf, H))
  }
  
  lik_entries = tnn_conv(v = rep(betahat_outer, each= K), 
                            s = rep(sebetahat, each = K), 
                            mu = 0, 
                            sigma = rep(mixsd, n),
                            a = rep(left_trunc_pts, n),
                            b = rep(right_trunc_pts, n))
  matrix_lik = matrix(lik_entries, byrow = T, nc = K)
  
  if (pointmass) {
    matrix_lik[, 1] = dnorm(x = betahat_outer, mean = 0, sd = sebetahat)
  }
  return(matrix_lik)
}



get_exclusions=function(data){
  return((data$s==0 | data$s == Inf | is.na(data$x) | is.na(data$s)))
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mode is the location about which inference is going to be centered
# mult is the multiplier by which the sds differ across the grid
# grange is the user-specified range of mixsd
autoselect.mixsd = function(data,mult,mode,grange,mixcompdist){
  if (data$lik$name %in% c("pois")){
    data$x = data$lik$data$y/data$lik$data$scale #estimate of lambda
    data$s = sqrt(data$x)/data$lik$data$scale #standard error of estimate
    # if the link is log we probably want to take the log of this?
  }
  if (data$lik$name %in% c("binom")){
    data$x = data$lik$data$y
  }
  
  betahat = data$x - mode
  sebetahat = data$s
  exclude = get_exclusions(data)
  betahat = betahat[!exclude]
  sebetahat = sebetahat[!exclude]
  
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }
  
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}


setprior=function(prior,k,nullweight, altweight, null.comp){
  if(!is.numeric(prior)){
    if(prior=="nullbiased"){ # set up prior to favour "null"
      prior = rep(altweight,k)  # (rep(1, k) is the standard prior values for non-nulls by M. Stephens)
      prior[null.comp] = nullweight 
    }else if(prior=="uniform"){
      prior = rep(1,k)
    } else if(prior=="unit"){
      prior = rep(1/k,k)
    }
  }
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  return(prior)
}

### not used at the moment, but intend to write a function that implements point-mass prior like KW did
# Lmat.kw = function(betahat, betehat_mirror, beta, mask_set, sebetahat){
#    
#    K <- length(beta)
#    m <- length(betahat) 
#    L <- matrix(0, nr = m, nc = K)
#    unmask_set <-  setdiff(1:m, mask_set)
# 
#    # compute rows in L for the unmasked betahat's, row by row
#     L[unmask_set, ] <- 
#       matrix(
#       dnorm(x = rep(betahat[unmask_set], each= K), 
#             mean= rep(beta, length(unmask_set)), 
#             sd = rep(sebetahat[unmask_set], each = K)),
#       byrow = T, 
#       nc = K)
#     
#     # compute rows in L for the masked betahat's
#     L[mask_set, ] <-
#       matrix(
#         dnorm(x = rep(betahat[mask_set], each= K),
#               mean = rep(beta, length(mask_set)), 
#               sd = rep(sebetahat[mask_set], each = K)) + 
#         dnorm(x = rep(betahat_mirror[mask_set], each= K),
#                 mean =  rep(beta, length(mask_set)), 
#                 sd = rep(sebetahat[mask_set], each = K)),
#         byrow = T, 
#         nc = K)
#     
#    return(L)
#   
# }





##### this is just a temporary function to give final estimate fo
# estimate_lfsr = function(betahat, betahat_mirror, sebetahat, 
#                          pihat, a, b, pointmass = T){
#   K = length(a)
#   k = (K-1)/2
#   n = length(betahat)
#   matrix_lik = matrix(0, nr = n, nc = K)
#   
#   entries = 
#     pnorm(q = rep(betahat, each= K),
#           mean  = rep(a, n), 
#           sd = rep(sebetahat, each = K))-
#     pnorm(q = rep(betahat, each= K) ,
#           mean  = rep(b, n), 
#           sd = rep(sebetahat, each = K))
#   
#   entries = entries/(rep(b, n) - rep(a, n))
#   matrix_lik = matrix(entries, byrow = T, nc = K)
#   
#   if (pointmass) {
#     matrix_lik[, 1] = dnorm(x = betahat, mean = 0, sd = sebetahat)
#     
#   }
#   lik = as.numeric(matrix_lik%*%pihat)  
#   neg_lik = as.numeric(matrix_lik[, 2:k]%*%pihat[2:k])
#   zero_lik = as.numeric(matrix_lik[ , 1]*pihat[1])
#   
#   lfsr = ashr::compute_lfsr(NegativeProb = neg_lik/lik,
#                                  ZeroProb = zero_lik/lik)
#   
#   
#   return(lfsr)
#   
# }


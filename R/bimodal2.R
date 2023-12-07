


# generate data with bimodal alternatives
# it will also compute lfsr based on ash, the oracle lfsr and oracle lfdr


bimodal= function(pi0 = 0.7,  sebetahat = rep(1, 5000),
                  g_mean_choice = c(-2, 2),
                  g_sd_choice = c(1, 1),
                  g_prob_choice = c(0.5, 0.5),
                  alpha = 0.1, nullweight = 10, rho = 0
){
  m = length(sebetahat)
  mean_pos = sample(x = 1:length(g_prob_choice),
                    size = m,
                    replace = TRUE,
                    prob = g_prob_choice)

  # generate beta from g1 and adding  zeros for nulls
  beta = rnorm(m, g_mean_choice[mean_pos], g_sd_choice[mean_pos])
  beta[sample(x = c(T, F) , size =  m,  replace = TRUE, prob =  c(pi0, 1 - pi0))] = 0
  # generate betahat
  betahat <- beta + as.numeric(arima.sim(list(order=c(1,0,0), ar=rho), n=m))
  #betahat = rnorm(m, beta, sebetahat)
  # get Stephen's lfsr
  ash_obj = ash(betahat= betahat, sebetahat=  sebetahat, mixcompdist= "halfuniform", nullweight = nullweight)
  lfsr_ash = get_lfsr(ash_obj)
  ash_direction = rep(-1, m)
  ash_direction[which(get_pp(ash_obj) >get_np(ash_obj)) ] = 1
  # oracles (both lfdr and lfsr)
  lfdr_numerator = pi0*dnorm(x= betahat,  sd = sebetahat)
  lfsr_numerator_left = pi0*dnorm(x= betahat,  sd = sebetahat)
  lfsr_numerator_right = pi0*dnorm(x= betahat,  sd = sebetahat)
  denominator = pi0*dnorm(x= betahat,  sd = sebetahat)

  # common denominator
  for (j in 1:length(g_mean_choice)){
    denominator = denominator + (1 - pi0)*g_prob_choice[j]*dnorm(x = betahat,
                                                                 mean = g_mean_choice[j],
                                                                 sd= sqrt(sebetahat^2 + g_sd_choice[j]^2))

    lfsr_numerator_left = lfsr_numerator_left+ (1 - pi0)*g_prob_choice[j]*
      pnorm(0, g_mean_choice[j],g_sd_choice[j]  )*
      tnn_conv(v= betahat,s=sebetahat,mu= g_mean_choice[j], sigma= g_sd_choice[j], a = - Inf, b= 0)

    lfsr_numerator_right = lfsr_numerator_right+ (1 - pi0)*g_prob_choice[j]*
      pnorm(0, g_mean_choice[j],g_sd_choice[j], lower.tail = F )*
      tnn_conv(v= betahat, s= sebetahat,mu= g_mean_choice[j], sigma= g_sd_choice[j], a = 0, b=  Inf)

  }

  lfdr =lfdr_numerator/denominator
  lfsr = pmin(lfsr_numerator_left, lfsr_numerator_right)/denominator
  # get the direction
  lfsr_direction = rep(-1, m)
  lfsr_direction[which(lfsr_numerator_left < lfsr_numerator_right)] = 1

  ## sorting
  lfsr_ash_sort  = sort(lfsr_ash)
  lfdr_sort  = sort(lfdr)
  lfsr_sort  = sort(lfsr)
  FDR_est_lfdr = cumsum(lfdr_sort)/(1:m)
  FDR_est_lfsr = cumsum(lfsr_sort)/(1:m)
  FDR_est_lfsr_ash = cumsum(lfsr_ash_sort)/(1:m)
  rej_set_lfdr = which(rank(lfdr) <= sum(FDR_est_lfdr < alpha))
  rej_set_lfsr = which(rank(lfsr) <= sum(FDR_est_lfsr < alpha))
  rej_set_ash = which(rank(lfsr_ash) <= sum(FDR_est_lfsr_ash < alpha))

  return(list(betahat = betahat,
              sebetahat = sebetahat,
              beta = beta,
              lfdr = lfdr,
              lfsr = lfsr,
              lfsr_direction = lfsr_direction,
              lfsr_ash  = lfsr_ash,
              ash_direction= ash_direction,
              rej_set_lfdr= rej_set_lfdr,
              rej_set_lfsr = rej_set_lfsr,
              rej_set_ash = rej_set_ash))


}

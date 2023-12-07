# a function for convolving a normal distribution N(0, s^2)
# with a truncated normal TN(mu, sigma, a, b)

# based on the notes by Sebastien Turban


tnn_conv = function(v, s, mu, sigma, a, b){
  
  alpha =  s^2*(v- mu)/(s^2 + sigma^2)
  beta =  s*sigma/sqrt(s^2 + sigma^2)
  c = (mu - b)/sigma
  d = (mu - a)/sigma
  gamma = sqrt(2*pi)*beta/(2*pi*s*sigma*(pnorm(d) - pnorm(c)))
  
  value = exp(- (v- mu)^2/(2*(s^2 + sigma^2)))*gamma*(pnorm( (v- a - alpha)/beta) - pnorm((v- b- alpha)/beta))
  return(value)
  
  # value = - (v- mu)^2/(2*(s^2 + sigma^2)) + log (gamma) + log((pnorm( (v- a - alpha)/beta) - pnorm((v- b- alpha)/beta)))
  # return(exp(value))
}
# zdirect
An R package that implements ZDIRECT and the directional Storey, Taylor and Siegmund procedure for directional false discovery rate control.

## Installation
```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("ninhtran02/zdirect")
```

## Usage Examples
```
set.seed(123)
theta <- rnorm(1000,2,1)
theta[1:500] <- 0
sigma_vec <- rep(1,1000)
zv <- theta + rnorm(1000, sd = sigma_vec)
pv <- 2*pnorm(-abs(zv))


zdirect_obj <- zdirect(betahat= zv,
                       sebetahat = sigma_vec,
                       s_l = rep(0.2,1000), s_r =  rep(0.8,1000),
                       alpha = 0.1,
                       mixcompdist = "halfuniform",
                       prior = "nullbiased",
                       nullweight = 0.8,
                       altweight = 1e+08,
                       nfits = 200,
                       epsilon = 1e-10)
# Rejections
zdirect_obj$rej_set

# STS
STS_rej <- STS(pv = pv, alpha = 0.1, lambdas = 0.5)
STS_rej
```

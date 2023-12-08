# Set the library path

.libPaths("<R LIBRARY LOCATION>")

## Save personal library path as a variable

lib = .libPaths()[1]
#
### Install a package
#
#install.packages("REBayes",
#                 lib = lib,
#                 repos = "https://cran.ms.unimelb.edu.au/")
#
#install.packages("ashr",
#                 lib = lib,
#                 repos = "https://cran.ms.unimelb.edu.au/")
#devtools::install_github("lihualei71/dbh")


library("ashr", lib.loc = lib)
library("Matrix", lib.loc = lib)
library("REBayes", lib.loc = lib)
library("dbh", lib.loc = lib)

source("<MSKHOME>/mosek/10.0/tools/platform/<PLATFORM>/rmosek/builder.R")
#attachbuilder()
#install.rmosek()

library("Rmosek", lib.loc = lib)

source("get_lfsr_mask.R")
source("tnn_conv.R")  # used by bimodal
source("zdirect.R")
source("bimodal2.R")
source("dStorey.R")

cmd_args <- commandArgs(TRUE)
nullprop = as.numeric(cmd_args[1])
mu = as.numeric(cmd_args[2])
symm = as.numeric(cmd_args[3])
rho = as.numeric(cmd_args[4])
if(rho == 0){
  outfilename = paste0(getwd(),
                       "/simResults/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}
if( rho == 0.50){
  outfilename = paste0(getwd(),
                       "/simResultsrhopos05/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}
if( rho == 0.80){
  outfilename = paste0(getwd(),
                       "/simResultsrhopos08/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}
if( rho == -0.50){
  outfilename = paste0(getwd(),
                       "/simResultsrhoneg05/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}
if( rho == -0.80){
  outfilename = paste0(getwd(),
                       "/simResultsrhoneg08/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}


set.seed(1111)
N = 1000 # number of repetitions
print(paste("rho:",rho))
m = 1000  # number of tests
alpha = 0.1

################# PREPPING VECTORS FOR STORING THE FDR AND POWER #############
method_names = c("BH", "zdirect","lfsr",  "ash", "GR", "Storey", "StoreyAdaptive", "dBH")
save_list = c()
for (method in method_names){
  FDP_vec_name = paste0("FDP_", method)
  TPP_vec_name = paste0("TPP_", method)
  save_list = c(save_list, c(FDP_vec_name,TPP_vec_name ))
  assign(x = FDP_vec_name, value = numeric(0))
  assign(x = TPP_vec_name, value = numeric(0))
}
#################################################################################

#################  INITIALIZATIONS FOR ZDIRECT  ################
prior = "nullbiased"
# # one good combo
# nullweight = 1e+08
# altweight = 0.8

mixcompdist = "halfuniform"
nullweight = 0.8
altweight = 1e+08 # 0.8
nfits = 200
cutline = 0.4 # should be some number in (0, 0.5), recommened 0.2
s_l = rep(cutline/2,m)  # left starting threshold
s_r = rep(1 - cutline/2,m)  # right starting threshold
epsilon = 1e-10 # epsilon adjustment for zdirect
##################################################################

system.time(
  for (iter in 1:N){

    print(paste0("this is the ", iter, "-th iteration"))
    ######################## GENERATE DATA #########################
    sebetahat = sqrt(rep(1,m)/(1-rho^2))  # sebetahat = runif(m, min = 0.5, max = 10)   # sebetahat = sample(c(1, 1, 2), m, TRUE, c(1/3, 1/3, 1/3))
    data = bimodal(pi0 = nullprop,  sebetahat = sebetahat,
                   g_mean_choice = c(-mu, mu),
                   g_sd_choice = c(1, 1),
                   g_prob_choice = c(1-symm, symm), alpha = alpha, nullweight = 10, rho = rho)

    betahat = data$betahat
    beta = data$beta
    zv= betahat/sebetahat
    pv = 2*pnorm(-abs(zv))

    ################GET ALL THE REJECTION SETS#####################
    # oracle  rejection set
    rej_set_lfsr = data$rej_set_lfsr
    lfsr_direction = data$lfsr_direction
    # ash rejection set
    rej_set_ash = data$rej_set_ash
    ash_direction = data$ash_direction
    # BH rejection set
    rej_set_BH = which(p.adjust(pv, method = "BH")  <= alpha)
    # GR rejection set
    rej_set_GR = which(p.adjust(pv, method = "BH")  <= 2*alpha)

    # zdirect rejection set
    obj_zdirect = zdirect(betahat= betahat,
                          sebetahat = sebetahat,
                          s_l = s_l, s_r =  s_r,
                          alpha = alpha,
                          mixcompdist =mixcompdist,
                          prior = prior,
                          nullweight = nullweight,
                          altweight = altweight,
                          nfits = nfits,
                          epsilon = epsilon)
    rej_set_zdirect = obj_zdirect$rej_set


    ##############################################################
    # Storey Adaptive rejection set
    rej_set_StoreyAdaptive = STS(pv = pv, alpha = alpha,
                                 B = 1000, lambdas = seq(0.01,0.99,0.01))


    # Storey rejection set
    rej_set_Storey = STS(pv = pv, alpha = alpha,
                         B = 1000, lambdas = 0.5)

    # dBH rejection set
    Sigma <- (1/(1-rho^2))*rho^(abs(outer(1:m, 1:m, "-")))
    res <- dBH_mvgauss(zvals = betahat, Sigma = Sigma, side = "two", alpha = alpha,
                       gamma = 0.9, niter = 1, avals_type = "BH")
    rej_set_dBH = res$rejs

    FDP_lfsr= append( FDP_lfsr ,  sum(sign(beta[rej_set_lfsr]) != sign(lfsr_direction[rej_set_lfsr]) )/
                        max(1, length(rej_set_lfsr)))

    TPP_lfsr = append(TPP_lfsr , sum(sign(beta[rej_set_lfsr]) == sign(lfsr_direction[rej_set_lfsr]) )/
                        max(1, sum(beta!=0)))

    FDP_ash= append( FDP_ash ,  sum(sign(beta[rej_set_ash]) != sign(ash_direction[rej_set_ash]) )/
                       max(1, length(rej_set_ash)))

    TPP_ash = append(TPP_ash , sum(sign(beta[rej_set_ash]) == sign(ash_direction[rej_set_ash]) )/
                       max(1, sum(beta!=0)))

    for (method in c("BH",
                     "zdirect",
                     "GR",
                     "Storey",
                     "StoreyAdaptive",
                     "dBH")){
      FDP_vec_name = paste0("FDP_", method)
      TPP_vec_name = paste0("TPP_", method)
      rej_set_name = paste0( "rej_set_" , method)
      FDP = sum(sign(beta[eval(as.symbol(rej_set_name))]) != sign(betahat[eval(as.symbol(rej_set_name))]) )/
        max(1, length(eval(as.symbol(rej_set_name))))
      TPP =  sum(sign(beta[eval(as.symbol(rej_set_name))]) == sign(betahat[eval(as.symbol(rej_set_name))]) )/
        max(1, sum(beta!=0))
      assign(x = FDP_vec_name, value = append(eval(as.symbol(FDP_vec_name)) , FDP))
      assign(x = TPP_vec_name, value = append(eval(as.symbol(TPP_vec_name)) , TPP))
    }
  }
)

save( list= save_list,   file = outfilename)

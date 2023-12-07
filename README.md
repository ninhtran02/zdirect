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

# ZDIRECT
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
# ZDIRECT rejections
zdirect_obj$rej_set

# STS
STS_rej <- STS(pv = pv, alpha = 0.1, lambdas = 0.5)
STS_rej
```
## How to reproduce the simulation results for "Adaptive Procedures for Directional False Discovery Rate Control"
To reproduce the simulation results in an efficient manner, we assume the reader has access to an  account in a high performance computing (HPC) system running the \emph{Slurm Workload Manager}. Follow the steps below:

1. Log onto your own HPC account.

2. Create an R script called *main_bimodal2.R* as below. Ensure the <TEXT> are addressed according to the user.
```
# Set the library path

.libPaths(<R LIBRARY LOCATION>)

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
if( (rho == 0.50)||(rho == 0.8)){
  outfilename = paste0(getwd(),
                       "/simResults2/bimodal/nullprop", nullprop,
                       "mu", mu,
                       "symm", symm,
                       ".RData")
}
if( (rho == -0.50)||(rho == -0.8)){
  outfilename = paste0(getwd(),
                       "/simResults3/bimodal/nullprop", nullprop,
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
                       gamma = 1, niter = 1, avals_type = "BH")
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
```

4. Create a slurm file called ** as follows:
```
#!/bin/bash
#SBATCH --job-name=dFDR_bimodal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3500
#SBATCH --time=00-80:00:00

nullprop=$1
mu=$2
symm=$3
rho=$4


module --force purge
module load foss/2022a R/4.2.2


Rscript --vanilla main_bimodal2.R $nullprop $mu $symm $rho

```

Furthermore, create a slurm file called *batch_submission_bimodal.slurm* as follows:
```
#!/bin/bash


for nullprop in 0 0.2 0.5 0.8
do

  for mu in 0.5 1 1.5 2 2.5 3 
  do
  
    for symm in 0.5 0.75 1
    do
    	for rho in 0 0.50 -0.50
    	do
  
              sbatch /data/gpfs/projects/punim1426/ZDIRECT/job_submission_bimodal.slurm $nullprop $mu $symm $rho
	
  	done
    
      done
  
  done

done


```


4. Make sure the  *sn*, *ashr*, *dbh* and *Rmosek* R packages are properly installed in the version of R you plan to use. In particular, 

    - The installation instructions of \texttt{ashr} can be found in Matt Stephen's [GitHub page](https://github.com/stephens999/ashr).

    - The installation instructions of \texttt{Rmosek} can be found in [here]({https://docs.mosek.com/latest/rmosek/install-interface.html) . The *Rmosek* in turn depends on *MOSEK*, which requires a license; a free one for 365 days can be obtained [here](https://www.mosek.com/products/academic-licenses/) 
    

5. In your HPC account, change your current directory to the ZDIRECT folder using the ``cd" command:
```
cd    (your own working directory)/ZDIRECT
```

6. To submit the jobs, run the two files \texttt{batch\_submission\_skewnormal.slurm} and \texttt{batch\_submission\_bimodal.slurm} with the two commands:
\begin{center}
\texttt{./batch\_submission\_skewnormal.slurm}
\end{center}
\begin{center}
\texttt{./batch\_submission\_bimodal.slurm}
\end{center}

\item The simulations will typically be finished in a few ($\leq 5$) hours, with the resulting data files saved to the two subfolders within \texttt{ZDIRECT/simResults}.

\end{enumerate}


\section{How to reproduce the figures}

Assuming the simulation results have been reproduced as in the prior section:


\begin{enumerate}
\item Log onto your own HPC account.




\item In your HPC account, load  R using the command, say, 
\begin{center}
\texttt{module load r/4.0.0}
\end{center}
(You may need to change ``4.0.0" to ``4.1.0", if that is your own version of R)

\item Make sure the R package \texttt{ggpubr} has been properly installed in R.

\item Change your current directory to the ZDIRECT folder using the ``cd" command:
\begin{center}
\texttt{cd} \quad  (your own working directory)\texttt{/ZDIRECT}
\end{center}

\item Run the two command lines:
\begin{center}
\texttt{Rscript plot\_skewnormal.R}
\end{center}
\begin{center}
\texttt{Rscript plot\_bimodal.R}
\end{center}

\item The two pdf files for the figures will then appear in the ZDIRECT folder. 
\end{enumerate}

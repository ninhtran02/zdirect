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
To reproduce the simulation results in an efficient manner, we assume the reader has access to an  account in a high performance computing (HPC) system running the *Slurm Workload Manager8. Follow the steps below:

1. Log onto your own HPC account.

2. Create a slurm file called *job_submission_bimodal.slurm* as follows:
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
    	for rho in -0.8	-0.5 0 0.50 0.80
    	do
  
              sbatch /data/gpfs/projects/punim1426/ZDIRECT/job_submission_bimodal.slurm $nullprop $mu $symm $rho
	
  	done
    
      done
  
  done

done


```


3. Make sure the  *REBayes*, *ashr*, *dbh* and *Rmosek* R packages are properly installed in the version of R you plan to use. In particular, 

    - The installation instructions of *ashr* can be found in Matt Stephen's [GitHub page](https://github.com/stephens999/ashr).

    - The installation instructions of *Rmosek* can be found in [here]({https://docs.mosek.com/latest/rmosek/install-interface.html) . The *Rmosek* in turn depends on *MOSEK*, which requires a license; a free one for 365 days can be obtained [here](https://www.mosek.com/products/academic-licenses/) 
    

4. In your HPC account, change your current directory to the ZDIRECT folder using the ``cd" command:
```
cd    (your own working directory)/ZDIRECT
```
Make sure *main_bimodal2.R*, *job_submission_bimodal.slurm* and *batch_submission_bimodal.slurm* are inside the ZDIRECT folder.

5. To submit the jobs, run the two files *batch\_submission\_bimodal.slurm* with the command:
```
sbatch batch_submission_bimodal.slurm
```

6. The simulations will typically be finished in about a day, with the resulting data files saved to the subfolders within *ZDIRECT/simResults*, *ZDIRECT/simResultsrhopos05*, *ZDIRECT/simResultsrhopos08*, *ZDIRECT/simResultsrhoneg05* and *ZDIRECT/simResultsrhoneg08*.


## How to reproduce the figures for "Adaptive Procedures for Directional False Discovery Rate Control"

Assuming the simulation results have been reproduced as in the prior section:

1. Log onto your own HPC account and download the folders *ZDIRECT/simResults*, *ZDIRECT/simResultsrhopos05*, *ZDIRECT/simResultsrhopos08*, *ZDIRECT/simResultsrhoneg05* and *ZDIRECT/simResultsrhoneg08* locally onto your PC in a folder called *ZDIRECT*.

2. Run the script *plot_bimodal2.R* locally on your PC to produce the plots.


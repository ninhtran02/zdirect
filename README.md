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

2. Depending on the version of R available in your own HPC system, you may want to change the line
```
module load r/4.0.0
```
in the \texttt{job\_submission\_bimodal.slurm} and \texttt{job\_submission\_skewnormal.slurm} files to, say,
```
module load r/4.1.0
```
You can check the versions of R available in your system with the command 
```
module avail
```
3. Make sure the  \texttt{sn}, \texttt{ashr} and \texttt{Rmosek} R packages are properly installed in the version of R you plan to use. In particular, 

    - The installation instructions of \texttt{ashr} can be found in Matt Stephen's [GitHub page](https://github.com/stephens999/ashr).

    - The installation instructions of \texttt{Rmosek} can be found in \url{https://docs.mosek.com/latest/rmosek/install-interface.html}. The \texttt{Rmosek} in turn depends on \texttt{MOSEK}, which requires a license; a free one for 365 days can be obtained at 
    
```
\url{https://www.mosek.com/products/academic-licenses/}.
```

\item In your HPC account, change your current directory to the ZDIRECT folder using the ``cd" command:
\begin{center}
\texttt{cd} \quad  (your own working directory)\texttt{/ZDIRECT}
\end{center}

\item To submit the jobs, run the two files \texttt{batch\_submission\_skewnormal.slurm} and \texttt{batch\_submission\_bimodal.slurm} with the two commands:
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

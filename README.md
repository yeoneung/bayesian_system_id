#  On Bayesian Identification of Linear Non-Gaussian Systems: Tractable Algorithm and Concentration Bound

This document provides Matlab codes together with explicit setup for the experiments 
whose result is reported in the the paper above titled. The code is written based on Matlab2021b.

# Scenarios covered

### True System parameters

- Stabilizable $(A,B)$ 

$T=300$ and $N=200$ are used.

(i) $n=m=3$

$A={\left\lbrack \matrix{0.5 & 0 & 0.6 \cr 0.2 & 0 & 0.1 \cr 0 & 0.5 & 0.3} \right\rbrack}$, $B={\left\lbrack \matrix{0.4 & 0 & 0.4 \cr 0 & 0.3 & 0.1 \cr 0.3 & 0.2 & 0.1} \right\rbrack}$

(ii) $n=m=5$

$A={\left\lbrack \matrix{0.5 & 0 & 0.6 & 0.3 & 0.1 \cr 0.2 & 0 & 0.1 & 0.4 & 0.2  \cr 0 & 0.2 & 0 & 0.5 & 0.3 \cr 0.6 & 0.1 & 0.2 & 0 & 0.3 \cr 0.1 & 0.4 & 0.3 & 0.5 & 0.1 } \right\rbrack}$, 
$B={\left\lbrack \matrix{0.4 & 0 & 0.5 & 0 & 0.4 \cr 0 & 0.1 & 0.2 & 0.3 & 0.1  \cr 0.3 & 0 & 0 & 0.2 & 0.1 \cr 0.1 & 0.5 & 0.1 & 0.2 & 0.4 \cr 0.4 & 0 & 0.5 & 0.3 & 0.1 } \right\rbrack}$

- Unstabilizable $(A,B)$ with $n=m=3$ and $n=m=5$

$T=50$ and $N=100$ are used.

(i) $n=m=3$

$A={\left\lbrack \matrix{0 & -1  & 1 \cr 1 & 0 & 1 \cr 0 & 0 & 1} \right\rbrack}$, $B={\left\lbrack \matrix{1 & 0 & 0 \cr 0 & 0 & 1 \cr 0 & 0 & 0} \right\rbrack}$

(ii) $n=m=5$

$A={\left\lbrack \matrix{1 & 0 & 0 & 0 & 0 \cr 1 & 0 & 1 & 0 & -1  \cr 0 & 0 & -1 & 0 & 1 \cr 1 & 1 & 1 & 0 & 0 \cr 1 & 0 & -1 & 0 & 0 } \right\rbrack}$, 
$B={\left\lbrack \matrix{0 & 0 & 0 & 0 & 0 \cr 1 & 0 & 0 & 0 & 0  \cr 0 & 0 & 0 & 0 & 0 \cr 1 & 0 & 1 & 0 & 0 \cr 1 & 0 & 0 & 0 & 0 } \right\rbrack}$

### Type of system noise
Three different types of system noise are considered to test the performance of the algorithm.
(i) Gaussian
(ii) Gaussian mixture
(iii) asymmetric noise.

For (i) and (ii), built-in random generators in Matlab is used while we synthesize (iii) based on Langevin Markov Chain Monte Carlo (MCMC).

# Running

#### - Asymmetric noise generators
Running 'ULA_asymmetric_1D.m' generates a set of synthetic asymmetric noises based on Lagnevin Markov Chain Monte Carlo. 
This file will generate a set of 1000000 independent noises and save it 'asymmetric_noise_1D.csv'. 

#### - Ours


#### - LSE
We compare the performance of our algorithm with those proposed in [On the Sample Complexity of the Linear Quadratic Regulator](
https://link.springer.com/article/10.1007/s10208-019-09426-y). 
Running each m-file will provide an error for estimation corresponding to the type of noise indicated in the filename.

For example, running 'bayesian_system_id/LSE/Gaussian/Gaussian_3D.m' will create a plot for the error between the true and estimated parameters for the linear system with Gaussian system noise.









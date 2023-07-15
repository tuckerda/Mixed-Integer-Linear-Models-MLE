## Abbreviated User Guide:  milm_mle.m


__[x_hat,k_hat,pre] = milm_mle(Ac,yc,Mc,V,Prec,Ns,pre)__


This is the main routine for solving maximum likelihood estimation in mixed integer linear models,
$$y=Ax + Mk + u$$
where $y \in \mathbb{R}^m$ is a vector of noisy observations, $x\in\mathbb{R}^n$ are real-valued unknowns, $k \in \mathbb{Z}^m$ are integer-valued unknowns, and $u$ is zero-mean additive noise with inverse covariance matrix given by Prec. The matrix $A$ has full column rank, and $M$ is invertible; both $A$ and $M$ are assumed to have rational entries in this code.

The input $V$ is an $n \times n$ real-valued matrix providing a basis for the fundamental parallelotope describing the un-aliased values for $x$.  From $A$ and $M$, a basis $V$ is constructed by the utility routine, __form_Lambda_basis.m__.

Multiple values of $x$ and $k$ providing large likelihood scores are returned by choosing integer-valued input, Ns, greater than 1.

For applications in which an estimator is sought for many instances with the same $A$, $M$, $V$, and Prec, an input data structure, pre, is computed only once via milm_mle_precompute.m then re-used for all calls to milm_mle.


A special case of the problem occurs for multi-variate congruence equations 
$$y_1  \equiv a_{11} x_1 + a_{12} x_2 + a_{1n} x_n \mod b_1$$
$$y_2  \equiv a_{21} x_1 + a_{22} x_2 + a_{2n} x_n \mod b_2$$ 
$$\vdots$$
$$y_m \equiv a_{m1} x_1 + a_{m2} x_2 + a_{mn} x_n \mod b_m$$
The matrix $A$ above is formed as $[a_{ij}]$.  $M$ is the diagonal matrix of rational moduli, $b_1, ..., b_m$.  The remainders, $y_1, ... , y_m$ are available only as noisy versions in the column vector, $y$. The perturbations in $y$ are assumed zero-mean Gaussian, but may have unequal variance and may be correlated.


The following is an annotated list of the routines in the package.
* __milm_mle.m__     Main routine
* __Figure3_DoA.simulation.m__     Reproduces direction of arrival estimation example found in Figure 3 of the referenced _IEEE Signal Processing Letters_ paper.
  * __doa_pue.m__     Estimator from phase of the data covariance matrix
  * __doa_mle.m__     Maximum likelihood estimator from raw IQ data
* __form_Lambda_basis.m__     Computes lattice bases given $A$ and $M$
* __sils_reduction_Q.m__     For sphere decoding; redistributed, with modification, from X-C Wang, et al. 
* __sils_search.m__     For sphere decoding: redistributed from X-C Wang, et al. 
* Selected utilities
  * __lcms.m__     Finds column-wise positive least common multiple of the non-infinity, non-zero elements in a rational matrix.
  * __LLLReduce.m__			Implementation of lattice reduction algorithm described in  Wubben et al., 2004. 
  * __awgn_phase_covariance.m__	Constructs the approximate covariance matrix for phase-difference pairs of $L$ complex-valued observation with iid complex-Gaussian noise
  * __lenum.m__			Find the shortest vector of a lattice using the method by Schnorr & Euchner 1994, as implemented by Christian Chapman 2018.
  * __milm_mle_precompute.m__	Pre-compute quantities used by milm_mle
  * __tril_vec.m__		Returns vector of the lower triangular elements of a square matrix.
  * __wrapping_error_bounds.m__	Computes upper and lower bounds on the probably of incorrectly detecting the integer unknowns, $k.$
		

Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
All rights reserved.

This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree. 

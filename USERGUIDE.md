## Abbreviated User Guide
This repository contains MATLAB codes and simulation scripts that accompany the letter "Maximum Likelihood Estimation in Mixed Integer Linear Models." The letter considers a linear model with an observation vector $\mathbf{y} \in \mathbb{R}^m$ that is related to parameter vectors $\mathbf{x} \in \mathbb{R}^n$ and $\mathbf{k} \in \mathbb{Z}^m$ as

$$  \mathbf{y} =  \mathbf{A x} + \mathbf{M} \mathbf{k} + \mathbf{u},$$

with full column rank $\mathbf{A}\in \mathbb{R}^{m \times n}$, nonsingular $\mathbf{M} \in \mathbb{R}^{m \times m}$, and zero-mean Gaussian noise $\mathbf{u}  \in \mathbb{R}^m$. For brevity and simplicity, the observation vector $\mathbf{y}$ in the accompanying letter is assumed to have been whitened, and thus $\mathbf{u}$ has covariance matrix $\mathbf{I}$. 

The primary function in the repository is ```milm_mle```, a routine for maximum likelihood (ML) parameter estimation in mixed integer linear models. The function is called as

```[x_hat,k_hat,pre] = milm_mle(Ac,yc,Mc,V,Prec,Ns,pre)```

In the code in this repository, the observations are not assumed to be pre-whitened, and we therefore use the suffix "c" in the input variables name ```yc, Ac, Mc``` to denote the unwhitened versions of $y$, $A$, and $M,$ respectively, from the linear model above. The unwhitened noise has inverse covariance matrix ```Prec```. The matrix input ```Ac``` has full column rank, and ```Mc``` is invertible; both ```Ac``` and ```Mc``` must have rational entries.

The $n \times n$ real-valued matrix $\mathbf{V}$ provides a basis for $\Lambda$, the lattice that describes the periodicity of the likelihood in $\mathbf{x}$. A basis for $\Lambda$ is constructed as ```V = form_Lambda_basis(Ac,Mc).```

Note that multiple values of $\mathbf{x}$ and $\mathbf{k}$ providing large likelihood scores can be returned from ```milm_mle``` by choosing integer-valued input ```Ns``` that is greater than one.

For applications in which an estimator is sought for many instances with the same ```Ac```, ```M_c```, ```Vc,``` and ```Prec```, an input data structure ```pre``` can be precomputed via ```milm_mle_precompute```. This structure can be reused for all calls to ```milm_mle```.


The following is an annotated list of the routines included in this repository.
* ```milm_mle``` The main routine for maximum likelihood estimation in mixed integer linear models,
* ```Figure3_DoA.simulation``` Reproduces direction of arrival estimation example found in Figure 3 of the referenced _IEEE Signal Processing Letters_ paper.
  * ```doa_pue```     Estimates direction of arrival (DoA) with maximum likelihood phase unwrapping. 
  * ```doa_mle```     Maximum likelihood DoA estimation for linear/planar array by grid search using the complex data
* ```form_Lambda_basis```     Computes a basis for the lattice $\Lambda$
* ```sils_reduction_Q```     For sphere decoding; redistributed, with modification, from X.-W Chang, et al. 
* ```sils_search```     For sphere decoding: redistributed from X.-W Cang, et al. 
* Selected utilities
  * ```lcms.m```     Finds column-wise positive least common multiple of the non-infinity, non-zero elements in a rational matrix.
  * ```LLLReduce``` Implementation of the Lenstra-Lenstra-Lov\'{a}sz (LLL) lattice basis reduction algorithm of Lenstra et al., 1982, as described in  Wubben et al., 2004. 
  * ```awgn_phase_covariance```	Constructs the approximate covariance matrix for phase-difference pairs of $L$ complex-valued observation with iid complex-Gaussian noise
  * ```lenum```			Find the shortest vector of a lattice using the method by Schnorr & Euchner 1994, as implemented by Christian Chapman 2018. See https://github.com/enthdegree/lenum.m
  * ```milm_mle_precompute```	Pre-compute quantities used by ```milm_mle```
  * ```tril_vec```		Returns vector of the lower triangular elements of a square matrix.
  * ```wrapping_error_bounds```	Computes upper and lower bounds on the probably of incorrectly detecting the integer unknowns, $\mathbf{k}.$
		

Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
All rights reserved.

This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree. 

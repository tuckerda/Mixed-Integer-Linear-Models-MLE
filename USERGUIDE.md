# Abbreviated User Guide

This repository contains MATLAB code and simulation scripts for the paper "Maximum Likelihood Estimation in Mixed Integer Linear Models." In the paper, we consider a linear model with an observation vector $\mathbf{y} \in \mathbb{R}^m$ that is related to parameter vectors $\mathbf{x} \in \mathbb{R}^n$ and $\mathbf{k} \in \mathbb{Z}^m$ as

$$  \mathbf{y} =  \mathbf{A x} + \mathbf{M} \mathbf{k} + \mathbf{u},$$

with full column rank $\mathbf{A}\in \mathbb{R}^{m \times n}$, nonsingular $\mathbf{M} \in \mathbb{R}^{m \times m}$, and zero-mean Gaussian noise $\mathbf{u}  \in \mathbb{R}^m$. For brevity, the observation vector $\mathbf{y}$ in the accompanying letter is assumed to have been whitened, and thus $\mathbf{u}$ has covariance matrix $\mathbf{I}$. 

## A Routine for Maximum Likelihood Parameter Estimation in Mixed Integer Linear Models
The primary function in the repository is `milm_mle`, a routine for maximum likelihood (ML) parameter estimation in mixed integer linear models. The function is called as

```[x_hat,k_hat,pre] = milm_mle(Ac,yc,Mc,V,Prec,Ns,pre)```

In this repository, the observations are not assumed to be pre-whitened, and we therefore use the suffix "c" in the input variable names `yc,` `Ac,` and `Mc` to denote correspondence to the unwhitened versions of $\mathbf{y}$, $\mathbf{A}$, and $\mathbf{M}$ from the linear model above, respectively. The input `Prec` is the inverse covariance matrix of the unwhitend noise. The matrix input `Ac` has full column rank, and `Mc` is invertible; both `Ac` and `Mc` must have rational entries. The input `V` is a rational $n \times n$ matrix that provides a basis for $\Lambda$, the lattice that describes the periodicity of the likelihood in $\mathbf{x}$. A basis for $\Lambda$ can be constructed as `V = form_Lambda_basis(Ac,Mc).`

Note that multiple values of $\mathbf{x}$ and $\mathbf{k}$ providing large likelihood scores can be returned from `milm_mle` by choosing integer-valued input `Ns` that is greater than one. 

For applications in which an estimator is sought for many instances with the same `Ac`, `Mc,` and `Prec`, an input data structure can be precomputed via `milm_mle_precompute`. This structure can then be used repeatedly for all calls to `milm_mle`.

## Function list
The following is an annotated list of the routines included in this repository.
* `milm_mle` Routine for maximum likelihood estimation in mixed integer linear models,
* `Figure3_DoA_simulation` Reproduces direction of arrival estimation example found in Figure 3 of the referenced _IEEE Signal Processing Letters_ paper.
  * `doa_pue`     Estimates direction of arrival (DoA) with maximum likelihood phase unwrapping. 
  * `doa_mle`     Maximum likelihood DoA estimation for linear/planar array by grid search using the complex data
* `form_Lambda_basis`     Computes a basis for the lattice $\Lambda$
* Selected utilities
  * `sils_reduction_Q`     Lattice basis reduction, redistributed, with modification, from [MILES package by X.-W Chang, et al.](https://www.cs.mcgill.ca/~chang/MILES_routine1.php)
  * `sils_search`     Depth-first sphere decoding, redistributed from [MILES package by X.-W Chang, et al.](https://www.cs.mcgill.ca/~chang/MILES_routine1.php) 
  * `lcms.m`     Finds column-wise positive least common multiple of the non-infinity, non-zero elements in a rational matrix.
  * `LLLReduce` Implementation of the Lenstra-Lenstra-Lovasz (LLL) lattice basis reduction algorithm of Lenstra et al., 1982, as described in  Wubben et al., 2004. 
  * `awgn_phase_covariance`	Constructs the approximate covariance matrix for phase-difference pairs of $L$ complex-valued observation with iid complex-Gaussian noise
  * `lenum`			Find the shortest vector of a lattice using the method by Schnorr & Euchner 1994, redistributed from [implementation by Christian Chapman](https://github.com/enthdegree/lenum.m). 
  * `milm_mle_precompute`	Pre-compute quantities used by `milm_mle`
  * `tril_vec`		Returns vector of the lower triangular elements of a square matrix.
  * `wrapping_error_bounds`	Computes upper and lower bounds on the probably of incorrectly detecting the integer unknowns, $\mathbf{k}.$

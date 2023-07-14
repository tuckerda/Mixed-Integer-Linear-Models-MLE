function [x_hat,k_hat,pre] = milm_mle(Ac,yc,Mc,V,Prec,Ns,pre)
% Finds maximum likelihood estimate of unknown {x,k} in the mixed integer linear model
%               yc = Ac*x + Mc*k + uc, 
% where uc follows a zero-mean Gaussian distribution.
%
% Input:
%   Ac: Unwhitened rational, full column rank m-by-n matrix
%   yc: Real m-by-1 vector of noisy, wrapped, unwhitened observations
%   Mc: Unwhitened rational, nonsingular p-by-p matrix
%   V: Basis for {x | Mc^(-1)*Ac*x is integer-valued }, e.g., from form_Lambda_basis
%   Prec: Precision matrix of zero-mean Gaussian noise uc
%   Ns: Number of solutions to return
%   pre: Optional struct of precomputed quantities.
%
% Output:
%   x_hat: Optimal solution for real n-by-1 x
%   k_hat: Optimal solution for integer m-by-1 k
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 

[m,n] = size(Ac);
% Compute quantitites to be reused
if (nargin < 7 || isempty(pre))
    pre = milm_mle_precompute(Ac,Mc,V,Prec);
end
y = pre.Whiten*yc;
% Solve min_k || ytilde - B*k ||, with  ytilde = B*M^(-1)*y
if m > n
    % Solve min_ell || Q'*ytilde - R*Z^(-1)*ell ||
    Qtytilde = pre.QtBMinv*y;
    ell_Babai = round(pre.Z*inv(pre.R)*Qtytilde);
    % Use the Babai point, if optimal according to Hassibi/Boyd 1998 eqn. (25). Otherwise, use sphere decoding
    if norm(Qtytilde - pre.R*ell_Babai) <= pre.d_min / 2
        ell_star = ell_Babai;
    else
        ell_star = pre.Z * sils_search(pre.R, Qtytilde, Ns);
    end
    k_star = pre.Gperp*ell_star;
else
    k_star = zeros([p,1]);
end
x_star =  pre.Apinv*(y - pre.M*k_star);
% Wrap x_star to parallelotope defined by V
x_hat = x_star - pre.V*floor(pre.Vinv*x_star);
k_hat = k_star - pre.Minv*pre.A*(x_hat - x_star);
end
function [x_mlpue,pre] = doa_pue(Ac,yc,Mc,V,Prec,pre)
% Estimates direction of arrival (DoA) with maximum likelihood phase unwrapping.
% Solves a closest lattice point problem and projects the solution to the set of feasible DoA.
%
% Input:
%   Ac: Unwhitened rational, full column rank m-by-n matrix
%   yc: Real m-by-1 vector of noisy, wrapped, unwhitened observations
%   Mc: Unwhitened rational, nonsingular p-by-p matrix
%   V: Basis for {x | Mc^(-1)*Ac*x is integer-valued }, e.g., from form_Lambda_basis
%   Prec: Precision matrix of zero-mean Gaussian noise vector uc
%   pre: Optional struct of precomputed quantities.
%
% Output:
%   x_mlpue: Optimal solution for real n-by-1 x
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 

[m,n] = size(Ac);
if ~(n==1||n==2||n==3)
    error("The parameter size n is invalid.");
end

% Find a particular solution to unconstrained problem
[x_hat,~,pre] = milm_mle(Ac,yc,Mc,V,Prec,1,pre);
if ~isfield(pre,'K')
    pre = doa_pue_precompute(pre);
end
% Wrap the solution to a parallelotope centered on the origin
x_hat = x_hat - pre.V*floor(pre.Vinv*x_hat + 1/2);

% Project the unconstrained solution to the feasible set
if n==1
    x_mlpue = x_hat;
    x_mlpue(abs(x_mlpue) > 1) = sign(x_mlpue(abs(x_mlpue) > 1));
else % (n==2||n==3)
    x_shifted = x_hat - V*pre.K;
    if n==2
        x_proj = x_shifted./max(1,vecnorm(x_shifted,2,1));
    else % n==3
        x_proj = x_shifted./vecnorm(x_shifted,2,1);
    end
    distance = x_shifted - x_proj;
    cost = sum(distance.^2,1);
    [~, idx] = min(cost);
    x_mlpue = x_proj(:,idx);
end
end

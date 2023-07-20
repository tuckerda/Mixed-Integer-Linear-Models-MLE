function [pre] = milm_mle_precompute(Ac,Mc,V,Prec)
% Compute quantities used by milm_mle

% Input:
%   Ac: Unwhitened rational, full column rank m-by-n matrix
%   Mc: Unwhitened rational, nonsingular p-by-p matrix
%   V: Basis for {x | Mc^(-1)*Ac*x is integer-valued }, e.g., from form_Lambda_basis
%   Prec: Precision matrix of zero-mean Gaussian noise u

% Output:
%   pre: Struct of precomputed quantities.
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 

[m,n] = size(Ac);

%[V, Gperp] = form_Lambda_basis(Ac,Mc);
D = diag(lcms(sym(1./(inv(Mc)*Ac))));
T = round(inv(Mc)*Ac*D); % Integer matrix
% H = U*T is Hermite normal form of T=Mc^(-1)*Ac*D, with upper-triangular H and unimodular U
[U,H] = hermiteForm(T);
Gperp = inv(U); Gperp = Gperp(:,(n+1):m);
%V = D*inv(H(1:n,:));

Whiten = sqrtm(Prec);
A = Whiten*Ac;
M = Whiten*Mc;
% B is m-by-m, generates (m-n)-dimensional lattice Lambda_B, range(B)=kern(A^T)
B = (eye(m) - A*pinv(A))*M;
% C is a lattice basis for Lambda_B
C = B*Gperp;
% Find a reduced basis for Lambda_B
[Q,R,Z] = sils_reduction_Q(C); % Q*R = C*Z = B*Gperp*Z
if m > n
    % Find the shortest vector in Lambda_B
    d_min = norm(lenum(R));
else
    d_min = Inf;
end
Prec_ml = A'*A;

pre = struct( ...
    'Q',Q, ...
    'R',R, ...
    'Z',Z, ...
    'C',C, ...
    'Gperp',Gperp, ...
    'V',V, ...
    'M',M, ...
    'A',A, ...
    'B',B, ...
    'QtBMinv', Q'*B*inv(M), ...
    'Vinv',inv(V), ...
    'Apinv',pinv(A), ...
    'Minv',inv(M), ...
    'd_min',d_min, ...
    'Whiten',Whiten, ...
    'Prec_ml',Prec_ml, ...
    'Rinv',inv(R), ...
    'GperpZ',Gperp*Z ...
    );
end

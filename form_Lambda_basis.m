function [V,Gperp] = form_Lambda_basis(Ac,Mc)
% Finds a basis for {x | Mc^(-1)*Ac*x is integer-valued}
%
% Input:
%   Ac: Rational m-by-n matrix with full column rank.
%   Mc: Rational, nonsingular p-by-p matrix
%
% Output:
%   V: Basis for {x | Mc^(-1)*Ac*x is integer-valued}
%   Gperp: Integer matrix, the last m-n columns of U^(-1)
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 

[m,n] = size(Ac);
D = diag(lcms(sym(1./(inv(Mc)*Ac))));
T = round(inv(Mc)*Ac*D); % Integer matrix

% H = U*T is Hermite normal form of T=Mc^(-1)*Ac*D, with upper-triangular H and unimodular U
[U,H] = hermiteForm(T);
Gperp = inv(U); Gperp = Gperp(:,(n+1):m);
V = D*inv(H(1:n,:));
end
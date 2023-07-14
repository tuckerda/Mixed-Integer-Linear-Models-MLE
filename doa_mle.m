function [omega_ml] = doa_mle(Ry,X,omega_scan)
% Maximum likelihood DoA estimation for linear/planar array by searching over increasingly fine grids
%
% Input:
%	Ry: Sample covariance matrix
%	X: Element positions of L-element array, in units of half-wavelengths, in n-dimensional space, [L,n]
%   omega_scan: Cell of search grids
%
% Output:
%	omega_ml: Maximum likelihood estimate of omega, [n,1]
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.

n = size(X,2);
if ~(n==1||n==2)
    error("The parameter size n is invalid.");
end
if ~iscell(omega_scan)
    Nres = 1;
    scan_current = omega_scan;
else
    Nres = length(omega_scan);
    scan_current = omega_scan{1};
end
for k = 1:Nres
    % Form steering vectors for candidate directions
    A = exp(1i*pi*X*scan_current);
    % Evaluate candidate directions
    NLL = -real(sum(conj(A).*(Ry*A)));
    [~,min_idx] = min(NLL);
    omega_ml = scan_current(:,min_idx);

    if k < Nres
        scan_current = omega_ml + omega_scan{k+1};
        if n==2
            scan_current(:, vecnorm(scan_current,2,1) > 1 ) = [];
        else %n==1
            scan_current( abs(scan_current) > 1 ) = [];
        end
    end
end
end

% grid_res = [ 12e-3, 4e-3, 1e-3];
% omega_scan = cell([length(grid_res),1]);
% for k=1:length(grid_res)
%     if k==1
%         omega_bounds = [-ones([n,1]), ones([n,1])];
%     else
%         omega_bounds = [-grid_res(k-1)*ones([n,1]), grid_res(k-1)*ones([n,1])];
%     end
%     omega_grid = cell(n,1);
%     for i = 1:n
%         omega_grid(i) = { linspace( omega_bounds(i,1), omega_bounds(i,2), ceil( ( omega_bounds(i,2) - omega_bounds(i,1))/grid_res(k) + 1 ) ) };
%     end
%     tmp = cell(size(omega_grid));
%     [tmp{:}] = ndgrid(omega_grid{:});
%     tmp = reshape(cat(n+1,tmp{:}),[],n).';
%     tmp(:, vecnorm(tmp,2,1) > 1 ) = [];
%     omega_scan{k} = tmp;
% end
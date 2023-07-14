function [CRB,CRB_s] = crb_c(X,N,aoa,s,n_var)
% Cramér–Rao bound for DoA estimation under the conditional (deterministic) signal model.
%
% Input:
%	X:   Sensor positions (in units of half-wavelengths) of L-element array in n-dimensional space (n=1,2,3), [L n]
%   N: Number of snapshots
%   aoa: Angles of arrival, in radians, of K sources, [1 K] for linear arrays or [2 K] for planar/volume arrays
%        % Angles specified as [ azimuth (from y-axis), elevation ]
%   s : Signal s, [K N]
%   n_var: Noise variance, [1 1]
%
% Output:
%   CRB: CRB for angle of arrival in radians
%   CRB_s: CRB for DoA parameter omega
%
% [1] Mailaender, L. Bounds for 2-D angle-of-arrival estimation with separate and joint processing.
%     Eurasip J. Adv. Signal Process. 2011. 1-11. 10.1186/1687-6180-2011-5.
% [2, p. 372] Petre Stoica and Randolph Moses, Spectral Analysis of Signals, Prentice- Hall, Upper Saddle River, NJ, 2005.
% [3] P. Stoica and A. Nehorai, "Performance study of conditional and unconditional direction-of-arrival estimation,"
%     IEEE Tran. Acoust. Speech and Signal Process., vol. 38, no. 10, pp. 1783-1795, Oct. 1990, doi: 10.1109/29.60109.

[L,n] = size(X);
K = size(aoa,2);

% If s is size [K,1], assume equal signal power for all snapshots
if (N>1 && numel(s)==K)
    s = reshape(s,[],1);
    s = repmat(s,[1 N]);
end

[A,dA,dA_s] = array_steering(X,aoa);
P = 1/N*(s*s');
Proj = eye(L) - A*inv(A'*A)*A';
CRB = n_var/N/2 * inv( real( dA'*Proj*dA .* (kron(ones([n,n]),P)) ) );
CRB_s = n_var/N/2 * inv( real( dA_s'*Proj*dA_s .* (kron(ones([n,n]),P)) ) );
end
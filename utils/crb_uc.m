function [CRB,CRB_s] = crb_uc(positions,N,aoa,P,n_var)
% Cramér–Rao bound for DoA estimation under the unconditional (stochastic) signal model.
%
% Input:
%   positions: Sensor positions (in units of half-wavelengths) of L-element array in n-dimensional space (n=1,2,3), [L n]
%   N: Number of snapshots
%   aoa: Angles of arrival, in radians, of K sources, [1 K] for linear arrays or [2 K] for planar/volume arrays
%        % Angles specified as [ azimuth (from y-axis), elevation ]
%   P: Signal covariance, [K K]
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

[L,n] = size(positions);
K = size(aoa,2);

[A,dA,dA_s] = array_steering(positions,aoa);
U = P*inv(A'*A*P + n_var*eye(K))*A'*A*P;
Proj = eye(L) - A*inv(A'*A)*A';
CRB = n_var/N/2 * inv( real( dA'*Proj*dA .* (kron(ones([n,n]),U)) ) );
CRB_s = n_var/N/2 * inv( real( dA_s'*Proj*dA_s .* (kron(ones([n,n]),U)) ) );
end

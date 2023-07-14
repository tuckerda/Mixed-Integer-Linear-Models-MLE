function [Sigma] = awgn_phase_covariance(pair_idx,SNR,N)
% Approximate covariance of phase-difference pairs from L-element array output with N i.i.d. snapshots.
% There is zero-mean i.i.d. Gaussian noise at the array outputs that is uncorrelated across snapshots.
%
% Input:
%   pair_idx: Index pairs from L-element array, e.g., from nchoosek(1:L,2)
%	SNR: Signal-to-noise ratio of the array output (linear)
%   N: Number of snapshots
%
% Output:
%   Sigma: Approximate phase noise covariance matrix for the L*(L-1)/2 sensor pairs

m = size(pair_idx,1);
Sigma = zeros([m,m]);
for k1=1:m
    for k2=k1+1:m
        if (pair_idx(k1,1)==pair_idx(k2,1) || pair_idx(k1,2)==pair_idx(k2,2))
            Sigma(k1,k2) = 1;
        elseif (pair_idx(k1,1)==pair_idx(k2,2) || pair_idx(k1,2)==pair_idx(k2,1))
            Sigma(k1,k2) = -1;
        end
    end
end
Sigma = Sigma + Sigma';
Sigma = 1/(2*N*SNR) * Sigma;
Sigma = Sigma + (1+2*SNR)/(2*N*SNR^2) * eye(m);
end
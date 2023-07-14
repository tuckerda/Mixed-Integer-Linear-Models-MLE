function [Q,R,H] = LLLReduce(Qin,Rin,delta)
% Lenstra–Lenstra–Lovász lattice basis reduction. See Table 1 in [1].
% The reduced basis B is given by B=Q*R=A*H, where A=Qin*Rin is the original basis.
%
% Input:
%    Qin - m-by-m orthonormal matrix with A=Qin*Rin
%    Rin - m-by-n upper triangular matrix with A=Qin*Rin
% Output:
%    Q - m-by-m LLL-reduced orthonormal matrix
%    R - m-by-n LLL-reduced upper triangular matrix
%    H - n-by-n unimodular matrix

% [1] D. Wubben, R. Bohnke, V. Kuhn and K.-D. Kammeyer, "Near-maximum-likelihood detection
% of MIMO systems using MMSE-based lattice-reduction," IEEE Int. Conf. Commun., 2004.

if nargin==1
    [Qin,Rin] = qr(Qin);
elseif nargin ==2
    delta = 0.75;
end

n = size(Rin,2); %[m,n] = size(A);
Q = Qin;
R = Rin;
H = eye(n);
k = 2;
while k <= n
    % Size reduction
    for j = (k-1):-1:1
        r = round(R(j,k) / R(j,j));
        if r ~=0
            % Update affected R values
            R(1:j,k) = R(1:j,k) - r * R(1:j,j);
            H(:,k) = H(:,k) - r*H(:,j);
        end
    end
    % Check Lovasz condition
    if delta*R(k-1,k-1)^2 > R(k-1,k)^2 + R(k,k)^2
        H(:,[k-1,k]) = H(:,[k,k-1]);
        R(:,[k-1,k]) = R(:,[k,k-1]);

        %  Givens rotation (from sils.m: https://www.cs.mcgill.ca/~chang/MILES_routine1.php)
        [G, R([k-1,k], k-1)] = planerot(R([k-1,k], k-1));
        R([k-1,k], k:n) = G*R([k-1, k],k:n);
        Q(:, [k-1,k]) = Q(:, [k-1, k])*G';

        k = max(2, k-1);
    else
        k = k+1;
    end
end
end

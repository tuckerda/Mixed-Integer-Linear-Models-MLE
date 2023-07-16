function [Q,R,Z] = sils_reduction_Q(B)
%
% [Q,R,Z] = sils_reduction_Q(B) reduces the general standard integer
% least squares problem to an upper triangular one by the LLL-QRZ
% factorization Q'*B*Z = [R; 0]. This version of sils_reduction
% has been modified to output the orthonormal Q matrix.
%
%
% Inputs:
%    B - m-by-n real matrix with full column rank
%
% Outputs:
%    Q - m-by-n orthonormal matrix
%    R - n-by-n LLL-reduced upper triangular matrix
%    Z - n-by-n unimodular matrix, i.e., an integer matrix with |det(Z)|=1

% Subfunction: qrmcp

% Main Reference:
% X. Xie, X.-W. Chang, and M. Al Borno. Partial LLL Reduction,
% Proceedings of IEEE GLOBECOM 2011, 5 pages.

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiaohu Xie, Tianyang Zhou
%
%
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% October 2006. Last revision: June 2016
%
% Modification of sils_reduction.m from https://www.cs.mcgill.ca/~chang/MILES_routine1.php
% Modified by David Tucker in May 2023

[~,n] = size(B);

% QR factorization with minimum-column pivoting
[Q,R,piv] = qrmcp(B);

% Obtain the permutation matrix Z
Z = zeros(n,n);
for j = 1 : n
    Z(piv(j),j) = 1;
end

% ------------------------------------------------------------------
% --------  Perfome the partial LLL reduction  ---------------------
% ------------------------------------------------------------------

k = 2;

while k <= n

    k1 = k-1;
    zeta = round(R(k1,k) / R(k1,k1));
    alpha = R(k1,k) - zeta * R(k1,k1);

    if R(k1,k1)^2 > (1 + 1.e-10) * (alpha^2 + R(k,k)^2)
        if zeta ~= 0
            % Perform a size reduction on R(k-1,k)
            R(k1,k) = alpha;
            R(1:k-2,k) = R(1:k-2,k) - zeta * R(1:k-2,k-1);
            Z(:,k) = Z(:,k) - zeta * Z(:,k-1);

            % Perform size reductions on R(1:k-2,k)
            for i = k-2:-1:1
                zeta = round(R(i,k)/R(i,i));
                if zeta ~= 0
                    R(1:i,k) = R(1:i,k) - zeta * R(1:i,i);
                    Z(:,k) = Z(:,k) - zeta * Z(:,i);
                end
            end
        end

        % Permute columns k-1 and k of R and Z
        R(1:k,[k1,k]) = R(1:k,[k,k1]);
        Z(:,[k1,k]) = Z(:,[k,k1]);

        % Bring R back to an upper triangular matrix by a Givens rotation
        [G,R([k1,k],k1)] = planerot(R([k1,k],k1));
        R([k1,k],k:n) = G * R([k1,k],k:n);

        Q(:,[k1,k]) =  Q(:,[k1,k])*G'; % DT

        if k > 2
            k = k - 1;
        end

    else
        k = k + 1;
    end
end
end


function [Q,R,piv] = qrmcp(B)
%
% [R,piv,y] = qrmcp(B,y) computes the QR factorization of B with
%             minimum-column pivoting:
%                  Q'BP = R (underdetermined B),
%                  Q'BP = [R; 0] (underdetermined B)
%             and computes Q'*y. The orthogonal matrix Q is not produced.
%
% Inputs:
%    B - m-by-n real matrix to be factorized
%    y - m-dimensional real vector to be transformed to Q'y
%
% Outputs:
%    R - m-by-n real upper trapezoidal matrix (m < n)
%        n-by-n real upper triangular matrix (m >= n)
%    piv - n-dimensional permutation vector representing P
%    y - m-vector transformed from the input y by Q, i.e., y := Q'*y

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiaohu Xie, Tianyang Zhou
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% October 2006; Last revision: June 2016


[m,n] = size(B);

% Initialization
colNormB = zeros(2,n);
piv = 1:n;
Q = eye(m); % DT

% Compute the 2-norm squared of each column of B
for j = 1:n
    colNormB(1,j) = (norm(B(:,j)))^2;
end

n_dim = min(m-1,n);

for k = 1 : n_dim
    % Find the column with minimum 2-norm in B(k:m,k:n)
    [~, i] = min(colNormB(1,k:n) - colNormB(2,k:n));
    q = i + k - 1;

    % Column interchange
    if q > k
        piv([k,q]) = piv([q,k]);
        colNormB(:,[k,q]) = colNormB(:,[q,k]);
        B(:,[k,q]) = B(:,[q,k]);
    end

    % Compute and apply the Householder transformation  I-tau*v*v'
    if norm(B(k+1:m,k)) > 0 % A Householder transformation is needed
        v = B(k:m,k);
        rho = norm(v);
        if v(1) >= 0
            rho = -rho;
        end
        v(1) = v(1) - rho; % B(k,k)+sgn(B(k,k))*norm(B(k:n,k))
        tao = -1 / (rho * v(1));
        B(k,k) = rho;
        if m < n
            B(k+1:m,k) = 0;
        end
        B(k:m,k+1:n) = B(k:m,k+1:n) - tao * v * (v' * B(k:m,k+1:n));
        Q(:,k:end) = Q(:,k:end)-(Q(:,k:end)*v)*(tao*v)'; % DT

    end

    % Update colnormB(2,k+1:n)
    colNormB(2,k+1:n) = colNormB(2,k+1:n) + B(k,k+1:n) .* B(k,k+1:n);
end

if m < n
    R = B;
else
    R = triu(B(1:n,1:n));
    Q = Q(:,1:n); % DT
end
end

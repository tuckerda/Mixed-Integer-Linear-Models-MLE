function [v] = tril_vec(A)
% Return vector v of lower triangular elements of NxN matrix A.
% v = [A_{2,1},A_{3,1},A_{4,1},...,A_{L,L-1}]^T
if ~(ismatrix(A) && size(A,1)==size(A,2))
    error('Input must be square matrix.')
end
N = size(A,1);
v = zeros([N*(N-1)/2,1]);
idx = 1;
for k2=1:N
    for k1 = k2+1:N(1)
        v(idx) = A(k1,k2);
        idx = idx + 1;
    end
end
end
function output = lcms(A)  
% Finds column-wise positive least common multiple of non-infinity,
% non-zero elements of rational matrix A.
%
% Inputs:
%   A       Rational p-by-q matrix.
% Outputs:
%   output  q-by-1 list
%
% SZ. May 2023

% Column-wise least common multiple
output = zeros(size(A,2),1);
for i = 1:size(A,2) %column by column
    Tmp = abs(A(:,i));
    Tmp(isinf(Tmp)) = [];%check for INF
    Tmp(Tmp==0) = [];%check for zero
    output(i) = Tmp(1);
    for n = 2:numel(Tmp)
        output(i) = lcm(output(i),Tmp(n));
    end
end

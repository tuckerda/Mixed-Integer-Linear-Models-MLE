function [pre] = doa_pue_precompute(pre)
[m,n] = size(pre.A);
if (n==2)
    pre.K = [ -1,0,1,-1,0,1,-1,0,1; -1,-1,-1,0,0,0,1,1,1];
elseif (n==3)
    C_Grid = cell(n,1);
    for i = 1:n
        C_Grid(i) = {[-1,1]};
    end
    C = cell(size(C_Grid));
    [C{:}] = ndgrid(C_Grid{:});
    C = reshape(cat(n+1,C{:}),[],n).';

    K_abs_max = max( ceil( abs( inv(1/2*pre.V)*C) ),[],2);
    K_Grid = cell(n,1);
    for i = 1:n
        K_Grid(i) = {-K_abs_max(i):K_abs_max(i)};
    end
    K = cell(size(K_Grid));
    [K{:}] = ndgrid(K_Grid{:});
    K = reshape(cat(n+1,K{:}),[],n).';
    pre.K = K;
end
%[pre.Q_M,pre.R_M,pre.Z_M] = sils_reduction_Q(pre.M); % Q*R = M*Z 
end

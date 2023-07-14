function [Pc_lb,Pc_ub] = wrapping_error_bounds(Ac,Mc,Prec)
% [1] A. Hassibi and S. Boyd, “Integer parameter estimation in linear models
% with applications to GPS," IEEE Trans. Signal Process., vol. 46, no. 11, pp. 2938–2952, 1998

[m,n] = size(Ac);
D = diag(lcms(sym(1./(inv(Mc)*Ac))));
T = round(inv(Mc)*Ac*D); % Integer matrix
% H = U*T is Hermite normal form of T=Mu^(-1)*Au*D, with upper-triangular H and unimodular U
[U,H] = hermiteForm(T);
Gperp = inv(U); Gperp = Gperp(:,(n+1):m);

Whiten = sqrtm(Prec);
A = Whiten*Ac;
M = Whiten*Mc;
B = (eye(m) - A*pinv(A))*M;
C = B*Gperp;

% Pc_lb lower-bounds Pc, the probability of detecting the integer parameter correctly (see [1,eqn. (20)])
d_min = norm(lenum(C));
Pc_lb = chi2pdf(d_min^2/4, m-n);

detC = sqrt( det(C'*C) );
Alpha = pi^((m-n)/2) / gamma(1+(m-n)/2);
% Pc_ub upper-bounds Pc (see [1,eqn. (18)])
Pc_ub = chi2pdf( (detC/Alpha)^(2/(m-n)), m-n);
end
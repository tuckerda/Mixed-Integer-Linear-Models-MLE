function [A,dA,dA_s] = array_steering(X,aoa)
% Generate array steering vectors
%
% Input
%	X: Sensor positions (in units of half-wavelengths) of L-element array in n-dimensional space (n=1,2,3), [L n]
%   aoa: Angles of arrival, in radians, of K sources, [1 K] for linear arrays or [2 K] for planar/volume arrays
%        % Angles specified as [ azimuth (from y-axis), elevation ]
%
% Output
%	A:   Steering vectors, [L K]
%   dA: Partial derivatives of A w.r.t. angles of arrival. Organized as in eqn. 3.8 in [1]
%   dA_s: Partial derivatives of A w.r.t. DoA parameter omega. Organized as in eqn. 3.8 in [1]

% David Tucker, June 2023
%
% [1] Mailaender, L. Bounds for 2-D angle-of-arrival estimation with separate and joint processing.
%     Eurasip J. Adv. Signal Process. 2011. 1-11. 10.1186/1687-6180-2011-5.

[L,n] = size(X);
K = size(aoa,2);
omega = ang2dir(aoa,n);
A = exp(1i*pi*X*omega)  ;

% Optionally, compute the steering vector partial derivatives
if nargout > 1
    switch n
        case 1
            dA = A.*( 1i*pi*X * cos(aoa) );
            dA_s = A.*(1i*pi*X);
        case 2
            dA = A.*(1i*pi*X*[-sin(aoa(1,:)).*cos(aoa(2,:)); -cos(aoa(1,:)).*cos(aoa(2,:)) ] );
            dA = [dA, A.*(1i*pi*X*[cos(aoa(1,:)).*(-sin(aoa(2,:))); sin(aoa(1,:)).*(-sin(aoa(2,:))) ] ) ];
            dA_s =  [A.*(1i*pi*X(:,1)), A.*(1i*pi*X(:,2))];
        case 3
            dA = A.*(1i*pi*X*[-sin(aoa(1,:)).*cos(aoa(2,:)); -cos(aoa(1,:)).*cos(aoa(2,:)); sin(aoa(2,:)) ] );
            dA = [dA, A.*(1i*pi*X*[cos(aoa(1,:)).*(-sin(aoa(2,:))); sin(aoa(1,:)).*(-sin(aoa(2,:))); cos(aoa(2,:)) ] ) ];
            dA_s =  [A.*(1i*pi*X(:,1)), A.*(1i*pi*X(:,2)), A.*(1i*pi*X(:,3))];
    end
end
end

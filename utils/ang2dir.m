function [omega] = ang2dir(ang,n)
% ang in radians with size [n,Nang]
% For n=1, ang is a row vector of azimuth angles. For n=2,3, the first row is azimuth and the second elevation.
switch n
    case 1 %  omega in [-1,1]
        omega = sin(ang);
    case 2 % omega in {omega s.t. ||omega|| <= 1 }
        omega = [ cos(ang(1,:)).*cos(ang(2,:)); sin(ang(1,:)).*cos(ang(2,:)) ];
    case 3 % omega in {omega s.t. ||omega||=1 }
        omega = [ cos(ang(1,:)).*cos(ang(2,:)); sin(ang(1,:)).*cos(ang(2,:)); sin(ang(2,:)) ];
end
end

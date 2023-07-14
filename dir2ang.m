function [ang] = dir2ang(omega)
n = numel(omega);
switch n
    case 1 % ang(1,:) in [-pi/2,pi/2]
        ang = asin(omega);
    case 2 % ang(1,:) in [0,2*pi], ang(2,:) in [0,pi/2]
        ang = zeros([2,size(omega,2)]);
        ang(1,:) = wrapTo2Pi(atan2(omega(2,:), omega(1,:)));
        ang(2,:) = acos(omega(1,:)./cos(ang(1,:)));
    case 3 % ang(1,:) in [0,2*pi], ang(2,:) in [-pi/2,pi/2]
        ang = zeros([2,size(omega,2)]);
        ang(2) = asin(omega(3,:));
        ang(1) = omega(1,:)./cos(ang(2,:));
end
end
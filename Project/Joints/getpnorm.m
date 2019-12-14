function [phi,nu,gamma,dphi_dr,dphi_dp] = getpnorm(p,pdot)
% function to calculate the Euler param normalization constraint
phi = p(4)^2 + p(1)^2 + p(2)^2 + p(3)^2 - 1;
nu = 0;
gamma = -2*(pdot).'*pdot;
dphi_dr = zeros(1,3);
dphi_dp = [2*p(1), 2*p(2), 2*p(3), 2*p(4)];
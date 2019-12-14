clear all
close all
clc




L = 1;
m = 10;
g = 9.81;
k = 25;
phiFree = 0;
fun = @(phi)[98.1*cos(phi)+ 25*phi ]
phiGuess = [2 3]
phi = fsolve(fun,phiGuess);
%   phi = fsolve(@(phi) L*m*g*cos(phi) + k*phi, [0 180]);
  
  phi = radtodeg(phi)
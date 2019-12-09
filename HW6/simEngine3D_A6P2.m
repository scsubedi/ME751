%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 6, Problem 2

clear all
close all
clc


q = [-1 0 0
    0 0 0
    1 0 0];

t= 0;
r_i = [0 0 0]';
ei = [1 0 0]';
eiDot = [1 0 0]';
ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot=-ei'*eiDot/ei0;
p_i = [ei0, ei']';
p_iDot = [ei0Dot,eiDot']';

A_i = Amatrix(ei0, ei);
a_iBar = ([2 0 0] - [0 0 0])';
a_i =A_i* a_iBar;

ft = cos((pi*cos(2*t))/4 - pi/2);
fDott = ((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2);
fDotDott = (pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4);

r_j = [0 0 0]';
r_jDot = [0 0 0]';
ej = [0 0 0]';
ejDot = [0 0 0]';
ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
ej0Dot=-ej'*ejDot/ej0;
p_j = [ej0, ej']';
p_jDot = [ej0Dot,ejDot']';
A_j = Amatrix(ej0, ej);
a_jBar = [0 0 0]';
a_j =A_j* a_jBar;
s_jBar = a_jBar;

theta = 0;
c = [2*cos(theta),0 0]';
CD = c'*(r_j + A_j*a_j - r_i - A_i*a_i)- ft;
VelocityCD = fDott;
AccelerationCD = c'*Bmatrix(p_iDot,a_iBar)*p_iDot  - c'*Bmatrix(p_jDot,a_jBar)*p_jDot +fDotDott;
PartialPhi_ri= -c';
PartialPhi_rj = c';
PartialPhi_pi= -c'*Bmatrix(p_i,a_iBar);
PartialPhi_pj = c'*Bmatrix(p_j,a_jBar);
PartialPhi_rCD = [PartialPhi_ri PartialPhi_rj];
PartialPhi_pCD = [PartialPhi_pi PartialPhi_pj];
phi_SJ = [CD;CD;CD];
mu_SJ = [VelocityCD;VelocityCD;VelocityCD];
gamma_SJ = [AccelerationCD;AccelerationCD;AccelerationCD];

DP1 = -a_iBar'*A_i'*A_j*a_jBar - ft;
VelocityDP1 = fDott;
AccelerationDP1 = -a_iBar'*Bmatrix(p_jDot,a_jBar)*p_jDot ...
    - a_j'*Bmatrix(p_iDot,a_iBar)*p_iDot...
    - 2*((Bmatrix(p_i,a_iBar)*p_iDot)'*(Bmatrix(p_j,a_jBar)*p_jDot))...
    + fDotDott;
PartialPhi_ri = [zeros(1,3)];
PartialPhi_rj = [zeros(1,3)];
PartialPhi_pi = [a_j'*Bmatrix(p_i,a_iBar)];
PartialPhi_pj = [a_i'*Bmatrix(p_j,a_jBar)];
PartialPhi_rDP1 = [PartialPhi_ri PartialPhi_rj];
PartialPhi_pDP1 = [PartialPhi_pi PartialPhi_pj];


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling the constraint matrix

phiG_q = [PartialPhi_rCD PartialPhi_pCD
    PartialPhi_rCD PartialPhi_pCD
    PartialPhi_rCD PartialPhi_pCD
    PartialPhi_rDP1 PartialPhi_pDP1];
% phiQ_t = [phiQ_G;phiQt_D]

phiP_p = (p_i'*p_i*ones(7,1) - ones(7,1))';
phiP_p = [phiP_p zeros(1,7)];
phiD_qt = [0   1.4142   -1.4142    0.6533    0.2706    0.6533    0.2706];
phiD_qt = [phiD_qt zeros(1,7)];
phi_qt = [phiG_q;phiD_qt]
phi_qtF = [phiG_q;phiD_qt; phiP_p]


Phi = [phi_SJ;DP1]
mu = [mu_SJ;VelocityDP1]
gamma = [gamma_SJ;AccelerationDP1]



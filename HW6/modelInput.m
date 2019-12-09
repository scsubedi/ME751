%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this file to define all the inputs and also to check the constraints
% The codes for the earlier defined constraints have also been updated
% The constraints have now been defined into a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


syms t;
c = [0.3,0.4,-6]';
r_i = [8, 6,-3]';
ei = [0.4201 -0.7001 0.1400]';

eiDot = [0.3566 0.9326 0]';
ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot=-ei'*eiDot/ei0;
p_i = [ei0, ei']';
p_iDot = [ei0Dot,eiDot']';

A_i = Amatrix(ei0, ei);
a_iBar = [-1.2, 1 ,0.3]';
a_i =A_i* a_iBar;

s_iBar = a_iBar;

f = pi()/4*cos(2*t);
df = diff(f,t);
dff = diff(df,t);

t0 = 2;
ft=vpa(subs(f,t,t0));
fDott=vpa(subs(df,t,t0));
fDotDott = vpa(subs(dff,t,t0));

r_j = [-0.5,1.6,-6.3]';
ej = [0.2 0.2 0.2]';
ejDot = [1 2 9]';
ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
ej0Dot=-ej'*ejDot/ej0;
p_j = [ej0, ej']';
p_jDot = [ej0Dot,ejDot']';
A_j = Amatrix(ej0, ej);
a_jBar = [3 5 2]';
a_j =A_j* a_jBar;

r_iDot = [7,8,9]';
r_jDot =[11, 12, 13]';

s_jBar = a_jBar;
dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
dijDot =  r_jDot - r_iDot + Bmatrix(p_j,s_jBar)*p_j - Bmatrix(p_i,s_iBar)*p_i;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs is an array with all the details about the bodies with their
% constraints defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
    A_i;a_iBar;a_i;s_iBar;f;df;dff;t0;ft;fDott;...
    fDotDott;r_j;ej;ejDot;ej0;ej0Dot;p_j;p_jDot;...
    A_j;a_jBar;a_j;s_jBar;dij;dijDot;r_iDot;r_jDot};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks return the specified output requested
% flag = 1 : Returns the value of the expression of the constraint
% flag = 2 : Returns the right-hand side of the velocity equation
% flag = 3 : Returns the right-hand side of the acceleration equation
% flag = 4 : Returns the expression of partial derivatives
% flag = 5 : Returns the expression of partial derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GroundStatus
% Grounded : 1
% UnGrounded: 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groundStatus = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying appropriate contraint to the body

% [CD, VelocityCD, AccelerationCD, PartialPhi_rCD, PartialPhi_pCD] = cons_cd(inputs,flag,groundStatus);

% [ DP1, VelocityDP1, AccelerationDP1, PartialPhi_rDP1, PartialPhi_pDP1] = cons_dp1(inputs,flag,groundStatus);

[ D, VelocityD, AccelerationD, PartialPhi_rD, PartialPhi_pD] = cons_d(inputs,flag,groundStatus);

[ DP2, VelocityDP2, AccelerationDP2, PartialPhi_rDP2, PartialPhi_pDP2] = cons_dp2(inputs,flag,groundStatus);


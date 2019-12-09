
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all the inputs for the constraint

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

[ D, VelocityD, AccelerationD, PartialPhi_rD, PartialPhi_pD] = cons_d(inputs,flag,groundStatus);

% [ DP1, VelocityDP1, AccelerationDP1, PartialPhi_rDP1, PartialPhi_pDP1] = cons_dp1(inputs,flag,groundStatus);

[ DP2, VelocityDP2, AccelerationDP2, PartialPhi_rDP2, PartialPhi_pDP2] = cons_dp2(inputs,flag,groundStatus);

%% Assignment 6 Part 2
% Defining the body 1
clear all
close all
clc;
body1 = [0 0 0
    0 1 0
    0 0 -1];

%defining the body 2
body2 = [0 0 0
    0 0 1
    -2 0 0
    1 0 0];
bodies = {body1;body2};
%no. of constraints
nc = 1;

%no. of bodies
nb = length(bodies)-1;

q = [0;0;0;1;0;0;0;0];
qDot = q;
nc = 1; % no. of joints
nb = 2-1;  %no. of bodies -1

%kinematic constraints
syms theta % t
% theta= 0:0.1: 2*pi;
% SphericalJoint = []
c = [2*cos(theta),0 0]';
r_i = [0 0 0]';
ei = [1 0 0]';

eiDot = [1 0 0]';
ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot=-ei'*eiDot/ei0;
p_i = [ei0, ei']';
p_iDot = [ei0Dot,eiDot']';
A_i = Amatrix(ei0, ei);


A_i = Amatrix(ei0, ei);
a_iBar = [-2 0 0]';
a_i =A_i* a_iBar;


r_j = [0 0 0]';
ej = [0 0 0]';
ejDot = [0 0 0]';
ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
ej0Dot=-ej'*ejDot/ej0;
p_j = [ej0, ej']';
p_jDot = [ej0Dot,ejDot']';
A_j = Amatrix(ej0, ej);

a_jBar = [0 0 0]';
a_j =A_j* a_jBar;
t = 0;
ft = pi()/4*cos(2*t);
CD = c'*(r_j + A_j*a_j - r_i - A_i*a_i)- ft;



%     VelocityCD = fDott;
%     AccelerationCD = c'*Bmatrix(p_iDot,a_iBar)*p_iDot + fDotDott;
%     PartialPhi_ri= -c';
%     PartialPhi_pi= -c'*Bmatrix(p_i,a_iBar);
%     PartialPhi_rCD = [PartialPhi_ri];
%     PartialPhi_pCD = [PartialPhi_pi];


%%

constraint_q(q,inputs,nc,nb);

kinematicConstraint(q,inputs,nb,nc);

velocity(inputs,q,nb,nc);

acceleration(inputs,q,qDot,nb,nc);

function phi_q = constraint_q(q,inputs,nc,nb)
q = reshape(q,length(q),1);
% phi_q is a matrix of nc constraints + euler constraints x
% variables which are nb*7
phi_q = zeros(nc + nb,7*nb);
for ii = 1:nc
    phi_q(ii,:) = cons_dp1(inputs,1,1);
end
% adding the euler parameterization constraints
% skipping first body since it is the ground body
for ii = 1 : nb
    indexing = nb*3 + ii*4 - 3 : nb*3 + ii*4;
    phi_q(nc + ii,indexing) = 2*q(indexing,1)';
end
end


function muF = velocity(inputs,q,nb,nc)
q = reshape(q,length(q),1);
%             obj = obj.setq(q);
muF = zeros(nc + nb - 1,1);
for ii = 1:nc
    muF(ii) = cons_dp1(inputs,2,1);
end
% adding the euler parameterization constraints
% skipping first body since it is the ground body
for ii = 2 : nb
    muF(nc + ii - 1) = 0;
end
end

function gammaF = acceleration(inputs,q,qDot,nb,nc)
q = reshape(q,length(q),1);
qDot = reshape(qDot,length(qDot),1);

gammaF = zeros(nc + nb - 1,1);
for ii = 1:nc
    gammaF(ii) = cons_dp1(inputs,3,1);
end
% adding the euler parameterization constraints
% skipping first body since it is the ground body
for ii = 2 : nb
    gammaF(nc + ii - 1) = -2*(p_iDot'*p_iDot);
end
end

function phiF = kinematicConstraint(q,inputs,nb,nc)
% This function computes all the kinematics constraints
% input - a vector of [q t]
% we need to remove the ground bodies euler parameterization
% constraint therefore : -1
q = reshape(q,length(q),1);
%
phiF = zeros(nc + nb - 1,1);
for i = 1:nc
    phiF(i) = cons_dp1(inputs,0,1);
end
% adding the euler parameterization constraints
% skipping first body since it is the ground body
for i = 2 : nb
    phiF(nc + i - 1) ...
        = p_i'*p_i - 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this file to define all the inputs and also to check the constraints
% The codes for the earlier defined constraints have also been updated
% The constraints have now been defined into a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


readInput;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs is an array with all the details about the bodies with their
% constraints defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
    A_i;a_iBar;a_i;s_iBar;f;df;dff;time;ft;fDott;...
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
% flag = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GroundStatus
% Grounded : 1
% UnGrounded: 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% groundStatus = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying appropriate contraint to the body

if Constraint ==1
    
    [CD, VelocityCD, AccelerationCD, PartialPhi_rCD, PartialPhi_pCD] = cons_cd(inputs,flag,groundStatus);
    
elseif Constraint == 2
    [ DP1, VelocityDP1, AccelerationDP1, PartialPhi_rDP1, PartialPhi_pDP1] = cons_dp1(inputs,flag,groundStatus);
elseif Constraint == 3
    [ DP2, VelocityDP2, AccelerationDP2, PartialPhi_rDP2, PartialPhi_pDP2] = cons_dp2(inputs,flag,groundStatus);
elseif Constraint ==4
    [ D, VelocityD, AccelerationD, PartialPhi_rD, PartialPhi_pD] = cons_d(inputs,flag,groundStatus);
    
else
    disp('Not a defined constraint')
end
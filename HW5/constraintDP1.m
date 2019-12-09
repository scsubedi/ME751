function [ DP1, VelocityDP1, AccelerationDP1, PartialPhi_rDP1, PartialPhi_pDP1] = cons_d(inputs,flags,groundStatus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric Constraint: DP1 Constraint
% Input:
% Manual input for the bodies Euler Parameters of body I and J
% Switch for Grounding body J

% ?DP1(i , a ? i , j , a ? j , f(t)) = a ? i T A i T A j a ? j - f(t) = 0
% Checks return the specified output requested
% check = 0 : Returns the value of the expression of the constraint
% check = 1 : Returns the right-hand side of the velocity equation
% check = 2: Returns the right-hand side of the acceleration equation
% check = 3: Returns the expression of partial derivatives
% check = 4:   Returns the expression of partial derivatives

% Output:
% Display of the quantities of interest
% Saves a file with all the values computed as Results.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
syms t;
str1  = 'DP1';
str2 = 'unGrounded';
str3 = 'Grounded';

% Body J can be grounded
bodyJGround = false;

% Change the value to change the output as indicated above
check = 4;

r_i = [8 6 3]';
ei = [0.4201 -0.7001 0.1400]';

eiDot = [0.3566 0.9326 0]';
ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot=-ei'*eiDot/ei0;
p_i = [ei0, ei']';
p_iDot = [ei0Dot,eiDot']';

A_i = Amatrix(ei0, ei);
a_iBar = [-1.2, 1 ,0.3]';
a_i =A_i* a_iBar;

f = t^2;
df = diff(f,t);
dff = diff(df,t);

t0 = 2;
ft=vpa(subs(f,t,t0));
fDott=vpa(subs(df,t,t0));
fDotDott = vpa(subs(dff,t,t0));


if(bodyJGround == false)
    
    r_j = [1 2 3]';
    ej = [0.2 0.2 0.2]';
    ejDot = [1 2 9]';
    ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
    ej0Dot=-ej'*ejDot/ej0;
    p_j = [ej0, ej']';
    p_jDot = [ej0Dot,ejDot']';
    
    A_j = Amatrix(ej0, ej);
    a_jBar = [3 5 2]';
    
    a_j =A_j* a_jBar;
    
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
    
    ResultsUnGroundedDP1 = {DP1 VelocityDP1...
        AccelerationDP1 PartialPhi_ri PartialPhi_rj ...
        PartialPhi_pi PartialPhi_pj};
    saveFile(ResultsUnGroundedDP1,str1,str2);
    
else
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
    
    DP1 = -a_iBar'*A_i'*A_j*a_jBar - ft;
    VelocityDP1 = fDott;
    AccelerationDP1 = - a_j'*Bmatrix(p_iDot,a_iBar)*p_iDot  + fDotDott;
    PartialPhi_ri = [zeros(1,3) zeros(1,3)];
    PartialPhi_pi = [a_j'*Bmatrix(p_i,a_iBar)];
    
    ResultsGroundedDP1 = {DP1 VelocityDP1 AccelerationDP1...
        PartialPhi_ri PartialPhi_pi};
    
    saveFile(ResultsGroundedDP1,str1,str3);
end


%%%%% Displaying Results based on the use input and choice of output %%%%%
if (bodyJGround== true)
    fprintf('The body J is grounded. \n');
    
else
    fprintf('The bodies I and J have a Co-ordinate Difference(CD) constraint. \n');
end
 fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
if (check == 1)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The value of the expression of the constraint is %f \n', DP1);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end
 
if (check == 2)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the velocity equation(mu) is %.5f \n', VelocityDP1);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end
 
if (check == 3)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the acceleration equation(gamma) is %.5f \n', AccelerationDP1);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end
 
if (check == 4)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The expression for partial derivatives phi_R \n');
    if (bodyJGround== true)
        fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
        PartialPhi_ri
        fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    else
        fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
        PartialPhi_ri
        PartialPhi_rj
        fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    end
end
if (check == 5)
    
    fprintf('The expression for partial derivatives phi_p ');
    if (bodyJGround== true)
        
        PartialPhi_pi
    else
        PartialPhi_pi
        PartialPhi_pj
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric Constraint 1: CD Constraint
% Input: 
% Manual input for the bodies Euler Parameters of body I and J
% Switch for Grounding body J

% ?CD(c,i,sip_,j,sjq_,f(t)) =c'dji - f(t) = 0
% Checks return the specified output requested
% check = 0 : Returns the value of the expression of the constraint
% check = 1 : Returns the right-hand side of the velocity equation
% check = 2 : Returns the right-hand side of the acceleration equation
% check = 3 : Returns the expression of partial derivatives
% check = 4 : Returns the expression of partial derivatives

% Output: 
% Display of the quantities of interest
% Saves a file with all the values computed as Results.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

syms t;
str1  = 'CD';
str2 = 'unGrounded';
str3 = 'Grounded';

%Body J can be grounded
bodyJGround = true;

% Change the value to change the output as indicated above
check = 3;
format short

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

f = t^2;
df = diff(f,t);
dff = diff(df,t);

t0 = 2;
ft=vpa(subs(f,t,t0));
fDott=vpa(subs(df,t,t0));
fDotDott = vpa(subs(dff,t,t0));

if (bodyJGround == false)
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
    
    s_jBar = a_jBar;

    CD = c'*(r_j + A_j*a_j - r_i - A_i*a_i)- ft;
    VelocityCD = fDott;
    AccelerationCD = c'*Bmatrix(p_iDot,a_iBar)*p_iDot  - c'*Bmatrix(p_jDot,a_jBar)*p_jDot +fDotDott;
    PartialPhi_ri= -c';
    PartialPhi_rj = c';
    PartialPhi_pi= -c'*Bmatrix(p_i,a_iBar);
    PartialPhi_pj = c'*Bmatrix(p_j,a_jBar);
    ResultsUnGroundedCD = {CD VelocityCD AccelerationCD PartialPhi_ri PartialPhi_rj PartialPhi_pi PartialPhi_pj};
    saveFile(ResultsUnGroundedCD,str1,str2);
    
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
    
    CD = c'*(r_j + A_j*a_j - r_i - A_i*a_i)- ft;
    VelocityCD = fDott;
    AccelerationCD = c'*Bmatrix(p_iDot,a_iBar)*p_iDot + fDotDott;
    PartialPhi_ri= -c';
    PartialPhi_pi= -c'*Bmatrix(p_i,a_iBar);
       
    ResultsGroundedCD = {CD VelocityCD AccelerationCD PartialPhi_ri PartialPhi_pi};
        
    saveFile(ResultsGroundedCD,str1,str3);
    
end

%%%%% Displaying Results based on the use input and choice of output %%%%%

if (bodyJGround== true)
    fprintf('The body J is grounded. \n');
    
else
    fprintf('The bodies I and J have a Co-ordinate Difference(CD) constraint. \n');
end

if (check == 1)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The value of the expression of the constraint is %f \n', CD);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end

if (check == 2)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the velocity equation(mu) is %.5f \n', VelocityCD);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end

if (check == 3)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the acceleration equation(gamma) is %.5f \n', AccelerationCD);
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

fileID = fopen('Results.txt','r');



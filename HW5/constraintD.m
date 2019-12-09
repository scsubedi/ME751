function A = cons_d(inputs,flags,groundStatus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric Constraint: Distance Constraint
% Input:
% Manual input for the bodies Euler Parameters of body I and J
% Switch for Grounding body J

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

syms t;
str1  = 'CD';

if (groundStatus == 0)
    str2 = 'Ungrounded';
    bodyJGround = false;
elseif (groundStatus == 1)
    str3 = 'Grounded';
    bodyJGround = true;
else
    disp('Invalid Ground status');
end

% Change the value to change the output as indicated above
check = flags;
format short

% % Body J can be grounded
% bodyJGround = true;
%
% % Change the value to change the output as indicated above
% check = 5;
%
% r_i = [8 6 3]';
% r_iDot = [7,8,9]';
% ei = [0.4201 -0.7001 0.1400]';
%
% eiDot = [0.3566 0.9326 0]';
% ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
% ei0Dot=-ei'*eiDot/ei0;
% p_i = [ei0, ei']';
% p_iDot = [ei0Dot,eiDot']';
%
% A_i = Amatrix(ei0, ei);
% a_iBar = [-1.2, 1 ,0.3]';
% a_i =A_i* a_iBar;
%
% s_iBar = a_iBar;
% a_iDot = Bmatrix(p_i,a_iBar)* p_iDot;
%
%
% f = t^2;
% df = diff(f,t);
% dff = diff(df,t);
%
% t0 = 2;
% ft=vpa(subs(f,t,t0));
% fDott=vpa(subs(df,t,t0));
% fDotDott = vpa(subs(dff,t,t0));

c = inputs{1};
r_i = inputs{2};
ei = inputs{3};
eiDot = inputs{4};
ei0 = inputs{5};
ei0Dot = inputs{6};
p_i = inputs{7};
p_iDot = inputs{8};
A_i = inputs{9};
a_iBar = inputs{10};
a_i = inputs{11};
s_iBar= inputs{12};
f = inputs{13};
df = inputs{14};
dff = inputs{15};
t0 = inputs{16};
ft = inputs{17};
fDott = inputs{18};
fDotDott = inputs{19};

if(bodyJGround == false)
    
    r_j = [1 2 3]';
    r_jDot = [11, 12, 13]';
    
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
    
    dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
    
    dijDot =  r_jDot - r_iDot + Bmatrix(p_j,s_jBar)*p_j - Bmatrix(p_i,s_iBar)*p_i;
    
    
    r_j = inputs{20};
    ej = inputs{21};
    ejDot = inputs{22};
    ej0 = inputs{23};
    ej0Dot= inputs{24};
    p_j = inputs{25};
    p_jDot = inputs{26};
    A_j = inputs{27};
    a_jBar = inputs{28};
    a_j = inputs{29};
    s_jBar = inputs{30};
    dij = inputs{31};
    dijDot = inputs{32};
    
    
    D = dij'*dij - ft;
    
    VelocityD = fDott;
    AccelerationD = -2*dij'*Bmatrix(p_jDot,s_jBar)*p_jDot...
        + 2*dij'*Bmatrix(p_iDot,s_iBar)*p_iDot - 2*(dijDot')*dijDot + fDotDott;
    
    PartialPhi_ri = -2*dij';
    PartialPhi_rj = 2*dij';
    PartialPhi_pi = -2*(dij')*Bmatrix(p_i,s_iBar);
    PartialPhi_pj = -2*dij'*Bmatrix(p_j,s_jBar);
    
    ResultsUnGroundedD = {D VelocityD...
        AccelerationD PartialPhi_ri PartialPhi_rj ...
        PartialPhi_pi PartialPhi_pj};
    saveFile(ResultsUnGroundedD,str1,str2);
    
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
    
    D = -a_iBar'*A_i'*A_j*a_jBar - ft;
    VelocityD = fDott;
    AccelerationD = - a_j'*Bmatrix(p_iDot,a_iBar)*p_iDot  + fDotDott;
    PartialPhi_ri = [zeros(1,3)];
    PartialPhi_pi = [zeros(1,4)];
    
    PartialPhi_rD = [PartialPhi_ri];
    PartialPhi_pD = [PartialPhi_pi];
    
    
    ResultsGroundedD = {D VelocityD AccelerationD...
        PartialPhi_ri PartialPhi_pi};
    
    saveFile(ResultsGroundedD,str1,str3);
end


%%%%% Displaying Results based on the use input and choice of output %%%%%
if (bodyJGround== true)
    fprintf('The body J is grounded. \n');
    
else
    fprintf('The bodies I and J have a Distance(D) constraint. \n');
end
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
if (check == 1)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The value of the expression of the constraint is %f \n', D);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end

if (check == 2)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the velocity equation(mu) is %.5f \n', VelocityD);
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end

if (check == 3)
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
    fprintf('The right-hand side of the acceleration equation(gamma) is %.5f \n', AccelerationD);
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

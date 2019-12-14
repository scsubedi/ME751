clear all
close all
clc

% syms t;
% c = [0.3,0.4,-6]';
% r_i = [8, 6,-3]';
% ei = [0.4201 -0.7001 0.1400]';
%
% eiDot = [0.3566 0.9326 0]';
% ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
% ei0Dot=-ei'*eiDot/ei0;
% p_i = [ei0, ei']';
% p_iDot = [ei0Dot,eiDot']';
%
% A_i = Amatrix(p_i);
% a_iBar = [-1.2, 1 ,0.3]';
% a_i =A_i* a_iBar;
%
% s_iBar = a_iBar;
%
% f = pi()/4*cos(2*t);
% df = diff(f,t);
% dff = diff(df,t);
%
% t0 = 2;
% ft=vpa(subs(f,t,t0));
% fDott=vpa(subs(df,t,t0));
% fDotDott = vpa(subs(dff,t,t0));
%
% r_j = [-0.5,1.6,-6.3]';
% ej = [0.2 0.2 0.2]';
% ejDot = [1 2 9]';
% ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
% ej0Dot=-ej'*ejDot/ej0;
% p_j = [ej0, ej']';
% p_jDot = [ej0Dot,ejDot']';
% A_j = Amatrix(p_j);
% a_jBar = [3 5 2]';
% a_j =A_j* a_jBar;
%
% r_iDot = [7,8,9]';
% r_jDot =[11, 12, 13]';
%
% s_jBar = a_jBar;
% dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
% dijDot =  r_jDot - r_iDot + Bmatrix(p_j,s_jBar)*p_j - Bmatrix(p_i,s_iBar)*p_i;

% readInput;
syms t;
str = fileread('inputTest.txt');
val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
val = strtrim(vertcat(val{:}));
num = cellfun(@(s)sscanf(s,'%f:'),val(:,2),'UniformOutput',false);
len = cellfun('length',num);
vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);
% val(len>0,2) = num(len>0);
% out = cell2struct(val(:,2),regexprep(val(:,1),' ','_'));

c = str2num(val{3,2})';
r_i = str2num(val{4,2})';
r_iDot = str2num(val{5,2})';
ei = str2num(val{6,2})';
eiDot = str2num(val{7,2})';
a_iBar = str2num(val{8,2})';

f = str2func(num2str(val{9,2}));

% time =str2double(val{10,2});
time = 2;
r_j = str2num(val{11,2})';
r_jDot = str2num(val{12,2})';
ej = str2num(val{13,2})';
ejDot = str2num(val{14,2})';
a_jBar = str2num(val{15,2})';



ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot=-ei'*eiDot/ei0;
p_i = [ei0, ei']';
p_iDot = [ei0Dot,eiDot']';
A_i = Amatrix(p_i);
a_i =A_i* a_iBar;
s_iBar = a_iBar;

ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
ej0Dot=-ej'*ejDot/ej0;
p_j = [ej0, ej']';
p_jDot = [ej0Dot,ejDot']';
A_j = Amatrix(p_j);
a_j =A_j* a_jBar;
s_jBar = a_jBar;
df = diff(f,t);
dff = diff(df,t);
ft=f(t);
fDott=vpa(subs(df,t,time));
fDotDott = vpa(subs(dff,t,time));
dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
dijDot =  r_jDot - r_iDot + Bmatrix(p_j,s_jBar)*p_j - Bmatrix(p_i,s_iBar)*p_i;




%%

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
flag = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GroundStatus
% Grounded : 1
% UnGrounded: 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groundStatus = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying appropriate contraint to the body
%% Spherical Joint

clear all
close all
clc
flag = 1;
groundStatus = 0;


str1  = 'Spherical Joint';
str4 = 'Results_SJ';

if (groundStatus == 0)
    str2 = 'Ungrounded';
    bodyJGround = false;
elseif (groundStatus == 1)
    str3 = 'Grounded';
    bodyJGround = true;
else
    disp('Invalid Ground status');
end
readInput;

 for w = 1:3
        c = zeros(3,1);
        c(w) = 1;
%         gcon_inp.c = c;
%         [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getcd(gcon_inp);




CD = c.'*(r_j + A_j*a_j - r_i - A_i*a_i)- ft;
VelocityCD = fDott;
AccelerationCD = c.'*Bmatrix(p_iDot,a_iBar)*p_iDot  - c.'*Bmatrix(p_jDot,a_jBar)*p_jDot +fDotDott;
PartialPhi_ri= -c.'*eye(3);
PartialPhi_rj = c.'*eye(3);
PartialPhi_pi= -c.'*Bmatrix(p_i,a_iBar);
PartialPhi_pj = c.'*Bmatrix(p_j,a_jBar);
PartialPhi_rCD = [PartialPhi_ri PartialPhi_rj];
PartialPhi_pCD = [PartialPhi_pi PartialPhi_pj];

ResultsUnGroundedCD{w} = {CD(w) VelocityCD(w) AccelerationCD(w) PartialPhi_ri(w) PartialPhi_rj(w) PartialPhi_pi(w) PartialPhi_pj(w)};
 end

% [CD, VelocityCD, AccelerationCD, PartialPhi_rCD, PartialPhi_pCD] = cons_cd(flag,groundStatus);
%%
for w = 1:3
    c = zeros(3,1);
    c(w) = 1;
    gcon_inp.c = c;
    [CD(w,1),VelocityCD(w,1),AccelerationCD(w,1),PartialPhi_rCD(w,:),PartialPhi_pCD(w,:)] = cons_cd(inputs,flag,groundStatus);
end
ResultsSJ= {CD VelocityCD AccelerationCD PartialPhi_rCD PartialPhi_pCD};
saveFile(ResultsSJ,str1,str2,str4);


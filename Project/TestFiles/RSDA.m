% TSDA

% Defining the bodies
% Setting up translational spring-damper-actuator acting between point Pi
% on body i,and Pj on body j
clear all
close all
clc

j = 1;
r_i = [3 2 0]';
r_iDot = [0 0 0]';
s_iBar = [0 4 0]';
ei = [0 1 0];
eiDot = [0 0 0];
ei0 = sqrt(ei(1)^2 + ei(2)^2 + ei(3)^2);
ei0Dot = [0];
p_i = [ei0; ei'];
p_iDot = [ei0Dot; eiDot'];
A_i = Amatrix(p_i);
a_iBar = s_iBar;
a_i = A_i*a_iBar;
if j == 0
    r_j = [0 0 0]';
    r_jDot = [0 0 0]';
    s_jBar= [0 0 0]';
    ej = [0 0 0]';
    ejDot = [0 0 0]';
    
elseif j == 1
    
    r_j = [3 5 0]';
    r_jDot = [0 0 0]';
    s_jBar= [0 3 0]';
    ej = [0 1 0]';
    ejDot = [0 0 0]';
else
    disp('Unrecognized option, check the status of body j.')
end

ej0 = sqrt(ej(1)^2 + ej(2)^2 + ej(3)^2);
ej0Dot = 0;
q_j = [ej0;ej];

a_jBar = s_jBar;
q_jDot = [ej0Dot; ejDot];
A_j = Amatrix(q_j);
a_j = A_j*a_jBar;

dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
dijDot =  r_jDot - r_iDot + Bmatrix(q_j,s_jBar)*q_j - Bmatrix(p_i,s_iBar)*p_i;


lij = norm(dij);

if lij ==0
    lij = lij + 0.01;
end

eij = dij./lij;

lijDot = eij.'*dijDot;


%%% new variables defined;
l0 = 0.75;  % Zero stress length of spring
k = 25;     % Spring stiffness
h = 10;     % actuation force as a function of (l,ldot, time)
c = 50;     % Damping coefficient
syms t f

time = 0:0.1:10;
for i=1:length(time)
    h = lij + lijDot +(2*sin(2*pi*time(i)));
    H(i) = h;
    f = real(vpa(k*(lij - l0) + c*lijDot' + h'));
    
    fTSDA = eij.*f;
    Ftsda{i,:} =fTSDA;
end


%%


% readInput;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs is an array with all the details about the bodies with their
% constraints defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
%     A_i;a_iBar;a_i;s_iBar;f;df;dff;time;ft;fDott;...
%     fDotDott;r_j;ej;ejDot;ej0;ej0Dot;p_j;p_jDot;...
%     A_j;a_jBar;a_j;s_jBar;dij;dijDot;r_iDot;r_jDot};

%
% ft = cos(theta(t) + pi/2);                       % function for driving constraint, and its derivatives
% df
syms t

ts = 0:0.01:1;
for n = 1:length(ts)
    t = ts(n);
    syms t
    f = pi()/4*cos(2*t);
    df = diff(f,t);
    dff = diff(df,t);
    
    t0 = 2;
    ft=vpa(subs(f,t,t0));
    fDott=vpa(subs(df,t,t0));
    fDotDott = vpa(subs(dff,t,t0));
    
    flags = 1;
    groundStatus = 0;
    for w = 1:3
        c = zeros(3,1);
        
        c(w) = 1;
        inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
            A_i;a_iBar;a_i;s_iBar;f;df;dff;time;ft;fDott;...
            fDotDott;r_j;ej;ejDot;ej0;ej0Dot;q_j;q_jDot;...
            A_j;a_jBar;a_j;s_jBar;dij;dijDot;r_iDot;r_jDot};
        [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_cd(inputs,flags,groundStatus);
    end
    
    % calling 2 DP1 GCon to add to SJ to form Revolute joint
    % defining such that x'-y' plane is perpendicular to X-axis
    a_j = [1,0,0].';                                % x-axis in G-RF
    a_i = [1,0,0].';                                % x-axis in L-RF
    b_i = [0,1,0].';                                % y-axis in L-RF
    
    % gcon_inp.a_2 = a_2;
    % gcon_inp.a_1 = a_1;
    inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
        A_i;a_iBar;a_i;s_iBar;f;df;dff;time;ft;fDott;...
        fDotDott;r_j;ej;ejDot;ej0;ej0Dot;q_j;q_jDot;...
        A_j;a_jBar;a_j;s_jBar;dij;dijDot;r_iDot;r_jDot};
    
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_dp1(inputs,flags,groundStatus);
    a_i = b_i;
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_dp1(inputs,flags,groundStatus);
    
    a_i = [0,1,0].';
    a_j = [0,0,-1].';
    
    
    inputs = {c;r_i;ei;eiDot;ei0;ei0Dot;p_i;p_iDot;...
        A_i;a_iBar;a_i;s_iBar;f;df;dff;time;ft;fDott;...
        fDotDott;r_j;ej;ejDot;ej0;ej0Dot;q_j;q_jDot;...
        A_j;a_jBar;a_j;s_jBar;dij;dijDot;r_iDot;r_jDot};
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_dp1(inputs,flags,groundStatus);
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getpnorm(p_i, p_iDot);
    
    Phi_q = [dphi_dr,dphi_dp];
    q(:,n) = [r_i;p_i];
    qbar(:,n) = A_i\q(1:3,n);      % in L-RF
    qdot(:,n) = Phi_q.\nu;
    q_dot(:,n) = A_i\qdot(1:3,n);   % in L-RF
    qddot(:,n) = Phi_q\gamma;
    q_ddot(:,n) = A1(theta0)\qddot(1:3,n);% in L-RF
end
% end
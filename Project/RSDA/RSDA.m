%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position, Velocity and Acceleration analysis of a single Pendulum with Revolute
% Joint with ground. The motion is simulated along with its position,
% velocity and acceleration. There's an option to save the video of the
% motion.
% Input: input all the attributes in the 'input.txt' file
% Output: Plot for position, velocity and acceleration, option to save
% video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace

clear;clc;
close all;

% disp('The two bodies are connected by a Revolute Joint.')
% Reading the input file by avoiding white spaces and delimiters
% syms t
% filename  ='input';
% str = fileread([filename '.txt']);
% val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
% val = strtrim(vertcat(val{:}));
% num = cellfun(@(s)sscanf(s,'%f:'),val(:,2),'UniformOutput',false);
% len = cellfun('length',num);
% vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
% num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);
% 
% i = str2num(val{3,2})';
% j = str2num(val{4,2})';
% L = str2num(val{5,2})';
% s1_p = L*str2num(val{6,2})';
% s2_q = L*str2num(val{7,2})';
% 
% ts = str2num(val{8,2});
% 
% joint = val{10,2};
% saveResults = val{11,2};
readInput;
A1 = @(th)[0,  0,      1;                        % Rotation matrix as a function of angle
    sin(th),  cos(th), 0;
    -cos(th), sin(th), 0];

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);

k = 25;
c = 10;

% RSDA Implementation
% variables derived from above input
for n = 1:length(ts)
    t = ts(n);
    omega = [dtheta(t),0,0].';                      % angular velocity obtained from f(t) = Acos(omega*t)
    theta0 = theta(t);
    r1 = A1(theta0)*[L,0,0].';
    r2 = zeros(3,1);
    d12 = r1- r2;
    thetaDot = theta0/t;
    h = theta;
    torque = k*(theta0-theta(1)) + c* thetaDot + (theta0 + thetaDot*sin(2*pi*t));
    
    %     torque = k*(theta0-theta(1)) + c* thetaDot - exp(t);
    % Torque acting on body 1 is given by
    torque1 = torque.*r1;
    torque2 = torque.*r2;
    nRSDA_1{n,1} =torque1(1);
    nRSDA_1{n,2} = torque1(2);
    nRSDA_1{n,3} = torque1(3);
    rsda = sqrt((cell2mat(nRSDA_1(:,2))).^2 + (cell2mat(nRSDA_1(:,2)).^2));
    
    
end
%%
rsda = sqrt((cell2mat(nRSDA_1(:,2))).^2 + (cell2mat(nRSDA_1(:,2)).^2));
figure(1)
plot(ts,rsda)
grid on
title('Concentrated Torque due to Rotational Spring-Damper-Actuator')
xlabel('Time(s)')
ylabel('Torque')

% output = {Phi_q, q,qbar, qdot q_dot qddot q_ddot ts};
%%
%
% subplot(2,3,1)
% plot(ts,q(2,:))
% subplot(2,3,2)
% plot(ts,qdot(2,:)) ;
% subplot(2,3,3)
% plot(ts,qddot(2,:))
% subplot(2,3,4)
% plot(ts,q(3,:)) ;
% subplot(2,3,5)
% plot(ts,qdot(3,:))
% subplot(2,3,6)
% plot(ts,qddot(3,:))

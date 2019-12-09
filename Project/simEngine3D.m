% ME751 A6P3 : pendulum, kinematic analysis, ts = 0:dt:10
% in this problem, body i/1 is the bar and body j/2 is ground
clear;clc;
close all;

%%
str = fileread('input.txt');
val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
val = strtrim(vertcat(val{:}));
num = cellfun(@(s)sscanf(s,'%d:'),val(:,2),'UniformOutput',false);
len = cellfun('length',num);
vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);
val(len>0,2) = num(len>0);
out = cell2struct(val(:,2),regexprep(val(:,1),' ','_'));


%%
% defining input variables
i = val{1,2};
j = val{2,2};                                  % setting second body to ground
L = val{3,2};                                  % half-span of pendulum


%%
% i = 1;
% j = 0;
% L = 2;
A1 = @(th)[0,        0,      1;                % A as a function of theta
    sin(th),  cos(th), 0;
    -cos(th), sin(th), 0];

% theta = @(t)val{4,2};

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);


starttime = val{9,2};
dt = 10^(val{8,2});
finalTime = val{10,2};
ts = starttime:dt:finalTime;
% dt = 1e-2;
% ts = 0:dt:10;
s1_p = [-L,0,0].';                              % defining P and Q for CD GCon
s2_q = [0,0,0].';


%%
% variables derived from above input
for n = 1:length(ts)
    t = ts(n);
omega = [dtheta(t),0,0].';                      % angular velocity obtained from f(t) = Acos(omega*t)
theta0 = theta(t);
r1 = A1(theta0)*[L,0,0].';
r2 = zeros(3,1);
ft = cos(theta(t) + pi/2);                       % function for driving constraint, and its derivatives
dft =-sin(theta(t) + pi/2)*dtheta(t);
ddft = -sin(theta(t) + pi/2)*ddtheta(t) - cos(theta(t) + pi/2)*dtheta(t)^2;


% calculating p and pdot from input
p1 = getp(A1(theta0));                           % calculate p from A
    unitycheck = p1(1)^2 + p1(2:4).'*p1(2:4);    % verifying that the A.'*A=eye(3) condition is met
    e1 = p1(2:4);    
    E1 = [-e1, tilda(e1)+p1(1)*eye(3)];
p1dot = 0.5*E1.'*omega;                          % using omega and E to get pdot
    orthocheck = p1dot.'*p1;
p2 = getp(eye(3));                               % euler paramters for G-RF
p2dot = zeros(4,1);

% calling 3 CD GCon to form the spherical joint (SJ)
gcon_inp.p1 = p1;
gcon_inp.p1dot = p1dot;
gcon_inp.p2 = p2;
gcon_inp.p2dot = p2dot;
gcon_inp.ft = 0;
gcon_inp.dft = 0;
gcon_inp.ddft = 0;
gcon_inp.r1 = r1;
gcon_inp.r2 = r2;
gcon_inp.s1_p = s1_p;
gcon_inp.s2_q = s2_q;
gcon_inp.i = i;
gcon_inp.j = j;

for w = 1:3
    c = zeros(3,1);
    c(w) = 1;
    gcon_inp.c = c;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getcd(gcon_inp);
end

% calling 2 DP1 GCon to add to SJ to form Revolute joint
% defining such that x'-y' plane is perpendicular to X-axis 
a_2 = [1,0,0].';                                % x-axis in G-RF
a_1 = [1,0,0].';                                % x-axis in L-RF
b_1 = [0,1,0].';                                % y-axis in L-RF

gcon_inp.a_2 = a_2;
gcon_inp.a_1 = a_1; 
w = w+1;
[phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getdp1(gcon_inp);
gcon_inp.a_1 = b_1; 
w = w+1;
[phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getdp1(gcon_inp);

% calling another DP1 for the driving constraint
a_1 = [0,1,0].';
a_2 = [0,0,-1].';
gcon_inp.a_2 = a_2;
gcon_inp.a_1 = a_1; 
gcon_inp.ft = ft;
gcon_inp.dft = dft;
gcon_inp.ddft = ddft;
w = w+1;
[phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getdp1(gcon_inp);

% euler parameter normalization constraint
w = w+1;
[phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getpnorm(p1);

% assembling the Jacobian matrix by combining dphi_dr and dphi_dp
Phi_q = [dphi_dr,dphi_dp];                        
q(:,n) = [r1;p1];
qbar(:,n) = A1(theta0)\q(1:3,n);      % in L-RF
qdot(:,n) = Phi_q\nu;
q_dot(:,n) = A1(theta0)\qdot(1:3,n); % in L-RF
qddot(:,n) = Phi_q\gamma;
q_ddot(:,n) = A1(theta0)\qddot(1:3,n);% in L-RF
end
%%
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
open(myVideo)
set(gcf,'WindowState','maximized')
O = [0 0];
y_pos = q(2,:);
z_pos = q(3,:);

subplot(3,3,4);
axis([min(ts) max(ts) min(q(2,:)) max(q(2,:))]);
xlabel('time(s)');ylabel('Displacement(m)')
title('Displacement along Y')
box on
grid on
ani4=animatedline('Color','r');
subplot(3,3,5)
axis([min(ts) max(ts) min(qdot(2,:)) max(qdot(2,:))])
xlabel('time(s)');ylabel('velocity(m/s)')
title('Velocity along Y')
grid on
box on
ani5=animatedline('Color','k');

subplot(3,3,6)
axis([min(ts) max(ts) min(qddot(2,:)) max(qddot(2,:))])
xlabel('time(s)');ylabel('acceleration(m^2/s)')
title('Acceleration along Y')
grid on
box on
ani6=animatedline('Color','k');
subplot(3,3,7)
axis([min(ts) max(ts) min(q(3,:)) max(q(3,:))]);
xlabel('time(s)');ylabel('Displacement(m)')
title('Displacement along Z')
grid on
box on
ani7=animatedline('Color','k');
subplot(3,3,8)
axis([min(ts) max(ts) min(qdot(3,:)) max(qdot(3,:))])
xlabel('time(s)');ylabel('velocity(m/s)')
title('Velocity along Z')
grid on
box on
ani8=animatedline('Color','k');
subplot(3,3,9)
axis([min(ts) max(ts) min(qddot(3,:)) max(qddot(3,:))])
xlabel('time(s)');ylabel('accleration(m^2/s)')
title('Acceleration along Z')
grid on
box on
ani9=animatedline('Color','k');

for k=1:length(ts)
    time = ts(k);
    subplot(3,3,2);
    title(['Pendulum with Revolute Joint, Time: ',num2str(time),' sec'])
    axis([min(y_pos)-0.5 max(y_pos)+0.5 min(z_pos)-0.5 1])
    P = [y_pos(k) z_pos(k)];
    O_circ = viscircles(O,0.02);
    pendulum = line([O(1) y_pos(k)], [O(2) z_pos(k)],'lineWidth',5);
    xlabel('Position along Y(m)');ylabel('Position along Z (m)')
    grid on
    box on
    hold on
    
    addpoints(ani4,ts(k),q(2,k))
    addpoints(ani5,ts(k),qdot(2,k)) ;
    addpoints(ani6,ts(k),qddot(2,k))
    addpoints(ani7,ts(k),q(3,k)) ;
    addpoints(ani8,ts(k),qdot(3,k))
    addpoints(ani9,ts(k),qddot(3,k))
    
    drawnow
    grid on
    box on
    pause(0.01);
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
    if k<length(ts)
        delete(pendulum)
    end
    
end
close(myVideo)
print('-dpng','-painters','PendulumSimulation.png');
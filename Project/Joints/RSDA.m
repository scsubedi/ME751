% Takes input for the Rotational Spring Damper Actuator
% Outputs the torque on the body on a plot as a function of time
% One of the body is grounded

A1 = @(th)[0,  0,      1;                        % Rotation matrix as a function of angle
    sin(th),  cos(th), 0;
    -cos(th), sin(th), 0];

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);

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
    torque = k*(theta0-theta(1)) + c0* thetaDot + (theta0 + thetaDot*sin(2*pi*t));
    % Torque acting on body 1 is given by
    torque1 = torque.*r1;
    torque2 = torque.*r2;
    nRSDA_1(n,1) =norm(torque1);
end

figure(1)
plot(ts,nRSDA_1)
grid on
title('Concentrated Torque due to Rotational Spring-Damper-Actuator(RSDA')
xlabel('Time(s)')
ylabel('Torque')
print('-dpng','-painters','RSDA.png');
disp('Plot for Torque generated')

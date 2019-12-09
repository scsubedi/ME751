%%%%%%%%%%%%%
%Assignment 6, Problem 3

clear all
close all
clc

points1 = [2 0 0;4 0 0];

q =  [0   1.4142   -1.4142    0.6533    0.2706    0.6533    0.2706]';
qdot = q;
qdotdot = q;
time = 1:10e-6:10;
time = linspace(1,10,3);
%generate Plots
for i = 1:size(points1,1)
    pos = zeros(3,size(q,2));
    mu = zeros(3,size(qdot,2));
    gamma = zeros(3,size(qdotdot,2));
    for t = 1:size(time,1)
        pos(:,t) = q(1:3,t)+Amatrix(q(4),q(5:7))*points1(i,:)';
        mu(:,t) = qdot(1:3,t)+Bmatrix(q(4:7),points1(i,:)')*qdot(4:7,t);
        gamma(:,t) = qdotdot(1:3,t)+Bmatrix(qdot(4:7),points1(i,:)')*qdot(4:7,t)+Bmatrix(q(4:7),points1(i,:)')*qdotdot(4:7,t);
    end
    
    if(i==1)
        titletext =  'Kinematic Analysis at the center of Pendulum, OPrime';
    else
        titletext = 'Kinematic Analysis at the center of rotation of Pendulum, Point Q';
    end
    
    figure()
    subplot(1,3,1);
    plot(time,pos,'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Position');
   
   sgtitle(titletext);
    subplot(1,3,2);
    plot(time,mu,'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Velocity');
  
%     title(titletext);
    subplot(1,3,3);
    plot(time,gamma,'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Acceleration');
%     title(titletext);
    
end

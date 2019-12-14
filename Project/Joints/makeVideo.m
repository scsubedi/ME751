%% Creates a video of the position, velocity, acceleration of th pendulum
function A = makeVideo(output)

Phi_q = output{1};
q = output{2};
qbar = output{3};
qdot = output{4};
q_dot = output{5};
qddot = output{6};
q_ddot = output{7};
ts = output{8};


myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 15; 
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
    P_circ = viscircles([y_pos(k),z_pos(k)],0.1);
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
        delete(P_circ)
    end
    
end
close(myVideo)
% print('-dpng','-painters','PendulumSimulation.png');
end
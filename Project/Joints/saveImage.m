function A = saveImage(L, output)
figure(2)
Phi_q = output{1};
q = output{2};
qbar = output{3};
qdot = output{4};
q_dot = output{5};
qddot = output{6};
q_ddot = output{7};
ts = output{8};


set(gcf,'WindowState','maximized')
y_pos = q(2,:);
z_pos = q(3,:);

subplot(3,3,2);
title('Pendulum with Revolute Joint')
axis([min(y_pos)-0.5 max(y_pos)+0.5 min(z_pos)-0.5 1])
O_circ = viscircles([0 0],0.05);
pendulum = line([0 0], [0 -L],'lineWidth',5);
Big_circ = viscircles([0 -L],0.1);
xlabel('Position along Y(m)');ylabel('Position along Z (m)')
grid on
box on


y_pos = q(2,:);
z_pos = q(3,:);

subplot(3,3,4);
axis([min(ts) max(ts) min(q(2,:)) max(q(2,:))]);
plot(ts,q(2,:))
xlabel('time(s)');ylabel('Displacement(m)')
title('Displacement along Y')
box on
grid on

subplot(3,3,5)
plot(ts,qdot(2,:))
axis([min(ts) max(ts) min(qdot(2,:)) max(qdot(2,:))])
xlabel('time(s)');ylabel('velocity(m/s)')
title('Velocity along Y')
grid on
box on

subplot(3,3,6)
plot(ts,qddot(2,:))
axis([min(ts) max(ts) min(qddot(2,:)) max(qddot(2,:))])
xlabel('time(s)');ylabel('acceleration(m^2/s)')
title('Acceleration along Y')
grid on
box on

subplot(3,3,7)
plot(ts,q(3,:))
axis([min(ts) max(ts) min(q(3,:)) max(q(3,:))]);
xlabel('time(s)');ylabel('Displacement(m)')
title('Displacement along Z')
grid on
box on

subplot(3,3,8)
plot(ts,qdot(3,:))
axis([min(ts) max(ts) min(qdot(3,:)) max(qdot(3,:))])
xlabel('time(s)');ylabel('velocity(m/s)')
title('Velocity along Z')
grid on
box on

subplot(3,3,9)
plot(ts,qddot(3,:))
axis([min(ts) max(ts) min(qddot(3,:)) max(qddot(3,:))])
xlabel('time(s)');ylabel('accleration(m^2/s)')
title('Acceleration along Z')
grid on
box on
print('-dpng','-painters','PendulumSimulation.png');
end
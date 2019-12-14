% TSDA
% Takes input from the 'input.txt' file
% Setting up translational spring-damper-actuator acting between point Pi
% on body i,and Pj on body j
% Outputs a plot on the variation of force acting on the ungrounded body as
% a function of time.


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

for i=1:length(ts)
    h = lij + lijDot + (sin(ts(i)));
    f = real(vpa(k*(lij - l0) + c0*lijDot' - h'));
    Ftsda = f.*eij;
    nTSDA_1(i,1) = norm(Ftsda);
end

figure(1)
plot(ts,nTSDA_1)
grid on
title('Concentrated Force due to Translational Spring-Damper-Actuator(TSDA)')
xlabel('Time(s)')
ylabel('Torque')
print('-dpng','-painters','TSDA.png');
disp('Plot for Forces generated')
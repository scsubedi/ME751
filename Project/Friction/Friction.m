clear all
close all
clc

% Friction by DEM-P approach

% Body Attributes
r_i = [1 3 0]';
r_j = [0 4 0]';
e_i = [1 0 0 ]';
e_iDot = [0 0 0]';
e_j = [0 1 0]';
e_jDot = [0 0 0]';
e_i0 = sqrt(e_i(1)^2 + e_i(2)^2 + e_i(3)^2);
e_i0Dot=-e_i'*e_iDot/e_i0;
e_j0 = sqrt(e_j(1)^2 + e_j(2)^2 + e_j(3)^2);
e_j0Dot=-e_j'*e_jDot/e_j0;
p_i = [e_i0;e_i];
p_iDot = [e_i0Dot;e_iDot];
p_j = [e_j0;e_j];
p_jDot = [e_j0Dot;e_jDot];
A_i = Amatrix(p_i);
A_j = Amatrix(p_j);
r_iDot = [0.5 0.5 0.5]';
r_jDot = [0.5 0.5 0.5]';

s_iBar = [1,2,0]';
s_jBar = [-1,3,4]';

mi = 1/1000; % 1gm
Radius_i = 1/1000; %1mm
I = 2/3*mi*Radius_i^2;

n = [1 3 0];
n_i = n./sqrt(n(1)^2 + n(2)^2 + n(3)^2);

mj = 1/1000;
Radius_j =  1/1000;

R_ij = Radius_j*Radius_i/(Radius_j+Radius_i);

m_ij = mi*mj/(mi + mj);

vi = r_iDot;
vj = r_jDot;


theta = 0:2*pi;
w_i = vi\r_i;
w_j = vj\r_j;

vij = vj + w_j*r_j - vj - w_i*r_i;

v_ij_n = (dot(vij,n_i)*n_i)';
v_ij_t = vij - v_ij_n;


q = [r_i;p_i;r_j;p_j];
qDot = [r_iDot;p_iDot;r_jDot;p_jDot];


v = [r_iDot;w_i;r_jDot;w_j];

G_i = getGfromP(p_i);
Gp_i = 1/2*G_i;
G_j = getGfromP(p_j);
Gp_j = 1/2*G_j;

Lq = diag([eye(3) Gp_i eye(3) Gp_j])
qdot = Lq.*v'




%% Friction and Contact

q = [r_i]






%%
w_2 = cos(theta)

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);

% variables derived from above input
% for n = 1:length(ts)
%     t = ts(n);
%     omega = [dtheta(t),0,0].';

%%
% Normal Friction Force
Fij_n = sqrt(del_ij/(2*R_ij))*(k_n*del_ij*n_ij -gamma_n*m_ij*v_ij) + Fc

% Tangential Friction Force
Fij_t = sqrt(del_ij/(2*R_ij))*(-k_t*del_ij_t - gamma_t*m_ij*v_ij_t)




function G = getGfromP(p)

G = [-p(2) p(1) p(4) -p(3)
    -p(3) -p(4) -p(1) p(2)
    -p(4) p(3) -p(2) p(1)];
end


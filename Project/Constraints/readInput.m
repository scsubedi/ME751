%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace. This file is used to display results for different constaints

clear;clc;
close all;


% Reading the input file by avoiding white spaces and delimiters
syms t
filename  ='input_Constraints';
str = fileread([filename '.txt']);
val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
val = strtrim(vertcat(val{:}));
num = cellfun(@(s)sscanf(s,'%f:'),val(:,2),'UniformOutput',false);
len = cellfun('length',num);
vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);


%%
c = str2num(val{3,2})';
r_i = str2num(val{4,2})';
r_iDot = str2num(val{5,2})';
ei = str2num(val{6,2})';
eiDot = str2num(val{7,2})';
a_iBar = str2num(val{8,2})';

f = str2func(num2str(val{9,2}));

time =str2num(val{10,2});
r_j = str2num(val{11,2})';
r_jDot = str2num(val{12,2})';
ej = str2num(val{13,2})';
ejDot = str2num(val{14,2})';
a_jBar = str2num(val{15,2})';
flag = str2num(val{16,2})';
groundStatus = str2num(val{17,2})';
Constraint = str2num(val{18,2}).';

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
ft=f(time);
fDott=vpa(subs(df,t,time));
fDotDott = vpa(subs(dff,t,time));
dij = r_j + A_j*s_jBar - r_i - A_i*s_iBar;
dijDot =  r_jDot - r_iDot + Bmatrix(p_j,s_jBar)*p_j - Bmatrix(p_i,s_iBar)*p_i;

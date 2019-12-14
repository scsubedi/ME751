clear all
close all
clc
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
syms t
vargin.pmin = val{1,2};
% vargin.pmax = val{}
vargin.theta =str2sym( val{7,2})
vargin.df = diff(vargin.theta,t)
vargin.dff = diff(vargin.df,t)

theta = @(t)vargin.theta
dtheta = @(t)vargin.df
ddtheta = @(t)vargin.dff
t = 5;
R =eval(ddtheta(t))
% vargin.A = val{8,2}
% t = 0:1:90;
% t = 90;
% theta = @t vargin.theta

%%

for i = 1:length(t)-1
    t = t(i+1)
theta(t)= vargin.theta
R(i,:) = eval(theta(90))
end
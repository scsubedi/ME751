%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace

clear;clc;
close all;


% Reading the input file by avoiding white spaces and delimiters
syms t
filename  ='input';
str = fileread([filename '.txt']);
val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
val = strtrim(vertcat(val{:}));
num = cellfun(@(s)sscanf(s,'%f:'),val(:,2),'UniformOutput',false);
len = cellfun('length',num);
vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);

i = str2num(val{3,2})';
j = str2num(val{4,2})';
L = str2num(val{5,2})';
s1_p = L*str2num(val{6,2})';
s2_q = L*str2num(val{7,2})';

ts = str2num(val{8,2});

type = val{12,2};

joint = val{10,2};
% saveResults = val{11,2};






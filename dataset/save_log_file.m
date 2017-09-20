clear all;close all;clc;
% Format for "filename.txt"
% | test number | sample number | methods: IQP = 0, ILP = 1 | solver: CPLEX=1, GUROBI=2,
% MOSEK=3, BARON=4,XPRESS=5 |
% e.g., 0800111.txt
%   | 08 | 001 | 0 | 1 | 
%   test number is 8, smaple number is 1, method is IQP, solver is CPLEX
test_num = 60;
sample_num = 100;


method = 0; % IQP = 0, ILP = 1

solver = 4; % solver: 1:CPLEX, 2:GUROBI, 3:MOSEK, 4:BARON, 5:EXPRESS

str1 = sprintf('%d',test_num);
str2 = sprintf('%d',sample_num);
if numel(str1) == 1
    str1 = sprintf('0%d',test_num);
end
switch(numel(str2))
    case 1
        str2 = sprintf('00%d',sample_num);
    case 2
        str2 = sprintf('0%d',sample_num);
    case 3
end
str = strcat(str1,str2,sprintf('%d',method),sprintf('%d',solver),'.txt');


diary(str)
diary on;
diary off;






clear all;close all;clc;

numData = 100;
for i = 5:60
    v = 2*i+1;
    for j = 1:numData
        pos{i}{j}= rand(v,2);
    end
end
clear all,close all,clc;
n = 5 ;             % number of custumers(n)
k = 2;              % capacity   
v = 2*n+1;

rng('shuffle');     % random seed: shuffle
vert = rand(2*n+1,2);
   
% if ~exist('map.mat'),
%    rng('shuffle');     % random seed: shuffle
%    vert = rand(2*n+1,2);
% else
%     load('map.mat');
% end
[ para,constr ] = pre_process( n,k,vert,[],1);
fun = @(x)quadcost(x,para.Q);
constraint = @(x)tsp_constr(x,constr);

% TEST_GRAD_MAIN( fun,constr,v );

% X_real = [     0     0     0     0     0     1     0     0     0     0     0;
%                0     0     0     0     0     0     0     0     0     0     1;
%                0     0     1     0     0     0     0     0     0     0     0;
%                0     0     0     0     0     0     0     1     0     0     0;
%                0     0     0     1     0     0     0     0     0     0     0;
%                0     1     0     0     0     0     0     0     0     0     0;
%                0     0     0     0     0     0     1     0     0     0     0;
%                0     0     0     0     0     0     0     0     1     0     0;
%                0     0     0     0     1     0     0     0     0     0     0;
%                0     0     0     0     0     0     0     0     0     1     0;
%                1     0     0     0     0     0     0     0     0     0     0];
% x0 = reshape(X_real',v^2,1);
x0 = reshape([zeros(v-1,1),eye(v-1);1,zeros(1,v-1)]',v^2,1);

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
options.MaxIterations = 600;
[x,fval,eflag,output,lambda,grad] = fmincon(fun,x0,...
    [],[],[],[],zeros(v^2,1),ones(v^2,1),constraint,options);
clear all,close all,clc;
n = 5 ;             % number of custumers(n)
k = 2;              % capacity   
v = 2*n+1;

if ~exist('map.mat'),
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load('map.mat');
end
[ para,constr ] = pre_process( n,k,vert,[],1);
fun = @(x)quadcost(x,para.Q);
constraint = @(x)tsp_constr(x,constr);

% TEST_GRAD_MAIN( fun,constr,v );

X_real = [     0     0     0     0     0     1     0     0     0     0     0;
               0     0     0     0     0     0     0     0     0     0     1;
               0     0     1     0     0     0     0     0     0     0     0;
               0     0     0     0     0     0     0     1     0     0     0;
               0     0     0     1     0     0     0     0     0     0     0;
               0     1     0     0     0     0     0     0     0     0     0;
               0     0     0     0     0     0     1     0     0     0     0;
               0     0     0     0     0     0     0     0     1     0     0;
               0     0     0     0     1     0     0     0     0     0     0;
               0     0     0     0     0     0     0     0     0     1     0;
               1     0     0     0     0     0     0     0     0     0     0];
x0 = reshape(X_real',v^2,1);
options = optimoptions(@fmincon,'Algorithm','active-set',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
[x,fval,eflag,output,lambda,grad] = fmincon(fun,x0,...
    [],[],[],[],[],[],constraint,options);
% solving tsppd using manifold  (This idea is doomed, since trust region only finds the local minimum base on the cost function with regularizer. For larger )
clear all,close all,clc;
n = 1 ;             % number of custumers(n)
k = 1;              % capacity   
v = 2*n+1;

file = strcat('map',num2str(n),'.mat');   
if ~exist(file),
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load(file);
end
[ para,A,b,Base1,Base2,totaldim,P,e1 ] = pre_process_new( n,k,vert);
manifold = symfixedrankYYfactory(totaldim,1);
problem.M = manifold;
N = 1000000;
problem.cost = @(y)cost( y,para.Q,A,b,N,para.aux );
problem.egrad = @(y)egrad(y,para.Q,A,Base1,Base2,P,N,para.aux);
problem.ehess = @(y,u)ehess(y,u,para.Q,A,Base1,Base2,P,N,para.aux);
checkgradient(problem); pause;
checkhessian(problem); pause;
y = trustregions(problem)
y_real = y*y(1);
solution = reshape(y_real(2:1+v^2),v,v)'





% truth = [1 0 1 0 0 0 1 1 0 0 1 1]';
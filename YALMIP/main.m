clear all,close all,clc;
n = 5 ;             % number of custumers(n)
k = 2;              % capacity   
v = 2*n+1;

file = strcat('map',num2str(n),'.mat');   
if ~exist(file),
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load(file);
end
[ para,constr ] = pre_process( n,k,vert);
x = binvar(v^2,1);
constraint = [];
obj = x'*para.Q*x;
% obj = trace(para.Q*x*x');

for i = 1:length(constr.eq.A),
%    constraint = [constraint;constr.eq.A{i}*(x.^2)==constr.eq.B{i}]; 
   constraint = [constraint;constr.eq.A{i}*x==constr.eq.B{i}]; 
end
for i = 1:length(constr.ineq.A),
%    constraint = [constraint;constr.ineq.A{i}*(x.^2)<=constr.ineq.B{i}]; 
   constraint = [constraint;constr.ineq.A{i}*x<=constr.ineq.B{i}]; 
end
M = [1 x'.^2;x.^2,x*x'];
% constraint = [constraint;M>=0];
% ops = sdpsettings('verbose',1,'solver','mosek');
% ops = sdpsettings('verbose',1,'solver','fminsdp');
ops = sdpsettings('verbose',1);
optimize(constraint,obj,ops)
solution = value(x);
solution_as_perm = reshape(solution,v,v)'

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

clear all,close all,clc;
n = 5 ;             % number of custumers(n)
k = 1;              % capacity   
v = 2*n+1;

file = strcat('map',num2str(n),'.mat');   
if ~exist(file),
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load(file);
end
[ para,constr ] = pre_process( n,k,vert);
fun = @(x)object(x,v,para.Q);
A = [];
b = [];
Aeq = [];
beq = [];
for i = 1:length(constr.eq.A),
    Aeq = [Aeq;constr.eq.A{i}];
    beq = [beq;constr.eq.B{i}];
end
for i = 1:length(constr.ineq.A),
    A = [A;constr.ineq.A{i}];
    b = [b;constr.ineq.B{i}];
end


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

%  x0_real = reshape(X_real',v^2,1);          
%  x = fmincon(fun,x0_real,A,b,Aeq,beq)

%% generate x0
x = sdpvar(v^2,1);
constraint = [x>=0;x<=1];
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
ops = sdpsettings('verbose',1);
optimize(constraint,obj,ops)
x0 = value(x);

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'CheckGradients',true,'SpecifyObjectiveGradient',true);

% x0 = [0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0]';
%  x = fmincon(fun,x0,A,b,Aeq,beq);
 x = fmincon(fun,x0+0.1*rand(size(x0)),A,b,Aeq,beq);
 x_fmincon = reshape(x,v,v)'
  
 
 
 
%% NLP with yalmip
 x = sdpvar(v^2,1);
 obj = x'*para.Q*x+10000000*ones(1,v^2)*((x.*(x-1)).^2);
 constraint = [];
 for i = 1:length(constr.eq.A),
%    constraint = [constraint;constr.eq.A{i}*(x.^2)==constr.eq.B{i}]; 
   constraint = [constraint;constr.eq.A{i}*x==constr.eq.B{i}]; 
end
for i = 1:length(constr.ineq.A),
%    constraint = [constraint;constr.ineq.A{i}*(x.^2)<=constr.ineq.B{i}]; 
   constraint = [constraint;constr.ineq.A{i}*x<=constr.ineq.B{i}]; 
end
ops = sdpsettings('verbose',1,'usex0',x0);
 optimize(constraint,obj,ops)
 x_yalmip = reshape(value(x),v,v)'
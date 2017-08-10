%% lassere
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
[ para,constr ] = pre_process( n,k,vert);

% create sosprog
prog = spotsosprog;
x = msspoly('x',v^2);
prog = prog.withIndeterminate(x);
f = x'*para.Q*x; % real cost funnction
degree = 2; % degre used in sosOnK
[prog,lambda]=prog.newFree(1);

% construct semi-algebraic set K containing x(x-1)=0 and Ax-b>=0
K = [x.*(x-1);x.*(1-x)];
for i=1:length(constr.eq.A),
   K = [K; constr.eq.A{i}*x-constr.eq.B{i};-constr.eq.A{i}*x+constr.eq.B{i}]; 
end
for i=1:length(constr.ineq.A),
   K = [K; -constr.ineq.A{i}*x+constr.ineq.B{i}]; 
end
[prog,tk,~] = sosOnK(prog,f+lambda,x,K,degree);
options = spot_sdp_default_options();
options.verbose = 2;
sol = prog.minimize(lambda, @spot_mosek_sos, options);
obj_sos = sol.eval(lambda)




%% YALMIP 
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
optimize(constraint,obj,ops);

%% results:
obj_YALMIP = value(obj)
obj_sos = -double(obj_sos)
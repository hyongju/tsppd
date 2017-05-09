function [ y,yeq,grady,gradyeq ] = tsp_constr_new( x,constr,real_length )
%This function should be used as a function handle for all the constraints.
%Typically for our problem, the first constraints are equality constraints,
%and the rest two are inequality constraints. 
%   input:  x      --  variable of the function handle.
%           constr --  structure that construced by 'pre_process'.
%   output: y      --  inequality constraints.
%           yeq    --  equality constraints.
%           grady  --  gradient of y.
%           gradyeq--  gradient of yeq.

%%%%%%%%%%%%%%%%%%%%%%% pre-process constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
% equality constraints
if ~exist('real_length', 'var') || isempty(real_length)
     real_length=length(x);
end
v = sqrt(real_length);% v has to be an integer
n = length(constr.eq.A);
yeq = [];
gradyeq = [];
for i = 1:n,
    for j = 1:size(constr.eq.A{i},1),
        yeq(end+1) = constr.eq.A{i}(j,:)*(x(1:real_length).^2)- constr.eq.B{i}(j,:);
        gradyeq(:,end+1) = [2*diag(x(1:real_length))*constr.eq.A{i}(j,:)';zeros(length(x)-real_length,1)];
    end
end

% inequality constraints
n = length(constr.ineq.A);
y = [];
grady = [];
for i = 1:n,
    for j = 1:size(constr.ineq.A{i},1),
        y(end+1) = constr.ineq.A{i}(j,:)*(x(1:real_length).^2)- constr.ineq.B{i}(j,:);
        grady(:,end+1) = [2*diag(x(1:real_length))*constr.ineq.A{i}(j,:)';zeros(length(x)-real_length,1)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% constraint: X(v,1)==1 %%%%%%%%%%%%%%%%%%%%%%%%%%
yeq(end+1) = x(v*(v-1)+1)-1;
gradyeq(:,end+1) = [zeros(v*(v-1),1);1;zeros(v-1,1);zeros(length(x)-real_length,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%% constraint: lambda==v+1 %%%%%%%%%%%%%%%%%%%%%%%%%
lambda = v+1;
eig_vec = x(real_length+1:end);
x_real = x(1:real_length);
M = [1,x_real'.^2;x_real.^2,x_real*x_real'];
yeq(end+1:end+v^2+1) = M*eig_vec - lambda*eig_vec;
gradyeq(:,end+1) = [2*x_real.*eig_vec(2:end);1-lambda;x_real.^2];
for i = 1:real_length,
   gradyeq(:,end+1) = [x_real(i)*eig_vec(2:end);x_real(i)^2;x_real(i)*x_real];
   gradyeq(i,end) = gradyeq(i,end)+sum(x_real.*eig_vec(2:end))+2*eig_vec(1)*x_real(i);
   gradyeq(real_length+1+i,end) = gradyeq(real_length+1+i,end)-lambda;
end

end


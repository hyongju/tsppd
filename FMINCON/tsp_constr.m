function [ y,yeq,grady,gradyeq ] = tsp_constr( x,constr )
%This function should be used as a function handle for all the constraints.
%Typically for our problem, the first constraints are equality constraints,
%and the rest two are inequality constraints. 
%   input:  x      --  variable of the function handle.
%           constr --  structure that construced by 'pre_process'.
%   output: y      --  inequality constraints.
%           yeq    --  equality constraints.
%           grady  --  gradient of y.
%           gradyeq--  gradient of yeq.

% equality constraints
n = length(constr.eq.A);
yeq = [];
gradyeq = [];
for i = 1:n,
    for j = 1:size(constr.eq.A{i},1),
        yeq(end+1) = constr.eq.A{i}(j,:)*(x.^2)- constr.eq.B{i}(j,:);
        gradyeq(:,end+1) = 2*diag(x)*constr.eq.A{i}(j,:)';
    end
end

% inequality constraints
n = length(constr.ineq.A);
y = [];
grady = [];
for i = 1:n,
    for j = 1:size(constr.ineq.A{i},1),
        y(end+1) = constr.ineq.A{i}(j,:)*(x.^2)- constr.ineq.B{i}(j,:);
        grady(:,end+1) = 2*diag(x)*constr.ineq.A{i}(j,:)';
    end
end

end


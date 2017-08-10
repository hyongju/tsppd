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


v = sqrt(length(x));
yeq = [];
gradyeq = [];
y = [];
grady = [];
% equality constraints
n = length(constr.eq.A);
for i = 1:n,
    for j = 1:size(constr.eq.A{i},1),
        yeq(end+1) = constr.eq.A{i}(j,:)*(x.^2)- constr.eq.B{i}(j,:);
        if nargout>2,
        gradyeq(:,end+1) = 2*diag(x)*constr.eq.A{i}(j,:)';
        end
    end
end

% inequality constraints
n = length(constr.ineq.A);
for i = 1:n,
    for j = 1:size(constr.ineq.A{i},1),
        y(end+1) = constr.ineq.A{i}(j,:)*(x.^2)- constr.ineq.B{i}(j,:);
        if nargout>2,
        grady(:,end+1) = 2*diag(x)*constr.ineq.A{i}(j,:)';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% constraint: X(v,1)==1 %%%%%%%%%%%%%%%%%%%%%%%%%%
yeq(end+1) = x(v*(v-1)+1)-1;
if nargout>2,
    gradyeq(:,end+1) = [zeros(v*(v-1),1);1;zeros(v-1,1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%% constraint: lambda==v+1 %%%%%%%%%%%%%%%%%%%%%%%%%
lambda = v+1;
eig_vec = [1;x]/sqrt(lambda);
M = [1,x'.^2;x.^2,x*x'];
yeq(end+1:end+v^2+1) = M*eig_vec - lambda*eig_vec;
if nargout>2,
    gradyeq(:,end+1) = 3*x.^2/sqrt(lambda);
    temp = (sum(x.^2)*ones(v^2,1)+2*x-lambda)/sqrt(lambda);
    temp = 2*x*x'/sqrt(lambda)+diag(temp);
    gradyeq = [gradyeq,temp'];
end



scale = 1;
y = scale*y;yeq = scale*yeq;grady = scale*grady;gradyeq = scale*gradyeq;

end


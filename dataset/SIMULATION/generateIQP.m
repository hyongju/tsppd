function [ ] = generateIQP( c,v,k,ops,n,sample_num )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

constr = [];

x = binvar(v^2,1);
for i = 1:v,
    x(v*(i-1)+1) = 0;
end
x(end-v+1:end) = [1;zeros(v-1,1)];
X = reshape(x,v,v)';

% construct cost function
Q = zeros(v^2,v^2);
for i = 1:v
    Q((i-1) * v + 1:(i-1) * v + v,(i-1) * v + 1:(i-1) * v + v) =c;
end
for i = 1:v-1
    Q(i*v+1:i*v+v,(i-1)*v+1:(i-1)*v+v) = c/2;
    Q((i-1)*v+1:(i-1)*v+v,i*v+1:i*v+v) = c/2;
end

Q((v-1)*v + 1:(v-1)*v + v,1:v) = c/2;
Q(1:v,(v-1)*v + 1:(v-1)*v + v) = c/2;
obj = x'*Q*x;
 
% d(i) = +1 (pickup), -1 (deliver)
d = [0;ones(n,1);-ones(n,1)];


% constraint: permutation matrix
for i = 1:v-1,
   constr = [constr; sum(X(i,:))==1];
   if i ~= v-1,
       constr = [constr; sum(X(:,i+1))==1];
   end
end

% constraint: capacity
L = tril(ones(v));
L(1:k,:) = [];
L(end-k+1:end,:) = [];
constr = [constr;L*X*d<=k];

% constraint: precedence
N = 1:v;
for i = 2:1+n,
   ei = zeros(v,1);
   eni = zeros(v,1);
   ei(i) = 1;
   eni(i+n) = 1;
   constr = [constr;N*X*(ei-eni)<=0];
end

% constraint: goes back to node 1 in the end
constr = [constr; X(v,1)==1];


% run IQP with 5 solvers
if n <= 50,           
    ops.solver = 'cplex';
    [txtfile,matfile] = generateFilename(n,sample_num,0,1);
    diary(txtfile);diary on
    IQP = optimize(constr,obj,ops);
    diary off
    save(matfile,'IQP','n','sample_num','k');
end
if n <= 25,
    ops.solver = 'gurobi';
    [txtfile,matfile] = generateFilename(n,sample_num,0,2);
    diary(txtfile);diary on
    IQP = optimize(constr,obj,ops);
    diary off
    save(matfile,'IQP','n','sample_num','k');
end
if n <= 25,
    ops.solver = 'mosek';
    [txtfile,matfile] = generateFilename(n,sample_num,0,3);
    diary(txtfile);diary on
    IQP = optimize(constr,obj,ops);
    diary off
    save(matfile,'IQP','n','sample_num','k');
end
if n <= 12,
    ops.solver = 'baron';
    [txtfile,matfile] = generateFilename(n,sample_num,0,4);
    diary(txtfile);diary on
    IQP = optimize(constr,obj,ops);
    diary off
    save(matfile,'IQP','n','sample_num','k');
end
if n <= 30,
    ops.solver = 'xpress';
    [txtfile,matfile] = generateFilename(n,sample_num,0,5);
    diary(txtfile);diary on
    IQP = optimize(constr,obj,ops);
    diary off
    save(matfile,'IQP','n','sample_num','k');
end   

end


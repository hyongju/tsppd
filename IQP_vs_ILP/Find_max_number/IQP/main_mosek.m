% iqp vs ilp
clear all,close all,clc;

n = 20;             % number of custumers(n)
k = 2;              % capacity   
v = 2*n+1;

% construct map and cost map
file = strcat('map',num2str(n),'.mat');   
if ~exist(file),
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load(file);
end

% cost matrix
c = zeros(v,v);
for i = 1:size(c,1)
    for j = 1:size(c,2)
        c(i,j) = norm(vert(i,:)-vert(j,:));
    end
end
% c = c + eye(v); % add identity to prevent zeros from appearing in the diagonal

% make c PSD...

for i = 1:size(c,1)
    sum_c(i) = 0;
    for j = 1:size(c,2)
        sum_c(i) = sum_c(i) + abs(c(i,j));
    end
    c(i,i) = c(i,i) + sum_c(i);
    cM(i) = c(i,i);
end
for i = 1:size(c,1)
    c(i,i) = max(cM)*0.5;
end
all(eig(c)>=0)  % check if c is PSD

%% IQP
constr = [];

x = binvar(v^2,1);
for i = 1:v,
    x(v*(i-1)+1) = 0;
end
x(end-v+1:end) = [1;zeros(v-1,1)];
X = reshape(x,v,v)';
x

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

% ops = sdpsettings('verbose',1,'solver','mosek');
% ops = sdpsettings('verbose',0,'solver','cplex');

ops = sdpsettings('solver','mosek','verbose',2,'debug',1,'savesolveroutput',1);
relative_gap = 0.035;
maxtime = 7200;

% setting cplex
ops.cplex.mip.tolerances.mipgap=relative_gap;
ops.cplex.MaxTime = maxtime;

% setting gurobi
ops.gurobi.MIPGap = relative_gap;
ops.gurobi.TimeLimit = maxtime;

% setting mosek
ops.mosek.MSK_DPAR_MIO_MAX_TIME = maxtime;
ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP = relative_gap;



% setting baron
ops.baron.cplexlibname = 'cplex1271.dll';  % probably we should not use cplex here
ops.baron.epsr = relative_gap;
ops.baron.maxtime = maxtime;

% setting xpress
ops.xpress.MIPRELSTOP = relative_gap;
ops.xpress.MAXTIME = maxtime;





% ops = sdpsettings('solver','baron','verbose',2);
% ops.baron.cplexlibname = 'cplex1271.dll';  % probably we should not use cplex here
% ops.baron.epsr = 0.035;

% ops = sdpsettings('solver','xpress','verbose',2);
% ops.xpress.MIPRELSTOP = 0.1;


i = 5;
j = 1;
% diary(sprintf('%d%d.txt',i,j));
diary('mosek_IQP.txt');
diary on
IQP = optimize(constr,obj,ops)
diary off
obj_IQP = value(obj) - c(1,1)*v
solution_IQP = value(X);
tour = round(solution_IQP*[1:v]');
tour = [1;tour];






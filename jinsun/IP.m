clear all;close all;clc

n = 10 ;             % number of custumers(n)
k = 2;              % capacity    
rng('shuffle');     % random seed: shuffle
alph = 1;
% generate random vehicle, customer pickup and delivery locations from [0,1]x[0,1]
% vehicle node: 1
% pick-up nodes: 2,...,n+1
% delivery nodes: n+2,...,2n+1
vert = rand(2*n+1,2);
load('map.mat'); % for case n=10
v = 2*n + 1;        % number of vertices |V| = 2n + 1

% cost matrix (c)
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
    c(i,i) = max(cM);
end
all(eig(c)>=0)  % check if c is PSD
% adjacency matrix for a default tour: 1->2->3->..->2n->2n+1
A0 = eye(2*n+1);
A0 = [A0(2:2*n+1,:);A0(v,:)] + alph*eye(v); % added identity to prevent zeros from appearing in the diagonal

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

%% integer programming
constr = [];
x = binvar(v^2,1);
X = [];
for i = 1:v,
    X = [X,x(v*i-v+1:v*i)];
end
X = X';
obj = x'*Q*x;
% d(i) = +1 (pickup), -1 (deliver)
d = [0;ones(n,1);-ones(n,1)];


% constraint: permutation matrix
for i = 1:v,
   constr = [constr; sum(X(i,:))==1];
   constr = [constr; sum(X(:,i))==1];
end

% constraint: capacity
L = tril(ones(v));
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
ops = sdpsettings('verbose',1);
optimize(constr,obj,ops);
obj_IP = value(obj);
X_IP = value(X);
% tour = round(X_IP*[1:v]');
% tour = [1;tour]

save IP X_IP obj_IP
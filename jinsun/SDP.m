clear all;close all;clc

n = 7 ;             % number of custumers(n)
k = 2;              % capacity    
rng('default');     % random seed: shuffle
alph = 1;
% generate random vehicle, customer pickup and delivery locations from [0,1]x[0,1]
% vehicle node: 1
% pick-up nodes: 2,...,n+1
% delivery nodes: n+2,...,2n+1
vert = rand(2*n+1,2);
load('map.mat');
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

Y = sdpvar(v^2,v^2);
x = diag(Y);  % imagin it as x = [x11 x21 ... xn1 x12 x22 ... xn2 ... x1n ... xnn]';
X = reshape(x,v,v)';

% sdp relaxation constraint
constr = [[1 x'; x Y]>=0];

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



% obj = x'*Q*x;
obj = trace(Y'*Q);
ops = sdpsettings('verbose',1,'solver','sdpt3');
optimize(constr,obj,ops);

disp('========================================')
obj_SDR = value(obj)

% extraction no.1: eigenvalue
y = value(Y);
[eigv,eige] = eig(y);
y_approx = eige(end,end)*eigv(:,end)*eigv(:,end)';
xx = diag(y_approx);
XX = reshape(xx,v,v);
XX(:,1) = [0*ones(v-1,1);1];
XX(end,:) = [1,0*ones(1,v-1)];
[ XX ] = construct_result( XX,v,n,k);
disp('========================================')
xx = reshape(XX,v^2,1);
obj_no1 = xx'*Q*xx
XX_no1 = XX;
tour = round(XX_no1*[1:v]');
tour = [1;tour];

save SDP obj_no1 XX_no1 tour
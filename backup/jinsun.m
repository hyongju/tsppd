% jinsun's version
% TSPPD (1-vehicle, 1-depot)
clear all;close all;clc

n = 5 ;             % number of custumers(n)
k = 2;              % capacity    
rng('shuffle');     % random seed: shuffle
alph = 1;
% generate random vehicle, customer pickup and delivery locations from [0,1]x[0,1]
% vehicle node: 1
% pick-up nodes: 2,...,n+1
% delivery nodes: n+2,...,2n+1
vert = rand(2*n+1,2);
load('matlab.mat');
v = 2*n + 1;        % number of vertices |V| = 2n + 1

constr = [];
% x = binvar(v^2,1);
x = sdpvar(v^2,1);
constr = [constr; 0<=x<=1];
X = [];
for i = 1:v,
    X = [X,x(v*i-v+1:v*i)];
end
X = X';

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

ops = sdpsettings('verbose',1,'solver','mosek');
% ops = sdpsettings('verbose',1);
optimize(constr,obj,ops)
solution = value(X)
tour = round(solution*[1:v]');
tour = [1;tour];




%% draw graph, and optimal tour
figure,
plot(vert(:,1),vert(:,2),'ok','MarkerSize',10,'LineWidth',2); hold on;

for i = 1:length(tour)-1
    dif = vert(tour(i+1),:)-vert(tour(i),:);
    quiver(vert(tour(i),1),vert(tour(i),2),0.1*dif(1)/norm(dif),0.1 *dif(2)/norm(dif),0, 'MaxHeadSize', 1/norm(dif),'LineWidth',2);hold on;
    line([vert(tour(i),1) vert(tour(i+1),1)],[vert(tour(i),2) vert(tour(i+1),2)],'Color','black');hold on;
end

for i = 1:size(vert,1)
    if (mod(i-1,n)) == 0 && i ~= 1
        prtVal = n;
    else
        prtVal = mod(i-1,n);
    end
    if i == 1
        str = 'I';
        col = 'black';
    elseif i <= n+1
        str = 'P';
        col = 'red';
    else
        str = 'D';
        col = 'blue';
    end
    text(vert(i,1)+0.02,vert(i,2)+0.02, sprintf('%.0f %s',prtVal,str),'Color',sprintf('%s',col),'FontSize',16);
end

axis('equal');
axis([0 1 0 1]);
set(gca,'FontSize',16);

% % plot constraint#1
% figure,
% plot(1:v,Acap*value(x),'bo-','MarkerSize',10,'LineWidth',2);
% hold on;
% line([0 v],[k,k],'Color','r','LineWidth',2);
% set(gca,'FontSize',16);
% xlabel('sequence');
% ylabel('# of customers in the vehicle')
% 
% % plot constraint#2
% figure,
% out = Apre*value(x);
% 
% plot(1:n,out,'s-b','MarkerSize',10,'LineWidth',2);
% line([0 n],[0,0],'Color','r','LineWidth',2);
% set(gca,'FontSize',16);
% xlabel('custumer ID');
% ylabel('precedence')
% 
% 









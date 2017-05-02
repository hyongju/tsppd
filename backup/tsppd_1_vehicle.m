% TSPPD (1-vehicle, 1-depot)
clear all;close all;clc

n = 7;             % number of custumers(n)
k = 3;              % capacity         
rng('shuffle');     % random seed: shuffle
% generate random vehicle, customer pickup and delivery locations from [0,1]x[0,1]
% vehicle node: 1
% pick-up nodes: 2,...,n+1
% delivery nodes: n+2,...,2n+1
vert = rand(2*n+1,2);
v = 2*n + 1;        % number of vertices |V| = 2n + 1


% cost matrix (c)
c = zeros(v,v);
for i = 1:size(c,1)
    for j = 1:size(c,2)
        c(i,j) = norm(vert(i,:)-vert(j,:));
    end
end
c = c + eye(v); % add identity to prevent zeros from appearing in the diagonal
% adjacency matrix for a default tour: 1->2->3->..->2n->2n+1
A0 = eye(2*n+1);
A0 = [A0(2:2*n+1,:);A0(1,:)] + eye(v); % added identity to prevent zeros from appearing in the diagonal

% constraint #0 (permutation): row, column sums are 1s
Aperm = ones(v,1);
bperm = ones(v,1);

% constraint #1 (capacity)
% d(i) = +1 (pickup), -1 (deliver)
d = [0;ones(n,1);-ones(n,1)];
L = tril(ones(v,v));

% constraint #2 (precedence), eini = ei - en+1
for i = 1:n
    eini(:,i) = zeros(v,1);
    eini(i+1,i) = 1;
    eini(i+1+n,i) = -1;
end

% integer programming
xVec = sdpvar(v*v,1);           
x = reshape(xVec,v,v);          % decision variables as a v x v matrix
F = [x*Aperm == bperm, x'*Aperm == bperm, x(v,1)== 1 ...,   % constraint #0
    L*x*d <= k*ones(v,1) ...,                               % constraint #1    
    (1:v)*x*eini <= zeros(1,n)];                            % constraint #2

% objective
obj = sum(sum(c.*x*A0*x'));
ops = sdpsettings('verbose',1,'solver','quadprog');
optimize(F,obj,ops)

% generate the optimal tour from x
solution = reshape(value(x),[size(c,1),size(c,2)]);
tour = round(solution*(1:v)');
tour = [1;tour];


% draw graph, and optimal tour
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

% plot constraint#1
figure,
plot(1:v,tril(ones(v,v))*value(x)*d,'bo-','MarkerSize',10,'LineWidth',2);
hold on;
line([0 v],[k,k],'Color','r','LineWidth',2);
set(gca,'FontSize',16);
xlabel('sequence');
ylabel('# of customers in the vehicle')

% plot constraint#2
figure,
for i = 1:size(eini,2)
    out(i) = (1:v)*value(x)*eini(:,i);
end
plot(1:n,out,'s-b','MarkerSize',10,'LineWidth',2);
line([0 n],[0,0],'Color','r','LineWidth',2);
set(gca,'FontSize',16);
xlabel('custumer ID');
ylabel('precedence')


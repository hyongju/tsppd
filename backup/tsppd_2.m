% TSPPD (1-vehicle, 1-depot)
clear all;close all;clc

n = 5;             % number of custumers(n)
k = 3;              % capacity    
rng('shuffle');     % random seed: shuffle
alph = 1;
% generate random vehicle, customer pickup and delivery locations from [0,1]x[0,1]
% vehicle node: 1
% pick-up nodes: 2,...,n+1
% delivery nodes: n+2,...,2n+1
vert = rand(2*n+1,2);
% load('matlab.mat');
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
A0 = [A0(2:2*n+1,:);A0(1,:)] + alph*eye(v); % added identity to prevent zeros from appearing in the diagonal


% d(i) = +1 (pickup), -1 (deliver)
d = [0;ones(n,1);-ones(n,1)];

% constraint #1: permutation

% row sums are 1s
for l = 1:2*n+1
    row1{l} = zeros(size(c));
    for i = 1:v
        for j = 1:v
            if j == l
                row1{l}(i,j) = 1;
            end          
        end
    end
end


for l = 1:length(row1)
    row2{l} = reshape(row1{l}',1,size(c,2)*size(c,1));
    
end

% column sums are 1s
for l = 1:2*n+1
    col1{l} = zeros(size(c));
    for i = 1:v
        for j = 1:v
            if i == l
                col1{l}(i,j) = 1;
            end
        end
    end
end

for l = 1:length(col1)
    col2{l} = reshape(col1{l}',1,size(c,2)*size(c,1));
end

% equality constraints... 
l = 0;
for i = 1:length(row2)
    l = l+1;
    Aperm(l,:) = row2{i};
end
for j = 1:length(col2)
    l = l+1;
    Aperm(l,:) = col2{j};
end
bperm = ones(size(Aperm,1),1);

% generate d's

for l = 1:length(d)
    for i = 1:size(c,1)
%         dRC = [dRC;d];
        for j = 1:size(c,2)
            if i <=l
                Acap(l,i+(j-1)*size(c,1)) = d(j);
            else
                Acap(l,i+(j-1)*size(c,1)) = 0;
            end
        end
    end
end

bcap = ones(size(Acap,1),1)*k;

% tour needs to end at node 1...

Ate = zeros(v,v);
for i = 1:v
    for j = 1:v
        if i == v && j ==1
            Ate(i,j) = 1;
        end
    end
end
Atend = reshape(Ate,size(c,2)*size(c,1),1)';

% constraint #2 (precedence), eini = ei - en+1
for i = 1:n
    eini(i,:) = zeros(1,v);
    eini(i,i+1) = 1;
    eini(i,i+n+1) = -1;
end

N = 1:1:2*n+1;
% N = [N(1,2:2*n+1),1];
Apre = zeros(size(eini,1),v*v);
ApreN = repmat(N,1,v);
% for l = 1:size(eini,1)
%     m = 0;
%     for i = 1:size(eini,2)
%         for j =1:length(N)
%             m = m+1;
%             Apre(l,m) = eini(l,i)*N(1,j);
%         end
% 
%     end
% end
einit = [];
for i = 1:v
    einil = repmat(eini(:,i),1,v);
    einit = [einit einil];
end
Apre = einit.*repmat(ApreN,size(eini,1),1);
bpre = zeros(size(Apre,1),1);

%% interger programming....
Aeq = [Aperm];
beq = [bperm];

bcapzero = zeros(size(bcap));

% Aieq = [Acap;Apre];
% bieq = [bcap;bpre];

x = binvar(v^2,1);
X = reshape(x,v,v)';
xx = reshape(X,v^2,1);

% F = [Acap*xx <=bcap, Apre * xx <= bpre, X * ones(v,1)==ones(v,1),X'*ones(v,1) == ones(v,1),X(v,1) == 1];
F = [];
F = [F,tril(ones(v,v))*X*d <= k];
F = [F, X * ones(v,1)==ones(v,1),X'*ones(v,1) == ones(v,1),X(v,1) == 1];
N = 1:1:v;
for i = 1:n
    eini = zeros(1,v);
    eini(i+1) = 1;
    eini(i+n+1) = -1;
    F = [F,N*X*eini' <= 0];
end
% F = [Acap*x <=bcap, Apre * x < bpre, Aeq*x == beq];
% F = [];

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



obj =x'*Q*x;
ops = sdpsettings('verbose',1);
optimize(F,obj,ops)

% generate the optimal tour from x
solution = reshape(value(xx),[size(c,1),size(c,2)]);
tour = round(solution*[1:v]');
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
plot(1:v,Acap*value(xx),'bo-','MarkerSize',10,'LineWidth',2);
hold on;
line([0 v],[k,k],'Color','r','LineWidth',2);
set(gca,'FontSize',16);
xlabel('sequence');
ylabel('# of customers in the vehicle')

% plot constraint#2
figure,
out = Apre*value(xx);

plot(1:n,out,'s-b','MarkerSize',10,'LineWidth',2);
line([0 n],[0,0],'Color','r','LineWidth',2);
set(gca,'FontSize',16);
xlabel('custumer ID');
ylabel('precedence')
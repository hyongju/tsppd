function [ para,A,b,Base1,Base2,totaldim,P,e1 ] = pre_process_new( n,k,vert,c,alpha)
%pre_process function (for single vehihicle problem), and it helps to
%calculate useful parameters and matrices in the linear constraints later.
%   input: n    --  number of custumers.
%          k    --  capacity of the vehicle.
%          vert --  a (2n+1)-by-2 matrix that stores the map information,
%                   i.e. pick-up, deliver location, and the defult vehicle
%                   location.
%          c    --  map weight(square matrix). Defult: distance between map
%                   nodes.
%          alpha--  auxillary variable to make the Q in 'quadcost'. Defult
%                   value is 1.
%   output: para   --  a structure that stored useful parameters including
%                      para.n, para.k, para.vert, para.c, para.v, para.Q,
%                      para.L, para.d, para.N, para.E, para.dim.
%           A,B    --  Set of cells that represent constraints A{i}*x<=b(i)
%                      Will be used to define the rank-1 constrained manifold

%%%%%%%%%%%%%%%%%%%%%%  PART ONE: para  %%%%%%%%%%%%%%%%%%%%%%%%%
    para.n = n;
    para.k = k;
    para.vert = vert;
    para.v = 2*n+1;   % number of vertices in the map   
    para.L = tril(ones(para.v));
    para.d = [0,ones(1,n),-ones(1,n)]';
    para.N = 1:para.v;  
    para.E = [zeros(1,n);eye(n);-eye(n)];
    para.dim = para.v^2+para.v+para.n+1; % real dimension of 'unknown' variables of the problem
    
    if ~exist('alpha', 'var') || isempty(alpha)
        alpha=1;
    end
    v = para.v;
    
    % construct weight matrix:
    if ~exist('c', 'var') || isempty(c) % if no weight given, use defult
        c = zeros(v,v);
        for i = 1:size(c,1)
            for j = 1:size(c,2)
                c(i,j) = norm(vert(i,:)-vert(j,:)); % calculate distance between two nodes
            end
        end
        
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
    end
    para.c = c;
    
    % construct Q matrix in quadratic cost function
    Q = zeros(v^2,v^2);
    for i = 1:v
        Q((i-1) * v + 1:(i-1) * v + v,(i-1) * v + 1:(i-1) * v + v) = alpha*(c+c')/2;
        if i<v,
            Q(i*v+1:i*v+v,(i-1)*v+1:(i-1)*v+v) = c'/2;
            Q((i-1)*v+1:(i-1)*v+v,i*v+1:i*v+v) = c/2;
        end
    end
    Q((v-1)*v + 1:(v-1)*v + v,1:v) = c'/2;
    Q(1:v,(v-1)*v + 1:(v-1)*v + v) = c/2;    
    para.Q = Q;

%%%%%%%%%%%%%%%%%%%%%%  END OF PART ONE  %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%  PART TWO: constr %%%%%%%%%%%%%%%%%%%%%%%%%
% This part construct Matrix A and b in which A*vec(uu')=b, where u belongs
% to the symfixedrankYYfactory (symmetric matrix manifold) with rank 1.

xdim = v^2;
mudim = (v-0*para.k*2)+n; % mu is the slack variable for capability and precedence constraints.
totaldim = xdim+mudim+1;
para.aux = zeros(totaldim,v^2);
para.aux(2:1+v^2,:) = eye(v^2);
list = zeros(totaldim);
for i = 1:totaldim
    list(:,i) = (i-1)*totaldim+[1:totaldim]';
end

xlist = (2:v^2+1)';
diaglist = diag(list(2:v^2+1,2:v^2+1));
mu2list = diag(list(v^2+2:end,v^2+2:end));
permlist = reshape(xlist,v,v);

A = [];
b = [];
constrcont = 1;
mu2cont = 1;

% 1st: formation of Z* = [1 diag(Y)' ; diag(Y) Y] where Y = xx', and x is
% the vectorization of our target permutation matrix.
A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,1) = 1; % for the 1 in Z*
b(constrcont,1) = 1; constrcont = constrcont+1;

for i = 1:length(xlist), % for diag(Y)=diag(Y)
    A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,xlist(i)) = 1;A(constrcont,diaglist(i)) = -1;
    b(constrcont,1) = 0; constrcont = constrcont+1;
end

% 2nd: permutation constraints --> row and column sums being 1
for i = 1:v,
   if i==v, % we can drop one column sum constraint
       A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,permlist(i,:)) = ones(1,v);
       b(constrcont,1) = 1; constrcont = constrcont+1;
       break
   end
   A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,permlist(i,:)) = ones(1,v);
   b(constrcont,1) = 1; constrcont = constrcont+1;
   A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,permlist(:,i)) = ones(1,v);
   b(constrcont,1) = 1; constrcont = constrcont+1;
end

% 3rd: capacity constraints
A3 = zeros(v,totaldim^2);
for i = 1:v,
    A3(i,2:1+i*length(para.d)) = repmat(para.d',1,i);
end
% A3([1:k,size(A3,1)-k+1:size(A3,1)],:)=[];
temp = size(A3,1);
for i = 1:temp,
   A3(i,mu2list(mu2cont)) = 1;
   mu2cont = mu2cont+1;
end
A(constrcont:constrcont+temp-1,:) = A3;
b(constrcont:constrcont+temp-1,1)=k*ones(temp,1);
constrcont = constrcont + temp;


% 4th: precedence constraints
A4 = zeros(n,totaldim^2);
for i = 1:v,
    A4(:,2+(i-1)*(2*n+1):1+i*(2*n+1)) = [zeros(n,1),i*eye(n),-i*eye(n)];
end
for i = 1:n,
    A4(i,mu2list(mu2cont)) = 1;
    mu2cont = mu2cont+1;
end
A(constrcont:constrcont+n-1,:) = A4;
b(constrcont:constrcont+n-1,1)=zeros(n,1);
constrcont = constrcont + n;

% 5th: return to the 1st point
A(constrcont,:) = zeros(1,totaldim^2); A(constrcont,1+v^2-v+1) = 1;
b(constrcont,1) = 1;


%%%%%%%%%%%%%%%%%%%%%%  END OF PART TWO  %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%  PART THREE: Base %%%%%%%%%%%%%%%%%%%%%%%%%
% This part construct Base1 and Base2 which help to vectorize Z*.
% This part also construct P which is used when calculating the
% gradient and hessian. 
II = eye(totaldim);
P = 0;
e1 = zeros(totaldim,1);
e1(1) = 1;
for i = 1:totaldim,
    Base1{i} = mat_col_switch(II,1,i)*e1; 
    Base2{i} = [zeros(totaldim*(i-1),totaldim);eye(totaldim);zeros(totaldim*(totaldim-i),totaldim)];
    P = P+Base2{i}'*A'*b*Base1{i}';
end
P = P+P';

%%%%%%%%%%%%%%%%%%%%%%  END OF PART THREE  %%%%%%%%%%%%%%%%%%%%%%%


end

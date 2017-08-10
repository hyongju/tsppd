function [ para,A,B ] = pre_process( n,k,vert,c,alpha)
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
% There are five linear constraints for the single vehicle case, and all of
% the linear constraints can be represented as x'*A*x=b somehow. Note that
% we assume x has the form as [1, xx', mu'] where xx represents the
% vectorized permutation matrix and mu is a vector of slack variables. 


% 0th: xx(i) is in {0,1}
% A = cell(1,para.v^2+3*para.v+para.n);
% B = zeros(1,para.v^2+3*para.v+para.n);
constrcont = 0;
% for i = 1:para.v^2,
%    temp = zeros(para.dim,para.dim);
%    temp(1,i+1) = -0.5;
%    temp(i+1,1) = -0.5;
%    temp(i+1,i+1) = 1;
%    A{i} = temp;
%    B(i) = 0;
% end
% constrcont = i;

% 1st & 2nd: the solution we are aiming for is a permutation matrix X, thus
% we need constraints: ones(1,v)*X=ones(1,v) and X*ones(v,1)=ones(v,1).
% Basically these two constraints say there should be only one element in 
% each row and column that is exact 1, since permutation matrix should only
% contain 0&1. 
% For the SDP relaxation, we replace the integer constraint by introducing
% Y=xx', then the permutation matrix X is relaxed to X_tilde=x_tilde*x_tilde'= 
% trans(reshape(diag(Y),para.v,para.v)). Notice that for x={0,1}, x^2=x.
% Thus all the linear constraints can be rewrite asfunctions of diag(Y),
% and this is also what we will do in the following part.
% Also, we assume all the linear constraints read: A{i}*x<=b{i} or A{i}*x==b{i}.
para.m = v+n;
constr.eq.A{1} = [zeros(v,1),repmat(eye(v),1,v),zeros(v,para.m)]; % constraint for row sum
for i = constrcont+1:constrcont+para.v,
    A{i} = diag(constr.eq.A{1}(i-constrcont,:));
    B(i) = 1;
end
constrcont = i;
% constr.eq.B{1} = ones(v,1);

temp = [];
for i = 1:v,
    temp = blkdiag(temp,ones(1,v)); 
end
constr.eq.A{2} = [zeros(v,1),temp,zeros(v,para.m)]; % constraint for column sum
for i = constrcont+1:constrcont+para.v-1,
    A{i} = diag(constr.eq.A{2}(i-constrcont,:));
    B(i) = 1;
end
% i = constrcont+1;
% A{i} = diag(constr.eq.A{2}(1,:));
% B(i) = 1;
constrcont = i;
% constr.eq.B{2} = ones(v,1);

% 3rd: the 3rd one is for capacity constraint, which says that at any time,
% there should not be more than k custumers in the vehicle. The original
% constraint reads: 0<=L*X*d<=k*ones(v,1)
temp = zeros(v,v^2);
for i = 1:v,
   temp(i,1:i*v) = repmat(para.d',1,i); 
end
constr.ineq.A{1} = [zeros(v,1),temp,eye(v),zeros(v,n)];
for i = constrcont+1:constrcont+para.v,
    A{i} = diag(constr.ineq.A{1}(i-constrcont,:));
    B(i) = k;
end
constrcont = i;
% constr.ineq.B{1} = k*ones(v,1);


% 4th: the 4th one is for precedence constraint, which says for every
% custumer, the vehicle should pick-up the custumer before deliver
% him(her). The original contraint reads: N*X*E<=ones(1,n)
temp = [];
for i = 1:v,
    temp = [temp,zeros(n,1),i*eye(n),-i*eye(n)];
end
constr.ineq.A{2} = [zeros(n,1),temp,zeros(n,v),eye(n)];
for i = constrcont+1:constrcont+para.n,
    A{i} = diag(constr.ineq.A{2}(i-constrcont,:));
    B(i) = 0;
end
constrcont = i;

% [r,c] = size(A{i});
% A{i+1} = rand(r,c);
% for ii = 1:r,
%    A{i+1}(ii,ii) = 0;
% end
% constr.ineq.B{2} = -1*ones(n,1);
%%%%%%%%%%%%%%%%%%%%%%  END OF PART TWO  %%%%%%%%%%%%%%%%%%%%%%%%%

end

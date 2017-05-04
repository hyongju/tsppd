function [ para,constr ] = pre_process( n,k,vert,c,alpha)
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
%                      para.L, para.d, para.N, para.E.
%           constr --  a structure that stored matrices for linear
%                      constraints. It contains constr.Z.A{i} and constr.Z.b{i},
%                      where i=1,2,3,4, and Z in {'eq','ineq'}
%                      (See code below for detailed discription for what
%                      those constraints are.)

%%%%%%%%%%%%%%%%%%%%%%  PART ONE: para  %%%%%%%%%%%%%%%%%%%%%%%%%
    para.n = n;
    para.k = k;
    para.vert = vert;
    para.v = 2*n+1;   % number of vertices in the map   
    para.L = tril(ones(para.v));
    para.d = [0,ones(1,n),-ones(1,n)]';
    para.N = 1:para.v;  
    para.E = [zeros(1,n);eye(n);-eye(n)];
    
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
% There are four linear constraints for the single vehicle case.

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
constr.eq.A{1} = repmat(eye(v),1,v); % constraint for row sum
constr.eq.B{1} = ones(v,1);

temp = [];
for i = 1:v,
    temp = blkdiag(temp,ones(1,v)); 
end
constr.eq.A{2} = temp; % constraint for column sum
constr.eq.B{2} = ones(v,1);

% 3rd: the 3rd one is for capacity constraint, which says that at any time,
% there should not be more than k custumers in the vehicle. The original
% constraint reads: L*X*d<=k*ones(v,1)
temp = zeros(v,v^2);
for i = 1:v,
   temp(1,1:i*v) = repmat(para.d',1,i); 
end
constr.ineq.A{1} = temp;
constr.ineq.B{1} = k*ones(v,1);

% 4th: the 4th one is for precedence constraint, which says for every
% custumer, the vehicle should pick-up the custumer before deliver
% him(her). The original contraint reads: N*X*E<=ones(1,n)
temp = [];
for i = 1:v,
    temp = [temp,zeros(n,1),i*eye(n),-i*eye(n)];
end
constr.ineq.A{2} = temp;
constr.ineq.B{2} = zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%%  END OF PART TWO  %%%%%%%%%%%%%%%%%%%%%%%%%
end

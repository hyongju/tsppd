function [ ] = generateILP( c,v,k,ops,n,sample_num )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


% EXTAND c first: add enging point that equals to the starting point
temp = [0,c(1,2:end),c(1,1)];
c = blkdiag(c,c(1,1));
c(end,:) = temp;
c(:,end) = temp';

% generate decision variable x that in A. Here # of x is larger than # of
% elements in A, but those extra elements are set to be 0s.
x = binvar(v+1,v+1,'full');
x(v+1,:) = zeros(1,v+1);
x(:,1) = zeros(v+1,1);
for i = 1:v+1,
    x(i,i) = 0;
end
% set constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = sdpvar(1);                      %
constr = [ y == c(1,1)*v];         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:v,
    constr = [constr; sum(x(:,i+1)) == 1]; % constraints w.r.t. (2)
    constr = [constr; sum(x(i,:)) == 1];   % constraints w.r.t. (3)
end

% capability constraint
Q = sdpvar(v+1,1);   % decision variable for load of vehicle (really need intvar?)
q = [d;0];
for i = 1:v+1,
    for j = 1:v+1,
        if ~isnumeric(x(i,j)),
            constr = [constr; implies(x(i,j),Q(j) == Q(i) + q(j))];
        end
    end
    if i~=v+1,
        constr = [constr; max(0,q(i)) <= Q(i) <= min(k,k+q(i))];
    end
end

% precedence constraint
B = sdpvar(v+1,1);B(1)=0;   % decision variable for begining of service at each vertex
dd = 0.1*ones(size(B));
% construct map of t_ij (travel time from vertex i to j)
t = c; % here we assume the vehicle always travels with constant speed 1
for i = 1:v+1,
    if 2<=i && i<=n+1,
        constr = [constr;B(i)<=B(n+i)];
    end
    for j = 1:v+1,
        if ~isnumeric(x(i,j)),
           constr = [constr;implies(x(i,j),B(j)>=B(i)+t(i,j)+dd(i))]; 
        end
    end
end


% construct object function
obj = sum(sum(c.*x)) + y;

% run ILP with 5 solvers
if n <= 10,           
    ops.solver = 'cplex';
    [txtfile,matfile] = generateFinename(n,sample_num,1,1);
    diary(txtfile);diary on
    ILP = optimize(constr,obj,ops);
    diary off
    save(matfile,'ILP','n','sample_num','k');
end
if n <= 12,
    ops.solver = 'gurobi';
    [txtfile,matfile] = generateFinename(n,sample_num,1,2);
    diary(txtfile);diary on
    ILP = optimize(constr,obj,ops);
    diary off
    save(matfile,'ILP','n','sample_num','k');
end
if n <= 12,
    ops.solver = 'mosek';
    [txtfile,matfile] = generateFinename(n,sample_num,1,3);
    diary(txtfile);diary on
    ILP = optimize(constr,obj,ops);
    diary off
    save(matfile,'ILP','n','sample_num','k');
end
if n <= 12,
    ops.solver = 'baron';
    [txtfile,matfile] = generateFinename(n,sample_num,1,4);
    diary(txtfile);diary on
    ILP = optimize(constr,obj,ops);
    diary off
    save(matfile,'ILP','n','sample_num','k');
end
if n <= 12,
    ops.solver = 'xpress';
    [txtfile,matfile] = generateFinename(n,sample_num,1,5);
    diary(txtfile);diary on
    ILP = optimize(constr,obj,ops);
    diary off
    save(matfile,'ILP','n','sample_num','k');
end   




end


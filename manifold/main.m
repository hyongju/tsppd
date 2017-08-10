clear all,close all,clc; % junk
% profile on
n = 1 ;             % number of custumers(n)
k = 1 ;              % capacity   
v = 2*n+1;

rng('shuffle');     % random seed: shuffle
vert = rand(2*n+1,2);

[ para,A,B ] = pre_process( n,k,vert);

AA = sym(zeros(length(A)));
y = round(rand(1+v^2+v+n,1));
eta = rand(1+v^2+v+n,1);
y = [1  0 1 0 0 0 1 1 0 0   0     1     1     1]';% n=k=1
% y = [1;0;1;0;1;0;0;0;0;1;0;0;1;2];
% y = [1 1 0 0 0 1 0 0 0 1 1 0 1 1]';
% y = [1;0;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;1;1;0;0;0;0;1;0;1;2;2;2;2];% n=k=2

%%%%%%%%%%%%%%%%%%%%%%%%%%%    symbolic
x = sym('x',[para.v^2,1],'real');
mu = sym('mu',[para.m,1],'real');
% y = [1;x;mu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(A)
    for j = 1:length(A)
        i;j;
        AA(i,j) = y'*A{i}*A{j}*y;
    end
    b(i) = eta'*A{i}*y;
end

% e = eig(AA);
% e(end)
rank(AA);

% 
% AA(1,:) = [];
% AA(:,1) = [];
rank(AA)
temp = rank(AA(2:end,2:end));
for i = 2:size(AA,1)-1,
    z = [AA(1:i-1,1:i-1),AA(1:i-1,i+1:end);AA(i+1:end,1:i-1),AA(i+1:end,i+1:end)];
%     if i == 8
%        simplify(AA) 
%     end
    temp = [temp,rank(z)];
end
temp = [temp,rank(AA(1:end-1,1:end-1))];
zzz = find(temp == rank(AA)-1)
% 
% temp = [];
% AA([10,11,12,16],:)=[]
% for i = 1:size(AA,1),
%    z = AA;
%    z(i,:)=[];
%    temp = [temp,rank(z)];
% end
% temp


% profile viewer
% profile off

[alpha,BL] = SolveAlpha(double(AA),b)
double(diag(AA))




function [res,BaseList ]= SolveAlpha(A,b)
    echo = rref(A);
    BaseList = [];
    [r,~] = size(A);
    for i = 1:r,
        BaseList = [BaseList,find(echo(i,:)~=0,1)];
    end
    unBaseList = setdiff(1:r,BaseList);
    Anew = A(BaseList,BaseList);
    bnew = b(BaseList);
    res = zeros(r,1);
    x = Anew\bnew';
    res(BaseList) = x;
    
    for i = 1:r,
       if sum(abs(echo(i,:)))~=0 && length(find(echo(i,:)~=0))>1,
          idx =  find(echo(i,:)~=0,1);
          temp = echo(i,:);
          temp = temp/sum(temp);
          temp = res(idx)*temp;
          res(idx) = 0;
          res = res+temp';
       end
    end
    

end


clear all;close all;clc
n = 3;                  % number of customers: (n-1)/2
C = sym('C',[n n]);
X = sym('X',[n n]);
A0 = eye(n);
alph = 2;
A0 = [A0(2:n,:);A0(1,:)]+alph*eye(n);
A0=sym(A0);
Q = expand(sum(sum(C.*(X*A0*X.'))));

var11 = repmat((1:n)',n,1);
var12 = repmat((1:n)',1,n);
var12 = reshape(repmat((1:n)',1,n)',n*n,1);
var1 = [var11 var12];

for i= 1:length(var1)
    for j = 1:length(var1)
        N(i,j) = jacobian(jacobian(Q,X(var1(i,1),var1(i,2))),X(var1(j,1),var1(j,2)));
        if i == j
            N(i,j) = N(i,j)/2;
        elseif length(coeffs(N(i,j)))==2
            N(i,j) = N(i,j)/2;
        elseif length(coeffs(N(i,j)))==1 && coeffs(N(i,j)) == 2
            N(i,j) = N(i,j)/2;
        end
    end
end

for i = 1:10
    rng('shuffle');
    rndsd = rand(n);
    symMat = rndsd'*rndsd;          % artificially created positive definite matrix...
    M = subs(N,C,symMat);
    dM = double(M);
    psd(i) = all(round(eig(dM),5) >= 0);    % check if the matrix is PSD...
end

length(find(psd)) % ? out of 10 is PSD...
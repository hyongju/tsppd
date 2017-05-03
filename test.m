clear all,close all,clc
% test
% B = [1 2 1;2 1 3;1 3 2];
B = [2 4 2;4 8 4;2 4 2];
x = binvar(3,1);
obj = x'*B*x;
ops = sdpsettings('verbose',1);
optimize([],obj,ops)
value(x)
value(obj)
%%
X = sdpvar(3,3);
diagX = [X(1,1);X(2,2);X(3,3)];
obj = trace(X'*B);
Xbar = [1 diagX';diagX X];
constr = [Xbar>=0];
ops = sdpsettings('verbose',1,'solver','mosek');
optimize([],obj,ops)
value(X)
value(obj)

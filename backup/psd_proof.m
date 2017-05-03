
clear all;close all;clc

n = 11;
Q = zeros(n,n);
Q = sym(Q);
for i = 1:size(Q,1)
    Q(i,i) = sym('c');
end
for i = 1:size(Q,1)-1
    Q(i+1,i) = 0.5*sym('c');
    Q(i,i+1) = 0.5*sym('c');
end
Q(n,1) = 0.5*sym('c');
Q(1,n) = 0.5*sym('c');
Q
for i = 1:n-1
    R = Q(1:(n-i+1)-1,1:(n-i+1)-1) - Q(1:(n-i+1)-1,(n-i+1))*(1/Q((n-i+1),(n-i+1)))*Q((n-i+1),1:(n-i+1)-1)
    clear Q;
    Q = R;
end
function [ f,gradf ] = generate_fun( x,A,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    f = A*x.^2-b;
    gradf = 2*diag(x)*A';
end

